from __future__ import division

import os
import sys

from   collections       import defaultdict
from   itertools         import groupby
from   operator          import itemgetter

import numpy as np
import pysam

from   glu.lib.utils     import Counter, namedtuple, iter_queue, unique
from   glu.lib.fileutils import table_reader, table_writer

from   glu.modules.seq.annotate import VariantAnnotator


Variant = namedtuple('Variant', 'chrom start end depth ref a1 d1 a2 d2 ratio pAA pAB pBB')


#REF_GENOME='/CGF/DCEGProjects/Exome/Pilot/build4/reference/hg18.fa'
REF_GENOME='/CGF/DCEGProjects/Exome/Pilot/build4/reference/hg19.fa'
GENE_DB='/home/jacobske/projects/genedb/out/genedb_hg19_snp131_rtree.db'


def geno_prob(start,end,loc,ref=None,prior=0.05,min_variant=None,collapse=True):
  names,bases,quals = zip(*loc)

  bases = np.array(bases, dtype=object)
  quals = np.array(quals, dtype=float)

  alleles = Counter(bases).most_common(2)

  if not alleles:
    return start,end,ref,'',0,'',0,0,0,0,0

  if ref:
    ref = ref[0]+'<' + '>'.join(ref[1:])
  else:
    ref = '<'

  a1,c1 = alleles[0]

  if len(alleles)==2:
    a2,c2 = alleles[1]
  elif a1==ref:
    a2,c2 = 'N<',0
  else:
    a2,c2 = ref,0

  m1 = bases==a1
  m2 = bases==a2

  a1 = a1.replace('-','')
  a2 = a2.replace('-','')

  if collapse:
    #print '  before collapse:',start,end,ref,a1,a2
    while ref and a1 and a2 and a1[0]==a2[0]==ref[0]:
      if ref[0] in '<>':
        start += 1
        #print '  bumping start location to',start
      ref    = ref[1:]
      a1     = a1[1:]
      a2     = a2[1:]

    while ref and a1 and a2 and a1[-1]==a2[-1]==ref[-1]:
      if ref[-1]=='>':
        end -= 1
      ref  = ref[:-1]
      a1   = a1[:-1]
      a2   = a2[:-1]

    #print '  after collapse:',start,end,ref,a1,a2

    a1  =  a1.replace('<','').replace('>','') or '-'
    a2  =  a2.replace('<','').replace('>','') or '-'
    ref = ref.replace('<','').replace('>','') or '-'

  if a2==ref:
    a1,a2 = a2,a1
    c1,c2 = c2,c1
    m1,m2 = m2,m1
  elif a1!=ref and c1==c2 and (len(a1),a1)>(len(a2),a2):
    a1,c1 = a2,c2
    m1,m2 = m2,m1

  ratio = c2/(c1+c2)*100.

  if min_variant is not None:
    if (a1==ref and c2<min_variant) or max(c1,c2)<min_variant:
      return start,end,ref,a1,c1,a2,c2,ratio,0,0,0

  if a1==ref:
    priorAA = 1-prior
    priorAB =   prior/2
    priorBB =   prior/2
  else:
    priorAA = priorBB = priorAB = 1

  perr  = 10**(-quals/10)
  pgood = 1-perr

  log_pD_AA = np.log(np.where(m1,    pgood,   perr  )).sum()+np.log(priorAA)
  log_pD_BB = np.log(np.where(m2,    pgood,   perr  )).sum()+np.log(priorBB)
  log_pD_AB = np.log(np.where(m1|m2, pgood/2, perr/2)).sum()+np.log(priorAB)

  norm = max(log_pD_AA,log_pD_BB,log_pD_AB)

  pD_AA = np.exp(log_pD_AA-norm)
  pD_BB = np.exp(log_pD_BB-norm)
  pD_AB = np.exp(log_pD_AB-norm)

  if not np.isfinite(pD_AA):
    pD_AA = 0
  if not np.isfinite(pD_BB):
    pD_BB = 0
  if not np.isfinite(pD_AB):
    pD_AB = 0

  pD  = pD_AA+pD_BB+pD_AB
  pAA = pD_AA/pD
  pAB = pD_AB/pD
  pBB = pD_BB/pD

  return start,end,ref,a1,c1,a2,c2,ratio,pAA,pAB,pBB


def parse_locations(filename):
  for row in table_reader(filename,want_header=False):
    if len(row)>=9:
      try:
        chrom,start,end,depth,ref_nuc,a1,c1,a2,c2 = row[:9]
        start = int(start)
        end   = int(end)
        c1    = int(c1)
        c2    = int(c2)
        yield chrom,start,end,a1,a2
      except ValueError:
        pass


def seed_variant(p):
  loc = []
  pileups = p.pileups
  for read in pileups:
    align = read.alignment
    if read.is_del or read.indel<0:
      bases     = '-'
      quals     = '+'
    else:
      pos_start = read.qpos
      pos_end   = read.qpos+read.indel+1
      bases     = align.seq[pos_start:pos_end]
      quals     = align.qual[pos_start:pos_end]

    if 'N' in bases:
      continue

    loc.append( (align.qname,bases,ord(quals[0])-33) )

  return loc


def strip_reference(pos, pileup, ref_nuc):
  stripped = []
  for name,bases,qual in pileup:
    if bases and bases[0]==ref_nuc:
      bases = bases[1:]
      stripped.append( (name,bases,qual) )
    else:
      return pileup,pos

  return stripped,pos+1


def start_variant(loc):
  new_loc = []
  for name,bases,qual in loc:
    if bases:
      bases = bases[0]+'<'+bases[1:]
    else:
      bases = '<'
    new_loc.append( (name,bases,qual) )
  return loc


def strip_one(pileup):
  stripped = []
  seen = set()
  for name,bases,qual in pileup:
    if bases:
      bases = bases[1:]

    # FIXME: wrong qual
    # FIXME: does not detect middle of read
    stripped.append( (name,bases,qual) )
    if bases:
      seen.add(bases)

  if not seen:
    return []

  return stripped


def extend_variant(pos, loc, ref_seq, pileup, exclude_start=False):
  start = pos
  end   = pos

  if not exclude_start:
    end += 1
  else:
    pos -= 1

  new_loc = []
  span    = False
  prefix  = defaultdict(list)
  for name,bases,qual in loc:
    prefix[bases].append(name)
    if bases:
      bases = bases[0]+'<'+bases[1:]
    else:
      bases = '<'
    new_loc.append( (name,bases,qual) )
  loc = new_loc

  # Do not extend if only reference bases are present
  if len(prefix)==1 and next(iter(prefix))==ref_seq[pos]:
    return loc,end

  groups = prefix.values()

  while 1:
    try:
      next_chrom,next_pos,next_depth,next_loc = pileup.peek(0)
      #print '  peeking at pile up data at',next_pos

    except ValueError:
      #print '  no more pileup data after',pos
      break

    if next_pos!=pos+1:
      #print '  missing extention point at',next_pos
      break

    next_loc = dict( (l[0],l) for l in next_loc )
    nucs = []
    for group in groups:
      n = set( next_loc[name][1] for name in group if name in next_loc )
      if len(n) != 1:
        #print '  variant and non-variant groups changed at',next_pos
        break
      nucs.append(next(iter(n)))

    ref = ref_seq[next_pos]
    if all(n==ref for n in nucs):
      #print '  all reads are reference at',next_pos
      break

    new_loc = []
    for name,bases,qual in loc:
      if name in next_loc:
        new_bases = next_loc[name][1]
        if new_bases:
          bases += new_bases[0] + '>' + new_bases[1:]
        else:
          bases += '>'
        if len(bases)>1:
          bases = bases.replace('-','') or '-'
      new_loc.append( (name,bases,qual) )

    # Advance pileup and mark new location
    next(pileup,None)
    loc = new_loc
    pos = next_pos
    end = pos+1

  return loc,end


def pileup_iter(inbam, locations=None):
  # Full scan over pileup locations
  references = inbam.references
  for p in inbam.pileup():
    if 0:
      seed = seed_variant(p)
    else:
      seed = p.bases

    if seed:
      yield references[p.tid],p.pos,p.n,seed


def location_iter(inbam, locations):
  # Scan over pre-defined locations
  locations = unique(sorted(locations))

  for chrom,start,end,a1,a2 in locations:
    pileup = inbam.pileup(chrom,start-1,end+1)
    locs   = []

    while 1:
      p = next(pileup,None)

      if p is None:
        break
      elif not start-1<=p.pos<=end+1:
        continue

      seed = p.bases

      if not seed:
        break

      locs.append( (chrom,p.pos,p.n,seed) )

    yield chrom,start,end,a1,a2,locs


def call_variants_all(inbam, reference):
  last_chrom = None

  pileup = pileup_iter(inbam)
  pileup = iter_queue(pileup)

  for chrom,start,depth,loc in pileup:
    if chrom!=last_chrom:
      last_chrom = chrom
      ref_seq = reference.fetch(chrom).upper()
      print >> sys.stderr, 'Chrom %s: len=%d' % (chrom,len(ref_seq))

    #print
    #print chrom,start,depth

    loc,end = extend_variant(start, loc, ref_seq, pileup)

    ref_nuc = ref_seq[start:end]
    geno    = geno_prob(start,end,loc,ref_nuc,min_variant=2)

    start,end,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB = geno

    # Here come the filters
    if not a1 or not a2:
      continue

    if ref==a1:
      if d2<2:
        continue

      if (a1=='-' or a2=='-') and pAA>0.25:
        continue

      if pAA>0.95:
        continue

    elif max(d1,d2)<2:
      continue

    #print '   ext=',loc
    #print '  geno=',geno

    ref_nuc = ref_seq[start:end]
    yield Variant(chrom,start,end,depth,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB)


def call_variants_locations(inbam, reference, locations):
  last_chrom = None

  pileup_locations = location_iter(inbam, locations)

  for chrom,start,end,a1,a2,pileup in pileup_locations:
    #print
    #print '!!! variant %s:%d-%d %s/%s' % (chrom,start,end,a1,a2)

    if not pileup:
      continue

    if chrom!=last_chrom:
      last_chrom = chrom
      ref_seq = reference.fetch(chrom).upper()
      print >> sys.stderr, 'Chrom %s: len=%d' % (chrom,len(ref_seq))

    #print chrom,start,end,a1,a2
    #for p in pileup:
    #  print '  ',p

    pileup = iter_queue(pileup)

    chrom,pos,depth,loc = next(pileup)

    exclude = False
    if pos==start-1:
      loc  = strip_one(loc)
      pos += 1

      if not loc:
        #print '  no variant insertions at start-1'
        try:
          chrom,pos,depth,loc = next(pileup)
        except StopIteration:
          depth,loc = 0,[]
          #print 'Nothing found after start-1'

      else:
        exclude = True
        #print '  using variant insertions at start-1'

    loc,end = extend_variant(pos, loc, ref_seq, pileup, exclude_start=exclude)
    #print '   extended end to %d: %s' % (end,loc)

    #if stripped:
    #  end -= 1

    ref_nuc = ref_seq[pos:end]
    geno    = geno_prob(pos,end,loc,ref_nuc)

    #print geno

    start,end,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB = geno

    ref_nuc = ref_seq[start:end]
    yield Variant(chrom,start,end,depth,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB)


def annotate_variants(variants, gene_db, reference):
  va = VariantAnnotator(gene_db, reference)

  for v in variants:
    evidence = []
    if v.ref!=v.a2 and (v.ref!=v.a1 or v.pAA<0.95):
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a2))
    if v.ref!=v.a1 and (v.ref!=v.a2 or v.pBB<0.95):
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a1))

    #if 'NON-SYNONYMOUS' not in [ e[6] for e in evidence ]:
    #  continue

    #evidence = [ e for e in evidence if e[6]=='NON-SYNONYMOUS' ]
    #if not evidence:
    #  evidence = [[]]

    for e in evidence:
      yield v+tuple(e[3:])


def option_parser():
  import optparse

  usage = 'usage: %prog [options] infile.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--locations', dest='locations', metavar='FILE',
                    help='Locations at which to call variants')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  reference = pysam.Fastafile(REF_GENOME)

  flags = 'rb' if args[0].endswith('.bam') else 'r'
  inbam = pysam.Samfile(args[0], flags)

  header = ['CHROM','START','END','DEPTH','REF_NUC',
            'ALLELE_A','DEPTH_A', 'ALLELE_B','DEPTH_B',
            'RATIO', 'PROB_AA', 'PROB_AB', 'PROB_BB']

  out = table_writer(options.output,hyphen=sys.stdout)

  if options.locations:
    locations = parse_locations(options.locations)
    variants = call_variants_locations(inbam, reference, locations)
  else:
    variants = call_variants_all(inbam, reference)

  if 1:
    header.extend(['INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                   'REF_NUC','VAR_NUC','REF_AA','VAR_AA', 'dbSNP_exact','dbSNP_inexact'])
    variants = annotate_variants(variants, GENE_DB, REF_GENOME)

  out.writerow(header)
  out.writerows(variants)


if __name__ == '__main__':
  main()
