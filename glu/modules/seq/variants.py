# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Simple probabilistic variant caller'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import sys

from   collections              import defaultdict
from   itertools                import groupby
from   operator                 import itemgetter

import numpy as np
import pysam

from   glu.lib.utils            import Counter, namedtuple, iter_queue, unique
from   glu.lib.fileutils        import table_reader, table_writer

from   glu.lib.genedb           import open_genedb
from   glu.modules.seq.annotate import VariantAnnotator


Variant = namedtuple('Variant', 'chrom start end depth ref a1 d1 a2 d2 ratio pAA pAB pBB')


#REF_GENOME='/CGF/DCEGProjects/Exome/Pilot/build4/reference/hg18.fa'
#REF_GENOME='/CGF/DCEGProjects/Exome/Pilot/build4/reference/hg19.fa'
#GENE_DB='/home/jacobske/projects/genedb/out/genedb_hg19_snp131_rtree.db'


def clean_allele(a):
  return a.replace('<','').replace('>','').replace('-','')


def geno_prob(start,end,loc,ref=None,prior=0.05,min_variant=None,collapse=True,alleles=None):
  names,bases,quals = zip(*loc)

  if not collapse:
    bases = [ clean_allele(b) for b in bases ]
  elif ref:
    ref = ref[0]+'<' + '>'.join(ref[1:])
  else:
    ref = '<'

  #print '  counts =',Counter(bases)

  if alleles:
    if len(alleles)>2:
      raise ValueError

    alleles = list(alleles)

    if len(alleles)==1 and ref not in alleles:
      alleles.append(ref)

    if len(alleles)!=2:
      raise ValueError

    valid_bases = [ b for b in bases if b in alleles ]
  else:
    valid_bases = bases
    alleles = ref or 'N','N'

  bases = np.array(bases, dtype=object)
  quals = np.array(quals, dtype=float)

  common_alleles = Counter(valid_bases).most_common(2)

  a1 = a2 = 'N'
  c1 = c2 = 0

  if not common_alleles and alleles:
    a1,a2 = alleles
  else:
    a1,c1 = common_alleles[0]

    if len(common_alleles)==2:
      a2,c2 = common_alleles[1]
    elif a1==alleles[0]:
      a2,c2 = alleles[1],0
    else:
      a2,c2 = alleles[0],0

  c = c1+c2

  if not c:
    return start,end,ref,a1,c1,a2,c2,0,0.,0.,0.

  if not c2 and a1==ref:
    return start,end,ref,a1,c1,a2,c2,0,1.,0.,0.

  ca1  = clean_allele(a1)  or '-'
  ca2  = clean_allele(a2)  or '-'
  cref = clean_allele(ref) or '-'

  if ca2==cref or (ca1!=cref and (len(ca1),ca1)>(len(ca2),ca2)):
    a1,a2 = a2,a1
    c1,c2 = c2,c1

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

  a1  = clean_allele(a1)  or '-'
  a2  = clean_allele(a2)  or '-'
  ref = clean_allele(ref) or '-'

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
        yield chrom,start,end,a1,a2
      except ValueError:
        pass


def parse_locations(filename):
  for row in table_reader(filename,want_header=False):
    if len(row)>=5:
      try:
        chrom,start,end,a1,a2 = row[:5]
        start = int(start)
        end   = int(end)
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

  #if not seen:
  #  return []

  return stripped


def extend_variant(pos, loc, ref_seq, pileup, exclude_start=False, fixed_end=None):
  start = pos
  end   = pos

  #print '  extending from pos=%d to %d' % (pos, fixed_end)

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
  if fixed_end is None and len(prefix)==1 and next(iter(prefix))==ref_seq[pos]:
    return loc,end

  groups = prefix.values()

  while 1:
    #print '  Extend loop'

    if fixed_end is not None and end>=fixed_end:
      #print '  exiting due to fixed end (%d>%d)' % (end,fixed_end)
      break

    try:
      next_chrom,next_pos,next_depth,next_loc = pileup.peek(0)
      #print '  peeking at pile up data at',next_pos

    except ValueError:
      #print '  no more pileup data after',pos
      break

    if next_pos!=pos+1:
      #print '  missing extension point at',next_pos
      break

    next_loc = dict( (l[0],l) for l in next_loc )

    if fixed_end is None:
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
    geno    = geno_prob(start,end,loc,ref_nuc,min_variant=2,collapse=True)

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

    if pos>start:
      continue

    exclude = False

    # Not an indel
    if pos==start-1:
      if len(a1.replace('-','')) == len(a2.replace('-','')) == end-start:
        try:
          chrom,pos,depth,loc = next(pileup)
        except StopIteration:
          depth,loc = 0,[]
          #print 'Nothing found after start-1'
      else:
        # Have to collect adjacent bases
        loc    = strip_one(loc)
        start += 1
        pos   += 1

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

    if 0: # start>=end:
      print '  stripping one more position...?'
      loc     = strip_one(loc)
      end    -= 1
      var_end = end
    else:
      loc,var_end = extend_variant(pos, loc, ref_seq, pileup, exclude_start=exclude, fixed_end=end)

    #print '   extended end to %d: %s' % (end,loc)

    if var_end!=end:
      #print '  inconsistent end location (%d!=%d)' % (var_end,end)
      continue

    ref_nuc = ref_seq[pos:end]
    geno    = geno_prob(pos,end,loc,ref_nuc,collapse=False,alleles=(a1,a2))

    #print loc
    #print
    #print geno

    start,end,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB = geno

    if end<start:
      end = start

    ref_nuc = ref_seq[start:end]
    yield Variant(chrom,start,end,depth,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB)


def annotate_variants(variants, gene_db, reference):
  va = VariantAnnotator(gene_db, reference)

  for v in variants:
    evidence = []
    if v.ref!=v.a2: #  and (v.ref!=v.a1 or v.pAA<0.95):
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a2, False))
    if v.ref!=v.a1 or (v.ref==v.a1==v.a2): #  and (v.ref!=v.a2 or v.pBB<0.95):
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a1, False))

    #if 'NON-SYNONYMOUS' not in [ e[6] for e in evidence ]:
    #  continue
    #evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[6] ]
    #if not evidence:
    #  continue
    if not evidence:
      evidence = [[]]

    for e in evidence:
      yield v+tuple(e[3:])


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('bamfile', help='Input BAM file')

  parser.add_argument('-a', '--annotate', action='store_true',
                    help='Annotate called variants')
  parser.add_argument('-g', '--genedb',   metavar='NAME',
                    help='Genedb genome annotation database name or file')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                    help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--locations', metavar='FILE',
                    help='Locations at which to call variants')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  reference = pysam.Fastafile(options.reference)

  flags = 'rb' if options.bamfile.endswith('.bam') else 'r'
  inbam = pysam.Samfile(options.bamfile, flags)

  header = ['CHROM','START','END','DEPTH','REF_NUC',
            'ALLELE_A','DEPTH_A', 'ALLELE_B','DEPTH_B',
            'RATIO', 'PROB_AA', 'PROB_AB', 'PROB_BB']

  out = table_writer(options.output,hyphen=sys.stdout)

  if options.locations:
    locations = parse_locations(options.locations)
    variants = call_variants_locations(inbam, reference, locations)
  else:
    variants = call_variants_all(inbam, reference)

  if options.annotate:
    header.extend(['INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                   'REF_NUC','VAR_NUC','REF_AA','VAR_AA', 'dbSNP_exact','dbSNP_inexact'])
    variants = annotate_variants(variants, options.genedb, options.reference)

  out.writerow(header)
  out.writerows(variants)


if __name__ == '__main__':
  main()
