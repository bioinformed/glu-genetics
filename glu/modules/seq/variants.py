from __future__ import division

import os
import sys

from   collections       import defaultdict

import numpy as np
import pysam

from   glu.lib.utils     import Counter, namedtuple, iter_queue
from   glu.lib.fileutils import table_reader, table_writer

from   glu.modules.seq.annotate import VariantAnnotator


Variant = namedtuple('Variant', 'chrom start end depth ref a1 d1 a2 d2 ratio pAA pAB pBB')


REF_GENOME='/CGF/DCEGProjects/Exome/Pilot/build4/reference/hg18.fa'


def option_parser():
  import optparse

  usage = 'usage: %prog [options] infile.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output MUT file')
  return parser


def geno_prob(start,end,loc,ref=None,prior=0.05):
  names,bases,quals = zip(*loc)

  bases = np.array(bases, dtype=object)
  quals = np.array(quals, dtype=float)
  perr  = 10**(-quals/10)

  alleles = Counter(bases).most_common(2)

  if not alleles:
    return None

  a1,c1 = alleles[0]

  if len(alleles)==2:
    a2,c2 = alleles[1]
  elif a1==ref:
    a2,c2 = 'N',0
  else:
    a2,c2 = ref,0

  if a2==ref:
    a1,a2 = a2,a1
    c1,c2 = c2,c1

  if a1==ref and c2<2:
    return None
  if a1!=ref and max(c1,c2)<2:
    return None

  m1 = bases==a1
  m2 = bases==a2

  a1 = a1.replace('-','')
  a2 = a2.replace('-','')

  while ref and a1 and a2 and a1[0]==a2[0]==ref[0]:
    ref    = ref[1:]
    a1     = a1[1:]
    a2     = a2[1:]
    start += 1

  while ref and a1 and a2 and a1[-1]==a2[-1]==ref[-1]:
    ref  = ref[:-1]
    a1   = a1[:-1]
    a2   = a2[:-1]
    end -= 1

  a1  = a1  or '-'
  a2  = a2  or '-'
  ref = ref or '-'

  if a1==ref:
    priorAA = 1-prior
    priorAB =   prior/2
    priorBB =   prior/2
  else:
    priorAA = priorBB = priorAB = 1

  pgood = 1-perr

  pD_AA = np.exp(np.log(np.where(m1,    pgood,   perr  )).sum())*priorAA
  pD_BB = np.exp(np.log(np.where(m2,    pgood,   perr  )).sum())*priorBB
  pD_AB = np.exp(np.log(np.where(m1|m2, pgood/2, perr/2)).sum())*priorAB

  pD = pD_AA+pD_BB+pD_AB

  AA  = a1+a1
  AB  = a1+a2
  BB  = a2+a2

  pAA = pD_AA/pD
  pAB = pD_AB/pD
  pBB = pD_BB/pD

  if ref==a1 and (a1=='-' or a2=='-') and pAA>0.25:
    return None

  if ref==a1 and pAA>0.95:
    return None

  ratio = c2/(c1+c2)*100.

  if 0:
    print
    print 'ref:',ref
    print 'bases:',bases
    print 'allele1:',bases[m1],perr[m1]
    print 'allele2:',bases[m2],perr[m2]
    print 'P(%s)=%.4f P(%s)=%.4f P(%s)=%.4f ratio=%.2f' % (AA,pAA,AB,pAB,BB,pBB,c2/(c1+c2)*100.)

  return start,end,ref,a1,sum(m1),a2,sum(m2),ratio,pAA,pAB,pBB


def pileup_iter(inbam, reference):
  for chrom in inbam.references:
    for p in inbam.pileup(chrom):
      loc = []
      pileups = p.pileups
      for read in pileups:
        align     = read.alignment
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

      if loc:
        yield chrom,p,loc


def extend_variant(pos, loc, ref_seq, pileup):
  start = pos
  end   = pos+1

  prefix = defaultdict(list)
  for name,bases,qual in loc:
    prefix[bases].append(name)

  groups = prefix.values()

  while 1:
    next_chrom,next_p,next_loc = pileup.peek()
    next_pos = next_p.pos

    if next_pos!=pos+1:
      break

    next_loc = dict( (l[0],l) for l in next_loc )
    nucs = []
    for group in groups:
      n = set( next_loc[name][1] for name in group if name in next_loc )
      if len(n) != 1:
        break
      nucs.append(next(iter(n)))

    ref = ref_seq[next_pos]
    if all(n==ref for n in nucs):
      break

    new_loc = []
    for name,bases,qual in loc:
      if name in next_loc:
        bases += next_loc[name][1]
        if len(bases)>1:
          bases = bases.replace('-','') or '-'
      new_loc.append( (name,bases,qual) )

    loc = new_loc
    pos = next_pos
    end = next_pos+1

  return loc,end


def call_variants(inbam, reference):
  last_chrom = None

  pileup = pileup_iter(inbam, reference)
  pileup = iter_queue(pileup)

  for chrom,p,loc in pileup:
    if chrom!=last_chrom:
      last_chrom = chrom
      ref_seq = reference.fetch(chrom).upper()
      print >> sys.stderr, 'Chrom %s: len=%d' % (chrom,len(ref_seq))

    ref_nuc = ref_seq[p.pos]
    seen = Counter(l[1] for l in loc)
    if not seen or (len(seen)==1 and ref_nuc in seen):
      continue

    #print chrom,p.pos # ,p.n,ref_nuc,loc
    loc,end = extend_variant(p.pos, loc, ref_seq, pileup)
    ref_nuc = ref_seq[p.pos:end]
    geno=geno_prob(p.pos,end,loc,ref_nuc)

    if geno:
      start,end,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB = geno
      ref_nuc = ref_seq[start:end]
      #ratio = '%.2f' % ratio
      #pAA   = '%.3g' % pAA
      #pAB   = '%.3g' % pAB
      #pBB   = '%.3g' % pBB
      yield Variant(chrom,start,end,p.n,ref,a1,d1,a2,d2,ratio,pAA,pAB,pBB)


def annotate_variants(variants):
  va = VariantAnnotator()

  for v in variants:
    evidence = []
    if v.ref!=v.a1 and max(v.pAA,v.pAB)>0.05:
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a1))
    if v.ref!=v.a2 and max(v.pAB,v.pBB)>0.05:
      evidence.extend(va.classify(v.chrom, v.start, v.end, v.a2))

    #if 'NON-SYNONYMOUS' not in [ e[6] for e in evidence ]:
    #  continue

    evidence = [ e for e in evidence if e[6]=='NON-SYNONYMOUS' ]
    #if not evidence:
    #  evidence = [[]]

    for e in evidence:
      yield v+tuple(e[3:])


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

  variants = call_variants(inbam, reference)

  if 1:
    header.extend(['INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                   'REF_NUC','VAR_NUC','REF_AA','VAR_AA', 'dbSNP_exact','dbSNP_inexact'])
    variants = annotate_variants(variants)

  out.writerow(header)
  out.writerows(variants)


if __name__ == '__main__':
  main()
