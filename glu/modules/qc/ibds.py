# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute IBS and IBD sharing for pairs of samples'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id: $'


import sys

import numpy as np

from   glu.lib.utils             import pair_generator
from   glu.lib.fileutils         import table_writer
from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.transform import _intersect_options, _union_options
from   glu.lib.genolib.genoarray import genoarray_ibs, genotype_count_matrix

from   glu.modules.qc.dupcheck   import file_pairs, progress_bar


def allele_counts(models,geno_counts):
  balleles = np.zeros( (len(models), 4), dtype=int)
  for i,model in enumerate(models):
    if len(model.alleles)==3:
      ref   = model.alleles[1]
      other = model.alleles[2]
      balleles[i,model[other,ref].index]   = 1
      balleles[i,model[other,other].index] = 2

  n = 2*geno_counts.sum(axis=1)
  x = (geno_counts*balleles).sum(axis=1)

  return x,n


def estimate_ibs_given_ibd(x,n):
  y  = n-x
  p  = x/n
  q  = 1-p

  # Constants to adjust for bias as per Nei M. Estimation of Average
  # Heterozygosity and Genetic Distance from a Small Number of Individuals.
  # PMID: 17248844, PMCID: PMC1213855
  n2 = n/(n-1) * n/(n-2)
  n3 =   n2    * n/(n-3)
  xx = (x-1)/x
  yy = (y-1)/y

  #      2*P(AA)*P(BB)   [bias constants]
  e00 = (2*p*p*q*q   *   xx * yy      * n3 )

  #      2*P(AA)*P(AB)   [bias constants]
  #      2*P(BB)*P(AB)
  e01 = (4*p*p*p*q   *   xx * (x-2)/x * n3
      +  4*p*q*q*q   *   yy * (y-2)/y * n3 )

  #      P(AA)P(AB)/P(A) [bias constants]
  #      P(BB)P(AB)/P(B)
  e11 = (2*p*p*q     *   xx           * n2
      +  2*p*q*q     *   yy           * n2 )

  mask = np.isfinite(e00) & np.isfinite(e01) & np.isfinite(e11)

  e00 = e00[mask].mean()
  e01 = e01[mask].mean()
  e11 = e11[mask].mean()

  e02 = 1-e00-e01
  e12 = 1-e11

  ibs_given_ibd = np.array([[ e00, e01, e02 ],
                            [   0, e11, e12 ],
                            [   0,   0,   1 ]], dtype=float)

  return ibs_given_ibd


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotypes'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True)

  parser.add_option('--frequencies', dest='frequencies', metavar='FILE', action='append',
                    help='Optional genotype file to estimate allele frequencies')
  parser.add_option('--includetest', dest='includetest', metavar='FILE', action='append',
                    help='List of samples to test')
  parser.add_option('--excludetest', dest='excludetest', metavar='FILE', action='append',
                    help='List of samples not to test')
  parser.add_option('--testpairs', dest='testpairs', metavar='FILE',
                    help='File containing a list of pairs to test')

  parser.add_option('-t', '--threshold', dest='threshold', metavar='N', type='float', default=0.90,
                    help='Output only pairs with estimated IBD0 sharing less than N (default=0.90)')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='output table file name')
  parser.add_option('-P', '--progress', dest='progress', action='store_true',
                    help='Show analysis progress bar, if possible')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    sys.exit(2)

  # Allow specification of a different genotype file with which to estimate allele frequencies
  if options.frequencies == args[0]:
    options.frequencies = None

  if options.frequencies:
    freqfile = options.frequencies
  else:
    freqfile = args[0]

  sys.stderr.write('Opening genotype data file...\n')
  genos = load_genostream(freqfile,format=options.informat,genorepr=options.ingenorepr,
                                   genome=options.loci,phenome=options.pedigree,
                                   transform=options, hyphen=sys.stdin)

  sys.stderr.write('Loading genotypes...\n')

  # Materialize only if genotypes for frequency estimation are also the test set
  if not options.frequencies:
    genos = genos.as_sdat().materialize()

  sys.stderr.write('Computing genotype frequencies...\n')
  loci,samples,geno_counts = genotype_count_matrix(genos)

  x,n = allele_counts(genos.models,geno_counts)
  ibs_given_ibd = estimate_ibs_given_ibd(x,n)

  # Load distinct genotype set for tests if needed
  if options.frequencies:
    # Merge in test includes and excludes
    options.includesamples.extend(options.includetest)
    options.excludesamples.extend(options.excludetest)

    genos = load_genostream(args[0],format=options.informat,genorepr=options.ingenorepr,
                                    genome=options.loci,phenome=options.pedigree,
                                    transform=options, hyphen=sys.stdin).as_sdat().materialize()

  # Apply test includes and excludes to the existing genotypes
  elif options.includetest or options.excludetest:
    genos = genos.transformed(include_samples=_intersect_options(options.includetest),
                              exclude_samples=    _union_options(options.excludetest)).materialize()

  if not len(genos):
    sys.stderr.write('No samples found.  Exiting...\n')
    return

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE1','SAMPLE2','COMPARISONS','IBS0','IBS1','IBS2','IBD0','IBD1','IBD2','PIHAT'])

  sys.stderr.write('Estimating pairwise IBS and IBD...\n')

  if options.testpairs:
    pairs,pair_count = file_pairs(genos, options.testpairs)
  else:
    pairs = pair_generator(genos)
    pair_count = len(genos)*(len(genos)-1)//2

  threshold = options.threshold

  e00 = float(ibs_given_ibd[0,0])
  e01 = float(ibs_given_ibd[0,1])
  e11 = float(ibs_given_ibd[1,1])

  if options.progress:
    pairs = progress_bar(pairs, pair_count)

  for (sample1,genos1),(sample2,genos2) in pairs:
    ibs0,ibs1,ibs2 = genoarray_ibs(genos1,genos2)
    n = ibs0+ibs1+ibs2

    if not n:
      continue

    ibs0 = ibs0/n
    ibd0 = min(1,ibs0/e00)

    if ibd0 >= threshold:
      continue

    ibs1  = ibs1/n
    ibs2  = ibs2/n
    ibd1  = max(0, min((ibs1 - ibd0*e01)/e11, 1))
    ibd2  = max(0,1-ibd0-ibd1)
    pihat = ibd1/2 + ibd2

    if sample1>sample2:
      sample1,sample2=sample2,sample1

    estimates = [ibs0,ibs1,ibs2,ibd0,ibd1,ibd2,pihat]
    out.writerow([sample1,sample2,n]+[ '%0.4f' % e for e in estimates ])


if __name__=='__main__':
  main()
