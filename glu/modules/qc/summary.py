# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute a battery of genotype summary statistics'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools                 import izip

import numpy as np

from   glu.lib.fileutils         import autofile,hyphen,list_reader,table_writer
from   glu.lib.hwp               import hwp_biallelic
from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.genoarray import locus_summary, sample_summary, \
                                        count_alleles_from_genocounts


LOCUS_HEADER  = ['LOCUS','CHROMOSOME','LOCATION','STRAND',
                 'NUM_ALLELES','ALLELES','ALLELE_COUNTS', 'MAF',
                 'NUM_GENOTYPES','GENOTYPES','GENOTYPE_COUNTS',
                 'MISSING_COUNT', 'INFORMATIVE_COUNT',
                 'ATTEMPTED_MISSING_RATE', 'OBSERVED_MISSING_RATE',
                 'NONEMPTY_MISSING_RATE', 'HW_PVALUE']

SAMPLE_HEADER = ['SAMPLE','MISSING_COUNT','HEMIZYGOTE_COUNT',
                 'HOMOZYGOTE_COUNT','HETEROZYGOTE_COUNT',
                 'INFORMATIVE_COUNT', 'ATTEMPTED_MISSING_RATE',
                 'OBSERVED_MISSING_RATE', 'NONEMPTY_MISSING_RATE',
                 'HETEROZYGOSITY']


def rate(a,b):
  return a/b if b!=0 else ''


def missing_rates(observed_missing,observed,empty,missing):
  return [ rate(observed_missing+missing,observed+missing),
           rate(observed_missing,        observed),
           rate(observed_missing-empty,  observed-empty) ]


def locus_row(lname,locus,counts,empty_samples,missing_samples,compute_hwp):
  model   = locus.model
  m       = len(model.genotypes)

  # Counts can be padded with extra zeros that must be removed
  assert not sum(counts[m:])
  counts  = counts[:m]

  # Count alleles
  acounts = count_alleles_from_genocounts(model,counts)
  acounts = izip(acounts,model.alleles) if acounts else []
  acounts = sorted(acounts,reverse=True)
  alleles = [ a for n,a in acounts if n if a and n ]
  acounts = [ n for n,a in acounts if n if a and n ]

  if not alleles:
    maf = ''
  elif len(alleles) == 1:
    maf = 0.0
  elif len(alleles) == 2:
    maf = acounts[-1]/sum(acounts)
  else:
    maf = ''

  alen    = max(len(a) for a in alleles) if alleles else 0
  delim   = '/' if alen>1 else ''
  lgenos  = sorted(izip(counts[1:],model.genotypes[1:]),reverse=True)
  lcounts = [ n for n,g in lgenos if n ]
  lgenos  = [ (delim.join( [g[0] or '',g[1] or ''])) for n,g in lgenos if n ]

  missing = counts[0]
  total   = counts.sum()
  rates   = missing_rates(missing,total,empty_samples,missing_samples)

  hwp = ''
  if compute_hwp:
    try:
      hwp = str(hwp_biallelic(model,counts))
    except ValueError:
      pass

  return [lname,locus.chromosome or '', str(locus.location or ''), locus.strand or '',
                str(len(alleles)),'|'.join(alleles),'|'.join(str(n) for n in acounts),maf,
                str(len(lgenos)), '|'.join(lgenos), '|'.join(str(n) for n in lcounts),
                missing, sum(counts[1:])]+rates+[hwp]


def locus_total(counts,empty_samples,missing_samples):
  missing = counts[0]
  total   = counts.sum()
  rates   = missing_rates(missing,total,empty_samples,missing_samples)

  return ['*','','','','','','','','','','',
          missing,sum(counts[1:])]+rates+['']


def sample_row(sample,counts,empty_loci,missing_loci):
  missing = counts[0]
  total   = counts.sum()
  rates   = missing_rates(missing,total,empty_loci,missing_loci)

  heterozygosity = rate(counts[3],counts[2]+counts[3])

  return ([sample]+counts.tolist()+[sum(counts[1:])]+rates+[heterozygosity])


def summarize(genos):
  if genos.format not in ('sdat','ldat'):
    genos = genos.as_ldat()

  if genos.format == 'ldat':
    loci          = []
    locus_counts  = []
    samples       = genos.samples
    sample_counts = None

    for (lname,geno) in genos:
      locus_count,sample_counts = locus_summary(geno,sample_counts)
      loci.append(lname)
      locus_counts.append(locus_count)

  else:
    samples       = []
    sample_counts = []
    loci          = genos.loci
    locus_counts  = None

    for sample,geno in genos:
      sample_count,locus_counts = sample_summary(geno,locus_counts)
      samples.append(sample)
      sample_counts.append(sample_count)

  if locus_counts is None:
    locus_counts  = []

  if sample_counts is None:
    sample_counts = []

  return loci,locus_counts,samples,sample_counts


def count_empty(counts,n):
  empty = 0
  for count in counts:
    if count[0]==n:
      empty += 1
  return empty


def format_rate(numerator, denominator):
  if denominator<=0:
    rate = 0
  else:
    rate = numerator/denominator

  return (numerator,denominator,rate)


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('--hwp', action='store_true',
                    help='Test for deviation from Hardy-Weinberg proportions')
  parser.add_argument('-s', '--summaryout', metavar='FILE', default='-',
                    help='Summary output file name')
  parser.add_argument('-o', '--locusout', metavar='FILE',
                    help='Locus output table file name')
  parser.add_argument('-O', '--sampleout', metavar='FILE',
                    help='Sample output table file name')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  # Include lists are used to communicate the universe of attempted samples/loci
  # Any not observed in the genotype data are classified as "missing"
  includeloci    = set(list_reader(options.includeloci))    if options.includeloci    else None
  includesamples = set(list_reader(options.includesamples)) if options.includesamples else None

  options.includeloci    = includeloci
  options.includesamples = includesamples

  genos = load_genostream(options.genotypes,format=options.informat,genorepr=options.ingenorepr,
                          genome=options.loci,phenome=options.pedigree,
                          transform=options, hyphen=sys.stdin)

  loci,locus_counts,samples,sample_counts = summarize(genos)
  sample_totals = sum(sample_counts, np.zeros(4))

  assert len(loci)*len(samples) == sample_totals.sum()

  # Locus statistics
  missing_loci        = len(includeloci-set(loci)) if includeloci else 0
  attempted_loci      = len(loci)+missing_loci
  empty_loci          = count_empty(locus_counts, len(samples))
  nonempty_loci       = len(loci)-empty_loci

  # Sample statistics
  missing_samples     = len(includesamples-set(samples)) if includesamples else 0
  attempted_samples   = len(samples)+missing_samples
  empty_samples       = count_empty(sample_counts,len(loci))
  nonempty_samples    = len(samples)-empty_samples

  # Observed genotype statistics
  missing_genotypes   = sample_totals[0]
  observed_genotypes  = len(loci)*len(samples)
  found_genotypes     = observed_genotypes - missing_genotypes

  # Non-empty genotype statistics
  nonempty_genotypes  = nonempty_loci*nonempty_samples
  empty_genotypes     = observed_genotypes - nonempty_genotypes

  # Attempted genotype statistics
  attempted_genotypes = attempted_loci*attempted_samples
  phantom_genotypes   = attempted_genotypes - observed_genotypes

  if options.summaryout:
    summaryout = autofile(hyphen(options.summaryout,sys.stdout),'w')
    summaryout.write('loci   : attempted=%8d, observed=%8d, empty=%8d, missing=%8d\n'
                         % (attempted_loci,len(loci),empty_loci,missing_loci))
    summaryout.write('samples: attempted=%8d, observed=%8d, empty=%8d, missing=%8d\n'
                         % (attempted_samples,len(samples),empty_samples,missing_samples))
    summaryout.write('\n')
    summaryout.write('Attempted genotypes: missing=%8d, total=%8d, rate=%g\n'
                         % format_rate(attempted_genotypes-found_genotypes,attempted_genotypes))
    summaryout.write('Observed  genotypes: missing=%8d, total=%8d, rate=%g\n'
                         % format_rate(missing_genotypes,observed_genotypes))
    summaryout.write('Non-empty genotypes: missing=%8d, total=%8d, rate=%g\n'
                         % format_rate(missing_genotypes-empty_genotypes,nonempty_genotypes))

  if options.sampleout:
    sampleout = table_writer(options.sampleout,hyphen=sys.stdout)
    sampleout.writerow(SAMPLE_HEADER)
    for sample,count in izip(samples,sample_counts):
      sampleout.writerow(sample_row(sample,count,empty_loci,missing_loci))
    sampleout.writerow(sample_row('*',sample_totals,empty_genotypes,phantom_genotypes))
    del sampleout

  if options.locusout:
    locusout  = table_writer(options.locusout,hyphen=sys.stdout)
    locusout.writerow(LOCUS_HEADER)
    for lname,locus_count in izip(loci,locus_counts):
      locus = genos.genome.get_locus(lname)
      locusout.writerow(locus_row(lname,locus,locus_count,empty_samples,missing_samples,options.hwp))
    locusout.writerow(locus_total(sample_totals,empty_genotypes,phantom_genotypes))
    del locusout


if __name__=='__main__':
  main()
