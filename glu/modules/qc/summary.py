# -*- coding: utf-8 -*-
'''
File:          summary.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       2008-05-22

Abstract:      Genotype summary statistics, fast and simple

Requires:      Python 2.5, glu

Revision:      $Id:
'''

from __future__ import division

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys

from   itertools                 import izip

from   glu.lib.fileutils         import table_writer
from   glu.lib.hwp               import hwp_biallelic
from   glu.lib.genolib           import load_genostream
from   glu.lib.genolib.genoarray import locus_summary, sample_summary, \
                                        count_alleles_from_genocounts


LOCUS_HEADER  = ['Locus','chromosome','location','strand',
                 'allele count','alleles','allele counts',
                 'minor allele frequency',
                 'genotype count','genotypes','genotype counts',
                 'missing count', 'informative count', 'missing rate','hw pvalue']

SAMPLE_HEADER = ['Sample','missing count','hemizygote count',
                 'homozygote count','heterozygote count',
                 'informative count', 'missing rate', 'heterozygosity']


def locus_row(lname,locus,model,counts):
  acounts = izip(count_alleles_from_genocounts(model,counts),model.alleles)
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

  n = counts.sum()
  missing_rate = str(counts[0]/n) if n else ''

  try:
    hwp = str(hwp_biallelic(model,counts))
  except ValueError:
    hwp = ''

  return [lname,locus.chromosome or '', str(locus.location or ''), locus.strand or '',
                str(len(alleles)),'|'.join(alleles),'|'.join(str(n) for n in acounts),maf,
                str(len(lgenos)), '|'.join(lgenos), '|'.join(str(n) for n in lcounts),
                counts[0], sum(counts[1:]), missing_rate, hwp]


def locus_total(counts):
  n = counts.sum()
  missing_rate = str(counts[0]/n) if n else ''

  return ['*','','','','','','','','','','',counts[0],sum(counts[1:]),missing_rate,'']


def sample_row(sample,counts):
  n = counts.sum()
  missing_rate = str(counts[0]/n) if n else ''

  n = counts[2]+counts[3]
  heterozygosity = str(counts[2]/(counts[2]+counts[3])) if n else ''

  return [sample]+counts.tolist()+[sum(counts[1:]),missing_rate,heterozygosity]


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-f', '--format', dest='format',
                    help='Input format for genotype or count data')
  parser.add_option('-g', '--genorepr',        dest='genorepr',      metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('-n', '--includesamples', dest='includesamples', metavar='FILE',
                    help='Include list for those samples to only use')
  parser.add_option('-u', '--includeloci', dest='includeloci', metavar='FILE',
                    help='Include list for those loci to only use')
  parser.add_option('-x', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='Exclude a list of samples')
  parser.add_option('-e', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='Exclude a list of loci')
  parser.add_option('-o', '--locusout', dest='locusout', metavar='FILE', default='-',
                    help='Locus output report')
  parser.add_option('-O', '--sampleout', dest='sampleout', metavar='FILE', default='-',
                    help='Locus output report')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                  genome=options.loci)

  genos = genos.transformed(include_loci=options.includeloci,
                            exclude_loci=options.excludeloci,
                            include_samples=options.includesamples,
                            exclude_samples=options.excludesamples)

  if genos.format not in ('sdat','ldat'):
    genos = genos.as_ldat()

  locusout  = table_writer(options.locusout,hyphen=sys.stdout)
  sampleout = table_writer(options.sampleout,hyphen=sys.stdout)

  if genos.format == 'ldat':
    sample_counts = None

    locusout.writerow(LOCUS_HEADER)
    for (lname,geno),model in izip(genos,genos.models):
      locus = genos.genome.get_locus(lname)
      locus_count,sample_counts = locus_summary(geno,sample_counts)

      locusout.writerow(locus_row(lname,locus,model,locus_count))

    sample_totals = sum(sample_counts)
    locusout.writerow(locus_total(sample_totals))

    sampleout.writerow(SAMPLE_HEADER)
    for sample,count in izip(genos.samples,sample_counts):
      sampleout.writerow(sample_row(sample,count))

    sampleout.writerow(sample_row('*',sample_totals))

  else:
    locus_counts  = None
    sample_counts = 0

    sampleout.writerow(SAMPLE_HEADER)
    for sample,geno in genos:
      sample_count,locus_counts = sample_summary(geno,locus_counts)
      sample_counts += sample_count
      sampleout.writerow(sample_row(sample,sample_count))

    sampleout.writerow(sample_row('*',sample_counts))

    locusout.writerow(LOCUS_HEADER)
    for lname,model,locus_count in izip(genos.loci,genos.models,locus_counts):
      locus = genos.genome.get_locus(lname)
      locusout.writerow(locus_row(lname,locus,model,locus_count))

    locusout.writerow(locus_total(sample_counts))


if __name__=='__main__':
  main()
