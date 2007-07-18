# -*- coding: utf-8 -*-
'''
File:          hwp.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-01

Abstract:      Test genotype data for deviations from Hardy-Weinbery proportions

Requires:      Python 2.5, biozilla

Revision:      $Id: hwp.py 463 2006-12-11 14:58:41Z staatsb $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv

from   itertools              import islice, izip, chain
from   operator               import itemgetter
from   textwrap               import fill

from   ill2raw                import load_list
from   biozilla.utils         import any, all, autofile, hyphen, percent, tally
from   biozilla.hwp           import hwp_exact_biallelic, hwp_chisq_biallelic, count_genos


# FIXME: Merge constants for maximum use of exact test with biozilla.hwp_biallelic
def test_loci_hwp(counts):
  for lname,count in counts:
    if count[0]*2 + count[1] < 10000:
      pexact = hwp_exact_biallelic(*count)
    else:
      pexact = None

    pasymp = hwp_chisq_biallelic(*count)

    yield lname,pexact,pasymp,count[0],count[1],count[2]


def hwp_output(out,results):

  def keyfunc(r):
    if r[1] is not None:
      return r[1],r[0]
    return r[2],r[0]

  results = list(results)
  results.sort(key=keyfunc)

  out.write('DEVIATIONS FROM HARDY-WEINBERG PROPORTIONS')
  out.write('\n')
  out.write('                                                     Exact       Aymptotic \n')
  out.write('  Rank  Locus            hom1s  hets   hom2s    n    p-value     p-value   \n')
  out.write('  ----  ---------------  -----  -----  -----  -----  ----------  ----------\n')

  for r,(l,p,p2,hom1,het,hom2) in enumerate(results):
    n = hom1+het+hom2
    pe = '%10.8f' % p if p is not None else '   N/A    '
    out.write('  %4d  %-15s  %5d  %5d  %5d  %5d  %s  %10.8f\n' % (r+1,l,hom1,het,hom2,n,pe,p2))

  out.write('\n')


# FIXME: Merge with genodata
# NOTE: Generic filter matrix row by label implementation
def filter_loci(loci, include):
  loci = iter(loci)
  yield loci.next()
  for locus in loci:
    if locus[0] in include:
      yield locus


# FIXME: Merge with genodata
# NOTE: Generic filter matrix column by label implementation
def filter_samples(loci, include):
  loci = iter(loci)
  header = loci.next()
  indices = [ i for i,h in enumerate(header) if not i or h in include ]
  loci = chain([header],loci)

  if len(indices) == len(header):
    return loci

  return ( [ locus[i] for i in indices ] for locus in loci )


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)
  parser.add_option('-s', '--samplesubset',   dest='samplesubset', metavar='FILE',
                    help='List of samples to include in the analysis')
  parser.add_option('-l', '--locussubset',    dest='locussubset',  metavar='FILE',
                    help='List of loci to include in the analysis')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-f', '--format', dest='format', type='string', default='ldat',
                    help='Input data format, either ldat or counts (default=ldat)')
  parser.add_option('-L', '--limit', dest='limit', metavar='N', type='int', default=0,
                    help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')
  return parser


def geno_counts(loci):
  for locus in loci:
    genos = iter(locus)
    lname = genos.next()
    yield lname,count_genos(genos)


def read_counts(loci):
  for locus in loci:
    if len(locus) != 4 or not locus[0]:
      continue

    hom1,het,hom2 = tuple(map(int,locus[1:]))
    hom1,hom2 = min(hom1,hom2),max(hom1,hom2)
    yield locus[0],(hom1,het,hom2)


def main():
  import sys,time

  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  loci = csv.reader(autofile(args[0]),dialect='excel-tab')

  if options.limit:
    loci = islice(loci, options.limit+1)

  if options.locussubset:
    locussubset = set(load_list(options.locussubset))
    loci = filter_loci(loci,locussubset)

  if options.format == 'ldat':
    if options.samplesubset:
      samplesubset = set(load_list(options.samplesubset))
      loci = filter_samples(loci, samplesubset)

    header = loci.next()
    counts = geno_counts(loci)

  elif options.format == 'counts':
    if options.samplesubset:
      raise ValueError,'Cannot specift a sample filter for count input'

    counts = read_counts(loci)

  print >> sys.stderr, 'Checking HWP...',

  results = test_loci_hwp(counts)

  hwp_output(out,results)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
