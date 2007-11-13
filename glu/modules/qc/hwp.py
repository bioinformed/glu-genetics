# -*- coding: utf-8 -*-
'''
File:          hwp.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-01

Abstract:      Test genotype data for deviations from Hardy-Weinbery proportions

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   textwrap              import fill

from   glu.lib.utils         import percent
from   glu.lib.fileutils     import autofile, hyphen, load_list
from   glu.lib.genolib       import load_genostream
from   glu.lib.hwp           import hwp_exact_biallelic, hwp_chisq_biallelic, count_genos
from   glu.lib.sections      import save_section, SectionWriter, save_metadata_section


# FIXME: Merge constants for maximum use of exact test with glu.lib.hwp_biallelic
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
  out.write('                                                     Exact       Asymptotic\n')
  out.write('  Rank  Locus            hom1s  hets   hom2s    n    p-value     p-value   \n')
  out.write('  ----  ---------------  -----  -----  -----  -----  ----------  ----------\n')

  for r,(l,p,p2,hom1,het,hom2) in enumerate(results):
    n = hom1+het+hom2
    pe = ('%10.8f' % p) if p is not None else '   N/A    '
    out.write('  %4d  %-15s  %5d  %5d  %5d  %5d  %s  %10.8f\n' % (r+1,l,hom1,het,hom2,n,pe,p2))

  out.write('\n')


def save_results(sw,results):
  results.sort()

  rows=[['locus', 'hom1s', 'hets', 'hom2s', 'num', 'p-value', 'aymptotic p-value']]
  for r,(l,p,p2,hom1,het,hom2) in enumerate(results):
    n = hom1+het+hom2
    pe = repr(p) if p is not None else ''
    rows.append([l,hom1,het,hom2,n,pe,repr(p2)])

  save_section(sw, 'hwp', rows)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of tests for deviation from Hardy-Weinberg proportions')
  parser.add_option('-f', '--format', dest='format', metavar='F',
                    help='Input data format, ldat, sdat, trip, or counts')
  parser.add_option('-g', '--genorepr',        dest='genorepr',        metavar='REPR', default='snp',
                    help='Input genotype representations. Values=snp (default), hapmap, or marker')
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
  parser.add_option(      '--limit', dest='limit', metavar='N', type='int',
                    help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')
  parser.add_option('--tablularoutput', dest='tablularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')
  return parser


def geno_counts(loci):
  for lname,genos in loci:
    yield lname,count_genos(genos)


def read_counts(filename):
  loci = csv.reader(autofile(filename),dialect='excel-tab')
  for locus in loci:
    if len(locus) < 4 or not locus[0]:
      continue

    hom1,het,hom2 = tuple(map(int,locus[1:4]))
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

  if options.format != 'counts':
    loci = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                   modelmap=options.loci,limit=options.limit).as_ldat()

    loci = loci.transformed(include_loci=options.includeloci,
                            exclude_loci=options.excludeloci,
                            include_samples=options.includesamples,
                            exclude_samples=options.excludesamples)

    counts = geno_counts(loci)

  else:
    if options.options.includesamples or options.excludesamples:
      raise ValueError,'Cannot specify a sample filter for count input'

    counts = read_counts(args[0])

    if options.includeloci:
      includeloci = set(load_list(options.includeloci))
      counts = [ c for c in counts if c[0] in includeloci ]

    if options.excludeloci:
      excludeloci = set(load_list(options.excludeloci))
      counts = [ c for c in counts if c[0] not in excludeloci ]

  print >> sys.stderr, 'Checking HWP...',

  results = test_loci_hwp(counts)

  hwp_output(out,results)

  if options.tablularoutput:
    sw = SectionWriter(options.tablularoutput)
    save_metadata_section(sw, analysis='hwp', analysis_version='0.1')
    save_results(sw, results)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
