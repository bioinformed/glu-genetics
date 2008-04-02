# -*- coding: utf-8 -*-
'''
File:          hwp.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-01

Abstract:      Test genotype data for deviations from Hardy-Weinbery proportions

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   glu.lib.utils             import izip_exact
from   glu.lib.fileutils         import autofile, hyphen, load_list
from   glu.lib.hwp               import hwp_exact_biallelic, hwp_chisq_biallelic, biallelic_counts, HWP_EXACT_THRESHOLD
from   glu.lib.sections          import save_section, SectionWriter, save_metadata_section

from   glu.lib.genolib           import load_genostream
from   glu.lib.genolib.genoarray import count_genotypes


def test_loci_hwp(data):
  for lname,_,counts in data:
    if counts[0]*2 + counts[1] < HWP_EXACT_THRESHOLD:
      pexact = hwp_exact_biallelic(*counts)
    else:
      pexact = None

    pasymp = hwp_chisq_biallelic(*counts)

    yield lname,pexact,pasymp,counts[0],counts[1],counts[2]


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
  parser.add_option('-g', '--genorepr',        dest='genorepr',        metavar='REP',
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
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=50, type='int',
                    help='Exclude loci with less than N non-missing genotypes (default=50)')
  parser.add_option('--mincompletion', dest='mincompletion', metavar='N', default=0, type='float',
                    help='Exclude loci with genotype completion rate less than N (default=0 for no no exclusion)')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0, type='float',
                    help='Exclude loci with minor allele frequency less than N (default=0 for no exclusion)')
  parser.add_option('--tablularoutput', dest='tablularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')
  return parser


def geno_counts(loci):
  for (lname,genos),model in izip_exact(loci,loci.models):
    counts = count_genotypes(model,genos)
    yield lname,counts,biallelic_counts(model,counts)


def read_counts(filename):
  loci = csv.reader(autofile(filename),dialect='excel-tab')
  for locus in loci:
    if len(locus) < 4 or not locus[0]:
      continue

    hom1,het,hom2 = tuple(map(int,locus[1:4]))
    hom1,hom2 = min(hom1,hom2),max(hom1,hom2)
    yield locus[0],None,(hom1,het,hom2)


def filter_counts(data, mingenos, minmaf, mincompletion):
  maxmissing = 1-mincompletion

  for locus,counts,bialleleic_counts in data:
    n = sum(counts)
    m = sum(bialleleic_counts)

    if sum(bialleleic_counts) < mingenos:
      continue

    # Completion does count hemizygotes
    if mincompletion>0 and (not n or float(counts[0])/n > maxmissing):
      continue

    if minmaf>0:
      if not m:
        continue

      # Compute MAF based on non-hemizygous counts only or else users will
      # be confused
      maf = float(min(2*bialleleic_counts[0]+bialleleic_counts[1],
                      2*bialleleic_counts[2]+bialleleic_counts[1])) / m

      if maf < minmaf:
        continue

    yield locus,counts,bialleleic_counts


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
                                   genome=options.loci).as_ldat()

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

  if options.mingenos or options.minmaf or options.mincompletion:
    counts = filter_counts(counts, options.mingenos, options.minmaf, options.mincompletion)

  results = test_loci_hwp(counts)

  hwp_output(out,results)

  if options.tablularoutput:
    sw = SectionWriter(options.tablularoutput)
    save_metadata_section(sw, analysis='hwp', analysis_version='0.1')
    save_results(sw, results)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
