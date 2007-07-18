# -*- coding: utf-8 -*-
'''
File:          hets.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-01

Abstract:      Compute genotype heterozygosity

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   textwrap              import fill
from   operator              import itemgetter

from   glu.lib.genoarray     import snp_marker
from   glu.lib.genodata      import load_list, load_genomatrixstream
from   glu.lib.utils         import percent, tally
from   glu.lib.fileutils     import autofile, hyphen
from   glu.lib.sections      import save_section, SectionWriter, save_metadata_section


def het_output(out,results):

  results = [ (l,hom,het,hom+het,percent(het,hom+het)) for l,(hom,het) in results ]

  results.sort(key=itemgetter(4,3,0))

  out.write('DEVIATIONS FROM HARDY-WEINBERG PROPORTIONS')
  out.write('\n')
  out.write('  Rank  Sample                     hom    hets     n     %het \n')
  out.write('  ----  -------------------------  -----  -----  -----  ------\n')

  for r,(l,hom,het,n,p) in enumerate(results):
    out.write('  %4d  %-25s  %5d  %5d  %5d  %6.2f\n' % (r+1,l,hom,het,n,p))

  out.write('\n')


def save_results(sw,results):
  results.sort()

  rows=[['sample', 'hom', 'hets']]
  for r,(l,(hom,het)) in enumerate(results):
    rows.append([l,hom,het])

  save_section(sw, 'het', rows)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-s', '--samplesubset',   dest='samplesubset', metavar='FILE',
                    help='List of samples to include in the analysis')
  parser.add_option('-l', '--locussubset',    dest='locussubset',  metavar='FILE',
                    help='List of loci to include in the analysis')
  parser.add_option('-m','--mincount', dest='mincount', metavar='N',default=None, type='int',
                    help='Minimum number of non-missing genotypes')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-f', '--format', dest='format', default='sdat',
                    help='Input data format, either sdat or counts (default=sdat)')
  parser.add_option('-L', '--limit', dest='limit', metavar='N', type='int',
                    help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')
  parser.add_option('--tablularoutput', dest='tablularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')
  return parser


def count_genos(genos):
  '''
  '''
  f = tally(g for g in genos if g and ' ' not in g)

  hom = het = 0

  for g,n in f.iteritems():
    if g[0] != g[1]:
      het += n
    else:
      hom += n

  return hom,het


def geno_counts(loci):
  # Skip header
  loci = iter(loci)
  loci.next()

  for lname,genos in loci:
    yield lname,count_genos(genos)


def read_counts(filename):
  loci = csv.reader(autofile(filename),dialect='excel-tab')
  for sample in loci:
    if len(sample) < 3 or not sample[0]:
      continue

    hom,het = tuple(map(int,sample[1:3]))
    yield sample[0],(hom,het)


def filter_counts(counts,mincount):
  for sample,(hom,het) in counts:
    if hom+het >= mincount:
      yield sample,(hom,het)


def main():
  import sys,time

  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  if options.format != 'counts':
    loci = load_genomatrixstream(args[0],options.format,limit=options.limit,genorepr=snp_marker).as_ldat()

    if options.locussubset:
      loci = loci.transformed(exclude_loci=options.locussubset)

    counts = geno_counts(loci)

  else:
    if options.locussubset:
      raise ValueError,'Cannot specift a locus filter for count input'

    counts = read_counts(args[0])

  if options.samplesubset:
    samplesubset = set(load_list(options.samplesubset))
    counts = [ c for c in counts if c[0] in samplesubset ]

  if options.mincount:
    counts = filter_counts(counts,options.mincount)

  het_output(out,counts)

  if options.tablularoutput:
    sw = SectionWriter(options.tablularoutput)
    save_metadata_section(sw, analysis='het', analysis_version='0.1')
    save_results(sw, results)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
