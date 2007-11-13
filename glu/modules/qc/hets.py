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

from   glu.lib.utils         import percent, tally
from   glu.lib.fileutils     import autofile, hyphen, load_list
from   glu.lib.genolib       import load_genostream
from   glu.lib.sections      import save_section, SectionWriter, save_metadata_section


def het_output(out,results):

  results = [ (l,hom,het,miss,percent(het,hom+het),percent(miss,hom+het+miss))
               for l,(hom,het,miss) in results ]

  results.sort(key=itemgetter(4,3,0))

  out.write('HETEROZYGOSITY')
  out.write('\n')
  out.write('  Rank  Sample                      homs     hets    missing   %het    %missing\n')
  out.write('  ----  -------------------------  -------  -------  -------  -------  --------\n')

  for r,(l,hom,het,miss,p1,p2) in enumerate(results):
    out.write('  %4d  %-25s  %7d  %7d  %7d  %7.3f  %7.3f\n' % (r+1,l,hom,het,miss,p1,p2))

  out.write('\n')


def save_results(sw,results):
  results.sort()
  rows  =[['sample', 'hom', 'hets', 'miss']]
  rows += ([l,hom,het,miss] for r,(l,(hom,het,miss)) in enumerate(results))
  save_section(sw, 'het', rows)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-f', '--format', dest='format',
                    help='Input data format, either hapmap,ldat,sdat or counts')
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
  parser.add_option('-m','--mincount', dest='mincount', metavar='N',default=None, type='int',
                    help='Minimum number of non-missing genotypes')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('--tablularoutput', dest='tablularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')
  return parser


def count_genos(genos):
  '''
  '''
  f = tally(genos)

  hom = het = miss = 0

  for g,n in f.iteritems():
    if not g or None in g:
      miss += n
    elif g[0] != g[1]:
      het  += n
    else:
      hom  += n

  return hom,het,miss


def geno_counts(loci):
  # Skip header
  for lname,genos in loci:
    yield lname,count_genos(genos)


def read_counts(filename):
  loci = csv.reader(autofile(filename),dialect='excel-tab')
  for sample in loci:
    if len(sample) < 3 or not sample[0]:
      continue

    hom,het,miss = tuple(map(int,sample[1:4]))
    yield sample[0],(hom,het,miss)


def filter_counts(counts,mincount):
  for sample,(hom,het,miss) in counts:
    if hom+het >= mincount:
      yield sample,(hom,het,miss)


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
                                   modelmap=options.loci).as_sdat()

    loci = loci.transformed(include_loci=options.includeloci,
                            exclude_loci=options.excludeloci,
                            include_samples=options.includesamples,
                            exclude_samples=options.excludesamples)

    counts = geno_counts(loci)

  else:
    if options.includeloci or options.excludeloci:
      raise ValueError,'Cannot specift a locus filter for count input'

    counts = read_counts(args[0])

    if options.includesamples:
      includesamples = set(load_list(options.includesamples))
      counts = [ c for c in counts if c[0] in includesamples ]

    if options.excludesamples:
      excludesamples = set(load_list(options.excludesamples))
      counts = [ c for c in counts if c[0] not in excludesamples ]

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
