# -*- coding: utf-8 -*-
'''
File:          transform.py

Authors:       Brian Staats  (staatsb@mail.nih.gov)
               Xiang Deng      (dengx@mail.nih.gov)
               Jun Lu          (lujun@mail.nih.gov)
               Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs various transformation on a genotype files

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   glu.lib.fileutils       import autofile,hyphen
from   glu.lib.genolib.merge   import get_genomerger, output_merge_statistics
from   glu.lib.genolib.streams import GenoTransform
from   glu.lib.genolib.io      import transform_files, save_genostream


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage)

  ioopts = optparse.OptionGroup(parser, 'Input/Output Options')

  ioopts.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output transformed data to FILE(default is "-" for standard out)')
  ioopts.add_option('-f','--informat',  dest='informat',
                    help='Input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  ioopts.add_option('-F','--outformat',  dest='outformat',
                    help='Output format for genotype data. Default is informat, cannot be hapmap.')

  ioopts.add_option('-g', '--ingenorepr', dest='ingenorepr', metavar='REP',
                    help='Input genotype representation')
  ioopts.add_option('-G', '--outgenorepr', dest='outgenorepr', metavar='REP', default=None,
                    help='Output genotype representation (see -g/--ingenorepr).  Default is ingenorepr')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  ioopts.add_option(      '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows of data to N for testing purposes')

  mopts = optparse.OptionGroup(parser, 'Genotype Merging and Reporting')

  mopts.add_option('--merge', dest='merge', metavar='METHOD:T', default='vote:1',
                    help='Genotype merge algorithm and optional consensus threshold used to form a consensus genotypes. '
                         'Values=vote,ordered,unique.  Value may be optionally followed by a colon and a threshold.  Default=vote:1')
  mopts.add_option('--samplemerge', dest='samplemerge', metavar='FILE',
                    help='Per sample concordance statistics output to FILE (optional)')
  mopts.add_option('--locusmerge',  dest='locusmerge',  metavar='FILE',
                    help='Per locus concordance statistics output to FILE (optional)')

  filter = optparse.OptionGroup(parser, 'Filtering and Renaming')

  filter.add_option('-c', '--filtermissing', action='store_true', dest='filtermissing',
                    help='Filters out the samples or loci with missing genotypes')

  filter.add_option('-n', '--includesamples', dest='includesamples', metavar='FILE',
                    help='Include list for those samples to only use in the transformation and output')
  filter.add_option('-u', '--includeloci', dest='includeloci', metavar='FILE',
                    help='Include list for those loci to only use in the transformation and output')

  filter.add_option('-x', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='Exclude a list of samples from the transformation and output')
  filter.add_option('-e', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='Exclude a list of loci from the transformation and output')

  filter.add_option('-m', '--renamesamples', dest='renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  filter.add_option('-r', '--renameloci', dest='renameloci', metavar='FILE',
                    help='Rename loci from a file containing rows of original name, tab, new name')

  trans = optparse.OptionGroup(parser, 'Transformations')

  trans.add_option('-d', '--ordersamples', dest='ordersamples', metavar='FILE',
                    help='Order samples based on the order of names in FILE')
  trans.add_option('-D', '--orderloci', dest='orderloci', metavar='FILE',
                    help='Order loci based on the order of names in FILE')

  trans.add_option('-a', '--renamealleles', dest='renamealleles', metavar='FILE',
                    help='Rename alleles based on file of locus name, tab, old alleles (comma separated), '
                         'tab, new alleles (comma separated)')

  parser.add_option_group(ioopts)
  parser.add_option_group(mopts)
  parser.add_option_group(filter)
  parser.add_option_group(trans)

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  if not options.output:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify output file'
    return

  infiles = sorted(set(args))
  merger  = get_genomerger(options.merge)

  transform = GenoTransform.from_options(options)

  transform_files(infiles, options.informat,  options.ingenorepr,
                  options.output, options.outformat, options.outgenorepr,
                  transform, mergefunc=merger, genome=options.loci,
                  limit=options.limit,inhyphen=sys.stdin,outhyphen=sys.stdout)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__ == '__main__':
  main()
