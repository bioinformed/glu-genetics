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

from   glu.lib.utils         import autofile,hyphen
from   glu.lib.genoarray     import get_genorepr
from   glu.lib.genomerge     import get_genomerger, output_merge_statistics
from   glu.lib.genodata      import transform_files, GenoTransform, save_genostream


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

  ioopts.add_option('-g', '--ingenorepr', dest='ingenorepr', metavar='REP', type='string', default='snp_acgt',
                    help='Input genotype representation.  Values=snp_acgt (default), snp_ab, snp_marker, or generic')
  ioopts.add_option('-G', '--outgenorepr', dest='outgenorepr', metavar='REP', type='string', default=None,
                    help='Output genotype representation (see -g/--ingenorepr).  Default is ingenorepr')
  ioopts.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
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

  ingenorepr  = get_genorepr(options.ingenorepr)
  outgenorepr = get_genorepr(options.outgenorepr) if options.outgenorepr else ingenorepr

  infiles = sorted(set(hyphen(arg,sys.stdin) for arg in args))
  outfile = hyphen(options.output,sys.stdout)

  merger = get_genomerger(options.merge,outgenorepr)

  transform = GenoTransform.from_options(options)

  genos = transform_files(infiles, options.informat,  ingenorepr,
                          outfile, options.outformat, outgenorepr,
                                   transform, mergefunc=merger,
                                   limit=options.limit)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__ == '__main__':
  main()
