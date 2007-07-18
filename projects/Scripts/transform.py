# -*- coding: utf-8 -*-
'''
File:          transform.py

Authors:       Brian Staats  (staatsb@mail.nih.gov)
               Xiang Deng      (dengx@mail.nih.gov)
               Jun Lu          (lujun@mail.nih.gov)
               Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs various transformation on a genotype files

Requires:      Python 2.5, biozilla

Revision:      $Id: $
'''

__version__ = '0.99'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv

from   biozilla.utils         import autofile,hyphen
from   biozilla.genoarray     import get_genorepr
from   biozilla.genomerge     import get_genomerger, output_merge_statistics
from   biozilla.genodata      import transform_files, GenoTransform, save_genostream


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  ioopts = optparse.OptionGroup(parser, 'Input/Output Options')

  ioopts.add_option('-f','--informat',  dest='informat', metavar='string',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  ioopts.add_option('-F','--outformat',  dest='outformat', metavar='string',
                    help='The file output format for genotype data. Values=ldat (default), sdat, trip or genotriple')

  ioopts.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of transformed data (default is "-" for standard out)')

  ioopts.add_option('-g', '--ingenorepr', dest='ingenorepr', metavar='REP', type='string', default='snp_acgt',
                    help='Input genotype representation.  Values=snp_acgt (default), snp_ab, snp_marker, or generic')
  ioopts.add_option('-G', '--outgenorepr', dest='outgenorepr', metavar='REP', type='string', default=None,
                    help='Output genotype representation (see -r/--ingenorepr).  Default is ingenorepr')

  ioopts.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows of data to N for testing purposes')

  mopts = optparse.OptionGroup(parser, 'Genotype Merge Options')

  mopts.add_option('--merge', dest='merge', metavar='METHOD:T', default='vote:1',
                    help='Genotype merge algorithm and optional consensus threshold used to form a consensus genotypes. '
                         'Values=vote,ordered.  Value may be optionally followed by a colon and a threshold.  Default=vote:1')
  mopts.add_option('--samplemerge', dest='samplemerge', metavar='FILE',
                    help='Output the concordance statistics by sample to FILE (optional)')
  mopts.add_option('--locusmerge',  dest='locusmerge',  metavar='FILE',
                    help='Output the concordance statistics by locus to FILE (optional)')

  filter = optparse.OptionGroup(parser, 'Filtering & Transformation Options')

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
                    help='A list of labels to rename samples')
  filter.add_option('-r', '--renameloci', dest='renameloci', metavar='FILE',
                    help='A list of labels to rename loci')

  filter.add_option('-d', '--ordersamples', dest='ordersamples', metavar='FILE',
                    help='A list of labels to reorder samples')
  filter.add_option('-D', '--orderloci', dest='orderloci', metavar='FILE',
                    help='A list of labels to reorder loci')

  parser.add_option_group(ioopts)
  parser.add_option_group(mopts)
  parser.add_option_group(filter)

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

  merger = get_genomerger(options.merge)

  transform = GenoTransform.from_options(options)

  genos = transform_files(infiles, options.informat,  ingenorepr,
                          outfile, options.outformat, outgenorepr,
                                   transform, mergefunc=merger,
                                   limit=options.limit)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__ == '__main__':
  main()
