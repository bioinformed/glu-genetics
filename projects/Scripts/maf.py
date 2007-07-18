# -*- coding: utf-8 -*-
'''
File:          transform2.py

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

from   operator               import itemgetter
from   itertools              import chain

from   biozilla.utils         import autofile,hyphen,tally
from   biozilla.genoarray     import get_genorepr,snp_marker
from   biozilla.genomerge     import get_genomerger, output_merge_statistics
from   biozilla.genodata      import load_genostream, guess_informat_list, guess_outformat


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  ioopts = optparse.OptionGroup(parser, 'Input/Output Options')

  ioopts.add_option('-f','--format',  dest='format', metavar='string',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  ioopts.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of transformed data (default is "-" for standard out)')

  ioopts.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', type='string', default='snp_acgt',
                    help='Input genotype representation.  Values=snp_acgt (default), snp_ab, snp_marker, or generic')

  ioopts.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows of data to N for testing purposes')

  parser.add_option_group(ioopts)

  return parser


def compute_maf(genos):
  counts = tally( chain(*(g for g in genos if g) ))
  counts.pop(' ',None)

  if not counts:
    return ['',''],[0,0]

  n = sum(counts.itervalues())
  counts = sorted(counts.iteritems(),key=itemgetter(1))
  as,fs = zip(*counts) + [('',0)]*(len(counts)-2)
  fs = [ float(f)/n for f in fs ]
  return list(as),fs


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  if not options.format:
    options.format = guess_informat_list(args)

  genorepr  = get_genorepr(options.genorepr)

  print >> sys.stderr, 'INPUT : format=%s,repr=%s' % (options.format, options.genorepr)

  data   = load_genostream(args[0],options.format,limit=options.limit,genorepr=genorepr).recoded(snp_marker).as_ldat()
  data   = iter(data)
  header = data.next()

  outfile = hyphen(options.output,sys.stdout)
  out = csv.writer(outfile,dialect='excel-tab')
  for lname,genos in data:
    as,fs = compute_maf(genos)
    out.writerow( [lname] + as + ['%6.4f' % f for f in fs] )


if __name__ == '__main__':
  main()
