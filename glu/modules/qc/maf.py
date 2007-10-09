# -*- coding: utf-8 -*-
'''
File:          maf.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs various transformation on a genotype files

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   operator           import itemgetter

from   glu.lib.utils     import tally
from   glu.lib.fileutils import autofile,hyphen
from   glu.lib.genolib   import load_genostream, get_genorepr


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f','--format',  dest='format', metavar='string',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of transformed data (default is "-" for standard out)')

  return parser


def compute_maf(genos):
  counts = tally(a for g in genos if g for a in g if a)

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

  genorepr  = get_genorepr(options.genorepr)
  data   = load_genostream(args[0],options.format,genorepr=genorepr).as_ldat()

  outfile = autofile(hyphen(options.output,sys.stdout),'w')
  out = csv.writer(outfile,dialect='excel-tab')
  for lname,genos in data:
    as,fs = compute_maf(genos)
    out.writerow( [lname] + as + ['%6.4f' % f for f in fs] )


if __name__ == '__main__':
  main()
