# -*- coding: utf-8 -*-

__abstract__  = 'format pairwise LD values into a matrix'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import csv
import sys

from   operator          import itemgetter
from   itertools         import islice, chain, izip

from   glu.lib.fileutils import load_table,table_writer
from   glu.lib.xtab      import xtab


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',   dest='outfile',   metavar='FILE',   default= '-',
                    help='Output file for formatted data')
  parser.add_option('-m', '--measure',  dest='measure', default='r2',
                    help="Measure of LD: r2 (default) or D'")

  return parser


def merge(i,j,x):
  if x:
    return x[0]
  else:
    return ''


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  out = table_writer(options.outfile,hyphen=sys.stdout)

  datain = [ islice(load_table(arg,hyphen=sys.stdin),1,None) for arg in args ]
  datain = chain(*datain)

  if options.measure.lower() == 'r2':
    col = 2
  elif options.measure.lower() == "d'":
    col = 3
  else:
    raise ValueError('Unknown or unsupported LD measure specified: %s' % options.measure)

  columns,rows,data = xtab(datain,itemgetter(1),itemgetter(0),itemgetter(col),merge)
  out.writerow(['']+columns)

  for label,row in izip(rows,data):
    out.writerow([label]+row)


if __name__=='__main__':
  main()
