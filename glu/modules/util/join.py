# -*- coding: utf-8 -*-
'''
File:          join.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2008-03-06

Abstract:      Utility to join two files on a common lookup key

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   glu.lib.fileutils import load_table,load_list,tryint,autofile,hyphen


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] table list1 [list2]...'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  parser.add_option('-1', '--key1', dest='key1',
                    help='Key column number or name')
  parser.add_option('-2', '--key2', dest='key2',
                    help='Key column number or name')
  parser.add_option('-s', '--skipkey', dest='skipkey', action='store_true',
                    help='Skip key2 in output')
  return parser


def strip(seq):
  return (i.strip() for i in seq)


def load(filename,key):
  table = load_table(hyphen(filename,sys.stdin),want_header=True)
  header = table.next()
  table = list(table)

  key = tryint(key or 0)
  if isinstance(key, basestring):
    try:
      key = header.index(key)
    except ValueError:
      raise ValueError("Cannot find header '%s'" % key)

  return header,table,key


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  lheader,ltable,key1 = load(args[0],options.key1)
  rheader,rtable,key2 = load(args[1],options.key2)

  n = len(rheader)
  if not options.skipkey:
    rmap = dict( (row[key2].strip(),row) for row in rtable )
  else:
    rmap = dict( (row[key2].strip(),row[:key2]+row[key2+1:]) for row in rtable )
    rheader = rheader[:key2]+rheader[key2+1:]
    n -= 1

  out = csv.writer(autofile(hyphen(options.output,sys.stdout),'w'),dialect='excel-tab')
  out.writerow(lheader + rheader)


  blank = ['']*n

  for lrow in ltable:
    key  = lrow[key1].strip()
    rrow = rmap.get(key,blank)

    out.writerow(lrow+rrow)


if __name__=='__main__':
  main()
