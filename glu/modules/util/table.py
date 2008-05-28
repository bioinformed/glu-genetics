# -*- coding: utf-8 -*-
'''
File:          join.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2008-03-06

Abstract:      Toy module to read and write GLU tables

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   glu.lib.fileutils import load_table,table_writer


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] table'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  table = load_table(args[0],hyphen=sys.stdin,want_header=True)
  out   = table_writer(options.output,hyphen=sys.stdout)

  out.writerows(table)


if __name__=='__main__':
  main()