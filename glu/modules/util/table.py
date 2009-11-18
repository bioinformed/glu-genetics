# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Manipulate delimited files, including filtering by column, value, and creating indicator variables based on a categorical variable'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   glu.lib.fileutils   import table_reader,table_writer,cook_table


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] table'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-c', '--categorical', dest='categorical', metavar='VAR', action='append',
                    help='Create indicator variables based on values of VAR')
  parser.add_option('--columnexpr', dest='columnexpr', metavar='VAR=EXPR', action='append',
                    help='Add a new column VAR with the value determined by expression EXPR')
  parser.add_option('--includevar', dest='includevar', metavar='VAR=VAL', action='append',
                    help='Include only records with variable VAR equal to VAL')
  parser.add_option('--excludevar', dest='excludevar', metavar='VAR=VAL', action='append',
                    help='Exclude all records with variable VAR equal to VAL')
  parser.add_option('--filterexpr', dest='filterexpr', metavar='EXPR', action='append',
                    help='Filter all records where EXPR is not true')
  parser.add_option('-s', '--sort', dest='sort', metavar='VAR', action='append',
                    help='Sort rows based on values in column VAR')
  parser.add_option('-u', '--uniq', dest='uniq', action='store_true',
                    help='Produce only unique rows by collapsing consecutive duplicate rows')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  table = table_reader(args[0],hyphen=sys.stdin,want_header=True)
  out   = table_writer(options.output,hyphen=sys.stdout)

  table = cook_table(table,options)

  out.writerows(table)


if __name__=='__main__':
  main()
