# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Manipulate delimited files, including filtering by column, value, and creating indicator variables based on a categorical variable'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   glu.lib.fileutils   import table_reader,table_writer,cook_table,table_options


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('table', help="File name or '-' for stdin", default='-')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')

  table_options(parser)

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  table = table_reader(options.table,hyphen=sys.stdin,want_header=True)
  out   = table_writer(options.output,hyphen=sys.stdout)

  table = cook_table(table,options)

  out.writerows(table)


if __name__=='__main__':
  main()
