# -*- coding: utf-8 -*-
'''
File:          preprocess.py

Authors:       Zhaoming Wang(wangzha@mail.nih.gov)
               Xiang    Deng(dengx@mail.nih.gov)

Created:       Tue Aug  1 14:45:03 EDT 2006

Abstract:      Generate the input file for find_snps script from SQLite database

Compatibility: Python 2.5 and above

Requires:      glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import sqlite3

from   glu.lib.fileutils          import hyphen,load_table,table_writer,resolve_column_headers

from   glu.modules.genedb.queries import query_snp,query_gene_neighborhood


HEADER = ['CHROMOSOME','LOCATION','GENE NEIGHBORHOOD']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-c', '--column',     dest='column',     default=0)
  parser.add_option('-u', '--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='the upstream margin in bases')
  parser.add_option('-d', '--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='the downstream margin in bases')
  parser.add_option('-o', '--outfile',    dest='outfile',    default='-',                metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def annotate(con,header,rows,options):
  up = options.upstream   or 0
  dn = options.downstream or 0

  column = resolve_column_headers(header,[options.column])

  if len(column) != 1:
    raise ValueError('Invalid SNP column')

  column = column[0]

  for row in rows:
    if len(row) <= column:
      yield row + ['']*3
      continue

    snps = query_snp(con,row[column])

    if not snps:
      info = ['UNKNOWN','','']
    elif len(snps) > 1:
      info = ['Multiple mappings','','']
    else:
      lname,chromosome,location,strand = snps[0]
      near = query_gene_neighborhood(con,chromosome,location,up,dn)
      info = [chromosome,location, ','.join(n[0] for n in near)]

    yield row + info


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)<2:
    parser.print_help(sys.stderr)
    return

  con    = sqlite3.connect(args[0])
  rows   = load_table(args[1],want_header=True,hyphen=sys.stdin)
  out    = table_writer(options.outfile,hyphen=sys.stdout)

  header = rows.next() or ['']

  out.writerow(header + HEADER)
  out.writerows( annotate(con,header,rows,options) )


if __name__=='__main__':
  main()
