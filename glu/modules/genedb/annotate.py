# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Add columns of genomic annotation to a file containing a list of SNP names'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.fileutils          import load_table,table_writer,resolve_column_headers

from   glu.modules.genedb         import open_genedb
from   glu.modules.genedb.queries import query_snps_by_name,query_gene_neighborhood


HEADER = ['CHROMOSOME','LOCATION','GENE NEIGHBORHOOD']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-g', '--genedb',   dest='genedb', metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_option('-c', '--column',     dest='column',     default=0,
                    help='Column name or number in which to find SNPs')
  parser.add_option('-u', '--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='upstream margin in bases (default=20000)')
  parser.add_option('-d', '--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='the downstream margin in bases (default=10000)')
  parser.add_option('-o', '--outfile',    dest='outfile',    default='-',                metavar='FILE',
                    help="name of the output file, '-' for standard out")
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

    snps = query_snps_by_name(con,row[column])

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

  if len(args)!=1:
    parser.print_help(sys.stderr)
    return

  con    = open_genedb(options.genedb)
  rows   = load_table(args[0],want_header=True,hyphen=sys.stdin)
  out    = table_writer(options.outfile,hyphen=sys.stdout)

  header = rows.next() or ['']

  out.writerow(header + HEADER)
  out.writerows( annotate(con,header,rows,options) )


if __name__=='__main__':
  main()
