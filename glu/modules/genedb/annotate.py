# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Add columns of genomic annotation to a file containing a list of SNP names'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.fileutils      import table_reader,table_writer,resolve_column_headers

from   glu.lib.genedb         import open_genedb
from   glu.lib.genedb.queries import query_snps_by_name,query_gene_neighborhood,query_cytoband_by_location


HEADER = ['CHROMOSOME','CYTOBAND','START','END','GENE NEIGHBORHOOD','dbSNP ANNOTATION']


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('table', help='Tabular or delimited input file')

  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('-c', '--column',     default=0,
                    help='Column name or number in which to find SNPs')
  parser.add_argument('-u', '--upstream',   default=20000, type=int,  metavar='N',
                    help='upstream margin in bases (default=20000)')
  parser.add_argument('-d', '--downstream', default=10000, type=int,  metavar='N',
                    help='the downstream margin in bases (default=10000)')
  parser.add_argument('-o', '--output',    default='-',                metavar='FILE',
                    help="name of the output file, '-' for standard out")
  return parser


def normalize(row,n):
  m = len(row)
  if m>n:
    row = row[:n]
  elif m<n:
    row += ['']*(n-m)
  return row


def annotate(con,header,rows,options):
  up = options.upstream   or 0
  dn = options.downstream or 0

  column = resolve_column_headers(header,[options.column])

  if len(column) != 1:
    raise ValueError('Invalid SNP column')

  column = column[0]
  n = len(header)

  for row in rows:
    if row.count('') == len(row):
      continue

    row = normalize(row,n)
    snp = row[column]
    results = query_snps_by_name(con,snp) if snp else None

    if results is None:
      info = ['']*3
    elif not results:
      info = ['UNKNOWN','','']
    elif len(results) > 1:
      info = ['Multiple mappings','','']
    else:
      name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight = results[0]
      near     = query_gene_neighborhood(con,chrom,start,end,up,dn)
      cytoband = query_cytoband_by_location(con,chrom,start)[0]

      info = [chrom,
              cytoband,
              start+1, end,
              ','.join(n[0] for n in near),
              func]

    yield row + info


def main():
  parser  = option_parser()
  options = parser.parse_args()
  con     = open_genedb(options.genedb)
  rows    = table_reader(options.table,want_header=True,hyphen=sys.stdin)
  out     = table_writer(options.output,hyphen=sys.stdout)

  try:
    header = rows.next()
  except StopIteration:
    header = []

  header = header or ['']

  out.writerow(header + HEADER)
  out.writerows( annotate(con,header,rows,options) )


if __name__=='__main__':
  main()
