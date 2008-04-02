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
import csv
import sqlite3

from   itertools         import chain,islice

from   glu.lib.fileutils import autofile,hyphen,load_table


HEADER = ['CHROMOSOME','LOCATION','GENE NEIGHBORHOOD']

def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-c', '--column',     dest='column', default=0)
  parser.add_option('-u', '--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='the upstream margin in bases')
  parser.add_option('-d', '--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='the downstream margin in bases')
  parser.add_option('-o', '--outfile',    dest='outfile',    default='-',                metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def annotate(con,header,rows,options):
  up = options.upstream
  dn = options.downstream

  try:
    column = int(options.column)
  except ValueError:
    column = header.index(options.column)

  for row in rows:
    if not row:
      continue

    snps = query_snp(con,row[column])

    if not snps:
      yield row + ['UNKNOWN','','']
    elif len(snps) > 1:
      yield row + ['Multiple mappings','','']
    else:
      lname,chromosome,location,strand = snps[0]
      near = query_gene_neighborhood(con,chromosome,location,up,dn)
      yield row + [chromosome,location, ','.join(n[0] for n in near)]


def query_gene_neighborhood(con,chromosome,location,up,dn):
  cur = con.cursor()
  sql = '''
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '+'
    AND    featureType='GENE'
    AND    submitGroup='reference'
    AND    chromosome=?
    AND    ? BETWEEN chrStart - ? AND chrEnd + ?
  UNION
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '-'
    AND    featureType='GENE'
    AND    submitGroup='reference'
    AND    chromosome=?
    AND    ? BETWEEN chrStart - ? AND chrEnd + ?
  ORDER BY chromosome,chrStart;
  '''
  cur.execute(sql, (chromosome,location,up,dn,chromosome,location,dn,up))
  return cur.fetchall()


def query_snp(con,lname):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     SNP
  WHERE    lname = ?
  '''
  cur = con.cursor()
  cur.execute(sql,(lname,))
  return cur.fetchall()


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)<2:
    parser.print_help(sys.stderr)
    return

  con = sqlite3.connect(args[0])
  out = hyphen(options.outfile,sys.stdout)
  out = csv.writer(autofile(out,'w'),dialect='excel-tab')

  rows = load_table(hyphen(args[1],sys.stdin),want_header=True)
  header = rows.next() or ['']
  out.writerow(header + HEADER)
  out.writerows( annotate(con,header,rows,options) )


if __name__=='__main__':
  main()
