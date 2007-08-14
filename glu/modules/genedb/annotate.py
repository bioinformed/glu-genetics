# -*- coding: utf-8 -*-
'''
File:          preprocess.py

Authors:       Zhaoming Wang(wangzha@mail.nih.gov)
               Xiang    Deng(dengx@mail.nih.gov)
G
Created:       Tue Aug  1 14:45:03 EDT 2006

Abstract:      Generate the input file for find_snps script from SQLite database

Compatibility: Python 2.5 and above

Requires:      glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import csv
import sqlite3

from   itertools        import chain,islice

from   glu.lib.fileutils    import autofile,hyphen,load_list


HEADER = ['CHROMOSOME','LOCATION','GENE NEIGHBORHOOD']

def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-u', '--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='the upstream margin in bases')
  parser.add_option('-d', '--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='the downstream margin in bases')
  parser.add_option('-o', '--outfile',    dest='outfile',    default='-',                metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def load_features(filename,limit=None):
  return


def process(con,rows,options):
  up = options.upstream
  dn = options.downstream

  for row in rows:
    if not row:
      continue

    snps = query_snp(con,row[0])

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

  rows = csv.reader(autofile(hyphen(args[1],sys.stdin)),dialect='excel-tab')
  out.writerow(rows.next() + HEADER)
  out.writerows( process(con,rows,options) )


if __name__=='__main__':
  main()
