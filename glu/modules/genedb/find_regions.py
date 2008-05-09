# -*- coding: utf-8 -*-
'''
File:          find_regions.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:      Front-end for running TagZilla for SNPPlex and Illumina assay design

Compatibility: Python 2.5 and above

Requires:      No external dependencies, yet...

Revision:      $Id$
'''

__program__   = 'find_regions'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys
import sqlite3

from   glu.lib.fileutils import load_table,table_writer

HAPMAP_PATH= '/home/jacobske/projects/CGEMS/hapmap'
GENOTYPES  = os.path.join(HAPMAP_PATH,'build20','non-redundant')
PEDIGREES  = os.path.join(HAPMAP_PATH,'pedigrees','peds')
POPS       = ['CEU','YRI','JPT','CHB','JPT+CHB']

def query_genes_by_name(con,gene):
  sql = '''
  SELECT   a.Alias,s.featureName,s.chromosome,s.chrStart,s.chrEnd,s.orientation
  FROM     alias a, gene s
  WHERE    s.geneID = a.geneID
    AND    a.Alias %s
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  if '%' in gene:
    sql = sql % 'LIKE ?'
  else:
    sql = sql % '= ?'

  cur = con.cursor()
  cur.execute(sql, [gene])
  return cur.fetchall()


def query_gene_by_name(con,gene):
  genes = query_genes_by_name(con,gene)
  if not genes:
    raise KeyError('Cannot find gene "%s"' % gene)
  elif len(genes) > 1:
    raise KeyError('Gene not unique "%s"' % gene)
  return genes[0]


def query_genes_by_location(con,chr,loc):
  sql = '''
  SELECT   s.featureName,s.featureName,s.chromosome,s.chrStart,s.chrEnd,s.orientation
  FROM     gene s
  WHERE    s.chromosome = ?
    AND    ? BETWEEN s.chrStart AND s.chrEnd
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, chr, loc)
  return cur.fetchall()


def gene_margin(gene, upstream=20000, downstream=10000):
  if None in tuple(gene[2:5]):
    return None,None
  if gene[5] == '+':
    return gene[-3]-upstream,gene[-2]+downstream
  elif gene[5] == '-':
    return gene[-3]-downstream,gene[-2]+upstream
  else:
    raise ValueError('Unknown gene orientation for %s=%s' % (gene[0],gene[5]))


def get_snps(con, chromosome, start, stop):
  sql = '''
  SELECT   lname, chromosome, location
  FROM     snp
  WHERE    chromosome=?
    AND    location BETWEEN ? AND ?;
  '''
  cur = con.cursor()
  cur.execute(sql,[chromosome,start,stop])
  return cur.fetchall()


def escape(s):
  return "'%s'" % s.replace("'","''")


def extend(s,n):
  m = len(s)
  if m < n:
    s = list(s)
    s.extend(['']*(n-m))
  return s


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome genofile...'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                          help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  return parser


def find_regions(options,args):
  con = sqlite3.connect(args[0])

  out = table_writer(sys.stdout)
  out.writerow( ['GENE','CHROMOSOME','GENE START','GENE END','LOCUS','LOCATION'] )

  for arg in args[1:]:
    genefile = load_table(arg,want_header=True)
    header= genefile.next()
    genes = [ line[0] for line in genefile if line and line[0] ]

    for gene in genes:
      try:
        row = query_gene_by_name(con, gene)
      except KeyError,e:
        print >> sys.stderr, 'Invalid gene: %s' % e.args[0]
        continue

      if gene != row[0]:
        gene = '%s (%s)' % (gene,row[0])

      chromosome = row[2]
      start,stop = gene_margin(row,20000,10000)
      snps = get_snps(con, chromosome, start, stop)

      for lname,chromosome,location in snps:
        out.writerow( [gene,chromosome,row[3],row[4],lname,location] )


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)<2:
    parser.print_help(sys.stderr)
    return

  find_regions(options,args)


if __name__ == '__main__':
  main()
