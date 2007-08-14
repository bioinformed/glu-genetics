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

from itertools        import chain,islice

from glu.lib.fileutils    import autofile,hyphen,load_list


HEADER = ['FEATURE_NAME','CHROMOSOME','STRAND','FEATURE_START','FEATURE_END','BASES_UP',
          'BASES_DOWN','SNPS_UP','SNPS_DOWN','FEATURE_TYPE']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-u', '--upstream',   dest='upstream',   default=20000, type='int',  metavar='N',
                    help='the upstream margin in bases')
  parser.add_option('-d', '--downstream', dest='downstream', default=10000, type='int',  metavar='N',
                    help='the downstream margin in bases')
  parser.add_option('-U', '--Upstream',   dest='Upstream',                  type='int',  metavar='N',
                    help='the maximum number of upstream SNPs')
  parser.add_option('-D', '--Downstream', dest='Downstream',                type='int',  metavar='N',
                    help='the maximum number of downstream SNPs')
  parser.add_option('-o', '--outfile',    dest='outfile',    default='-',                metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def load_features(filename,limit=None):
  return islice(csv.reader(autofile(filename),dialect='excel-tab'),1,limit)


def process(con,features,options):
  for feature in features:
    feature += [None]*(5-len(feature))

    name = feature[0]

    marginup1   = marginup2   = None
    margindown1 = margindown2 = None
    chr   = feature[1] or None
    start = feature[2] or None
    end   = feature[3] or None
    strand = feature[4] or '+'

    #FIXME: handle only one part of margins is missing
    marginup1,marginup2 = options.upstream,options.Upstream
    margindown1,margindown2 = options.downstream,options.Downstream
    if start and end:
      geneinfo = query_gene(con,name)
      if not geneinfo:
        yield name,chr,strand,start,end,marginup1,margindown1,marginup2,margindown2,'REGION'
      elif (chr,int(start),int(end),strand) == geneinfo[0][2:6]:
        yield name,chr,strand,start,end,marginup1,margindown1,marginup2,margindown2,'GENE'
      else:
        print >> sys.stderr, 'Warning: conflict information for %s' % name

    elif start and not end:
      #FIXME: verify the SNP information by quering snp table
      yield name,chr,strand,start,end,marginup1,margindown1,marginup2,margindown2,'SNP'

    else:
      geneinfo = query_gene(con,name)
      if geneinfo:
        yield list(chain(geneinfo[0][1:6],[marginup1,margindown1,marginup2,margindown2,'GENE']))
        continue

      snp_info = query_snp(con,name)
      if snp_info:
        lname,chr,location,strand = snp_info[0]
        yield lname,chr,strand,location,None,None,None,None,'SNP'
        continue

      yield name,None,None,None,None,None,None,None,None,'UNKNOWN'


def query_gene(con,name):
  cur = con.cursor()
  sql = '''
  SELECT   a.Alias,s.featureName,s.chromosome,s.orientation,s.chrStart,s.chrEnd
  FROM     alias a, gene s
  WHERE    s.geneID = a.geneID
    AND    a.Alias = ?
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  cur.execute(sql, (name,))
  return cur.fetchall()


def query_snp(con,lname):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     snp
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
  out.writerow(HEADER)

  for infile in args[1:]:
    features = load_features(hyphen(infile,sys.stdin))
    results = process(con,features,options)
    out.writerows(results)


if __name__=='__main__':
  main()
