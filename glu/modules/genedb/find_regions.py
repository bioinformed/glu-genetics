# -*- coding: utf-8 -*-
'''
File:          find_regions.py

Authors:       Zhaoming Wang(wangzha@mail.nih.gov)
               Xiang    Deng(dengx@mail.nih.gov)

Created:       Tue Aug  1 14:45:03 EDT 2006

Abstract:      Resolve genomic metadata given feature names for SNPs, genes,
               and bounded regions

Compatibility: Python 2.5 and above

Requires:      glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import sqlite3

from   glu.lib.fileutils          import load_table,table_writer,tryint

from   glu.modules.genedb.queries import query_gene, query_snp


HEADER = ['FEATURE_NAME','CHROMOSOME','STRAND','FEATURE_START','FEATURE_END','BASES_UP',
          'BASES_DOWN','SNPS_UP','SNPS_DOWN','FEATURE_TYPE']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-u', '--upbases',   dest='upbases',   default=20000, type='int',  metavar='N',
                    help='upstream margin in bases')
  parser.add_option('-d', '--downbases', dest='downbases', default=10000, type='int',  metavar='N',
                    help='downstream margin in bases')
  parser.add_option('-U', '--upsnps',    dest='upsnps',                   type='int',  metavar='N',
                    help='maximum number of upstream SNPs')
  parser.add_option('-D', '--downsnps',  dest='downsnps',                 type='int',  metavar='N',
                    help='maximum number of downstream SNPs')
  parser.add_option('-o', '--outfile',   dest='outfile',   default='-',                metavar='FILE',
                    help="output file name, '-' for standard out")
  return parser


def coalesce(*items):
  for i in items:
    if i not in ('',None):
      return i
  return None


def is_int(i):
  return isinstance(i, (int,long))


def resolve_feature(con,feature,options):
  name        = feature[0]
  chr         = feature[1] or None
  strand      = feature[2] or '+'
  start       = coalesce(tryint(feature[3]))
  end         = coalesce(tryint(feature[4]))
  upbases     = coalesce(tryint(feature[5]), options.upbases)
  downbases   = coalesce(tryint(feature[6]), options.downbases)
  upsnps      = coalesce(tryint(feature[7]), options.upsnps)
  downsnps    = coalesce(tryint(feature[8]), options.downsnps)

  if (start is not None and end is None) or (is_int(start) and is_int(end) and start+1==end):
    feature = 'SNP'
  elif chr and start and end:
    geneinfo = query_gene(con,name)
    if any( (chr,start,end) == (gi[2],gi[4],gi[5]) for gi in geneinfo):
      feature = 'GENE'
    else:
      feature = 'REGION'
  else:
    geneinfo = query_gene(con,name)
    if len(geneinfo) == 1:
      name,chr,strand,start,end = geneinfo[0][1:6]
      feature = geneinfo[0][6]
    elif len(geneinfo) > 1:
      feature = 'AMBIGUOUS'
    else:
      snpinfo = query_snp(con,name)
      if len(snpinfo)==1:
        lname,chr,start,strand = snpinfo[0]
        end = start+1
        feature = 'SNP'
      elif len(snpinfo) > 1:
        feature = 'AMBIGUOUS'
      else:
        feature = 'UNKNOWN'

  return name,chr,strand,start,end,upbases,downbases,upsnps,downsnps,feature


def resolve_features(con,features,options):
  for feature in features:
    feature += [None]*(10-len(feature))
    yield resolve_feature(con,feature,options)


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)<2:
    parser.print_help(sys.stderr)
    return

  con = sqlite3.connect(args[0])
  out = table_writer(options.outfile,hyphen=sys.stdout)
  out.writerow(HEADER)

  for infile in args[1:]:
    features = load_table(infile,want_header=True,hyphen=sys.stdin)
    results  = resolve_features(con,features,options)
    out.writerows(results)


if __name__=='__main__':
  main()
