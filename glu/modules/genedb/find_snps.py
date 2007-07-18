# -*- coding: utf-8 -*-
'''
File:          find_regions.py

Authors:       Zhaoming Wang(wangzha@mail.nih.gov)
               Xiang    Deng(dengx@mail.nih.gov)

Created:       Thr Aug  3 14:45:03 EDT 2006

Abstract:      Find the snps for the region around a feature

Compatibility: Python 2.5 and above

Requires:      No external dependencies, yet...

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import csv
import sqlite3

from   itertools     import islice

from   glu.lib.utils import autofile,hyphen
from   preprocess    import load_features


HEADER = ['SNP_NAME','CHRMOSOME','LOCATION','STRAND','DISTANCE','DISTANCE_RANK',
          'REGION_START','REGION_END','FEATURE_NAME','FEATURE_STRAND','FEATURE_START','FEATURE_END','FEATURE_TYPE']


def query_snps_by_location(con,chr,start,end):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     snp
  WHERE    chromosome = ?
    AND    location BETWEEN ? AND ?
  ORDER BY location;
  '''
  cur = con.cursor()
  cur.execute(sql,(chr,start,end))
  return cur.fetchall()


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genome_database'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',  dest='outfile', default='-', metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  parser.add_option('-i', '--infile',   dest='infile',  default='-', metavar='FILE',
                    help="the name of the feature file containing list of features, '-' for standard in")
  parser.add_option('-l', '--limit',     dest='limit',   type='int',
                    help='limit the number of features being explored for testing purpose.')
  return parser


def feature_margin(start,end,strand,mup,mdown):
  if strand == '+':
    return start-mup,end+mdown
  elif strand == '-':
    return start-mdown,end+mup
  else:
    raise ValueError, 'Unknown feature orientation!'


def process_results(results,start,end,strand,nup,ndown):
  def _calc_rank_dist():
    istart = locs.index(start)
    iend = locs.index(end)
    for i,loc in enumerate(locs):
      if loc < start:
        distance = start-loc
        rank = i - istart
      elif loc > end:
        distance = loc-end
        rank = i - iend
      else:
        rank=distance=0
      if strand == '-':
        rank = -rank
      result = locdict.get(loc,None)
      if result is not None:
        result.extend([distance,rank])
        yield result

  locdict = dict( (result[2],list(result)) for result in results )
  locs= set(locdict)
  locs.add(start)
  locs.add(end)
  locs = list(locs)
  locs.sort()

  if strand == '+':
    left  = max(0,locs.index(start)-nup) if nup else 0
    right = min(len(locs),locs.index(end)+ndown+1) if ndown else len(locs)
  else:
    left  = max(0,locs.index(start)-ndown) if ndown else 0
    right = min(len(locs),locs.index(end)+nup+1) if nup else len(locs)

  locs = list(islice(locs,left,right))
  return _calc_rank_dist()


def main():
  parser = option_parser()
  options,args = parser.parse_args()
  con = sqlite3.connect(args[0])

  if len(args)<2:
    parser.print_help(sys.stderr)
    return

  out = hyphen(options.outfile,sys.stdout)
  out = csv.writer(autofile(out,'w'),dialect='excel-tab')
  out.writerow(HEADER)

  for infile in args[1:]:
    features = load_features(hyphen(infile,sys.stdin),options.limit)
    for name,chr,strand,start,end,mup,mdown,nup,ndown,featuretype in features:
      if featuretype == 'UNKNOWN':
        continue
      start = int(start)
      end = int(end or start)
      nup = nup or 0
      ndown = ndown or 0
      chrStart,chrEnd = feature_margin(start,end,strand,int(mup or 0),int(mdown or 0))
      results = query_snps_by_location(con,chr,chrStart,chrEnd)
      results = process_results(results,start,end,strand,int(nup or 0),int(ndown or 0))
      for result in results:
        result.extend([chrStart,chrEnd,name,strand,start,end,featuretype])
        out.writerow(result)


if __name__ == '__main__':
  main()
