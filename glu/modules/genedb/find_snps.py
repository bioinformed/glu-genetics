# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Find SNPs near a set of genomic features (SNPs, genes or regions)'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import bisect

from   operator                        import itemgetter

from   glu.lib.fileutils               import list_reader,table_reader,table_writer

from   glu.lib.genedb                  import open_genedb
from   glu.lib.genedb.queries          import query_snps_by_location

from   glu.modules.genedb.find_regions import resolve_features


HEADER = ['LOCUS','CHROMOSOME','START','END','STRAND','REF_ALLELE','OBSERVED_ALLELES','dbSNP_VCLASS','dbSNP_FUNC',
          'DISTANCE','DISTANCE_RANK','REGION_START','REGION_END',
          'FEATURE_NAME','FEATURE_STRAND','FEATURE_START','FEATURE_END','FEATURE_TYPE']


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('table', nargs='+', help='Tabular or delimited input file')

  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('--includeloci', metavar='FILE',
                    help='List of loci to include')
  parser.add_argument('--excludeloci', metavar='FILE',
                    help='List of loci to exclude')
  parser.add_argument('-u', '--upbases',   default=20000, type=int,  metavar='N',
                    help='upstream margin in bases (default=20000)')
  parser.add_argument('-d', '--downbases', default=10000, type=int,  metavar='N',
                    help='downstream margin in bases (default=10000)')
  parser.add_argument('-U', '--upsnps',    type=int,  metavar='N',
                    help='maximum number of upstream SNPs (default=0 for no limit)')
  parser.add_argument('-D', '--downsnps',  type=int,  metavar='N',
                    help='maximum number of downstream SNPs (default=0 for no limit)')
  parser.add_argument('-F', '--outformat', default='GLU',              metavar='NAME',
                    help='Output format (GLU or BED)')
  parser.add_argument('-o', '--output',  default='-', metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def feature_margin(start,end,strand,mup,mdown):
  if strand == '+':
    return start-mup,end+mdown
  elif strand == '-':
    return start-mdown,end+mup
  else:
    raise ValueError('Unknown feature orientation')


def annotate_results(results,start,end,strand,nup,ndown):
  # Sort by chrom, start, end, name
  m1,m2  = (nup,ndown) if strand=='+' else (ndown,nup)
  loci   = sorted(results,key=itemgetter(1,2,3,0))
  starts = [ loc[2] for loc in loci ]

  istart = bisect.bisect_left(starts,  start)
  iend   = bisect.bisect_right(starts, end)

  if m2 and len(loci)-iend>m2:
    loci = loci[:iend+m2]

  if m1 and istart > m1:
    loci   = loci[istart-m1:]
    iend  -= istart-m1
    istart = m1

  def _calc_rank_dist():
    for i,loc in enumerate(loci):
      var_start,var_end = loc[2:4]

      rank = distance = 0
      if var_end < start:
        distance = var_end-start
        rank     = i-istart
      elif var_start > end:
        distance = var_start-end
        rank     = i-iend+1

      if strand == '-':
        rank = -rank

      yield list(loc)+[distance,rank]

  return _calc_rank_dist()


def as_set(f):
  '''
  Return f, based on the following rules.  If f is:
    1. None,  return None
    2. a set, return f
    3. a dict, return set(f)
    4. otherwise, pass f through a list_reader and
       return the results as a set

  >>> as_set(None) is None
  True
  >>> as_set(set('abc')) == set('abc')
  True
  >>> as_set([1, 2, 3]) == set([1,2,3])
  True
  >>> as_set(iter(['a\\t1','b\\t2','c\\t3'])) == set('abc')
  True
  '''
  if f is None:
    return None
  elif isinstance(f, set):
    return f
  elif isinstance(f, (dict,list,tuple)):
    return set(f)
  else:
    return set(list_reader(f))


class ResultFilter(object):
  '''
  Filter SNPs based on an inclusion or exclusion list
  '''
  def __init__(self,options):
    include = as_set(options.includeloci)
    exclude = as_set(options.excludeloci)

    if include is not None and exclude is not None:
      include = include - exclude
      exclude = None

    self.include = include
    self.exclude = exclude

  def filter(self,results):
    if self.include is None and self.exclude is None:
      return results

    if self.include is not None:
      def _filter():
        include = self.include
        for result in results:
          if result[0] in include:
            yield result
    else:
      def _filter():
        exclude = self.exclude
        for result in results:
          if result[0] not in exclude:
            yield result

    return _filter()


def main():
  parser  = option_parser()
  options = parser.parse_args()

  options.outformat = options.outformat.lower()

  if options.outformat not in ('glu','bed'):
    raise ValueError('Unknown output format selected: %s' % options.outformat)

  con = open_genedb(options.genedb)
  out = table_writer(options.output,hyphen=sys.stdout)

  if options.outformat not in ('glu','bed'):
    out.writerow(HEADER)

  filter_results = ResultFilter(options).filter

  for filename in options.table:
    features = table_reader(filename,want_header=True,hyphen=sys.stdin)
    features = resolve_features(con,features,options)

    for name,chrom,bands,strand,start,end,mup,mdown,nup,ndown,featuretype in features:
      if featuretype == 'UNKNOWN' or start is None:
        continue

      start = int(start)
      end   = int(end or start)
      nup   = nup or 0
      ndown = ndown or 0
      chromStart,chromEnd = feature_margin(start,end,strand,int(mup or 0),int(mdown or 0))

      results = query_snps_by_location(con,chrom,chromStart,chromEnd)
      results = [ r[:-1] for r in results ]
      results = filter_results(results)
      results = annotate_results(results,start,end,strand,int(nup or 0),int(ndown or 0))

      if options.outformat=='glu':
        for result in results:
          row = [ result[0], result[1], result[2]+1, result[3], result[4], result[5],
                  result[6], result[7], result[8], chromStart, chromEnd, name, strand,
                  start, end, featuretype ]
          out.writerow(row)
      else:
        for result in results:
          row = [ result[1], result[2], result[3], result[0] ]
          out.writerow(row)


if __name__ == '__main__':
  main()
