# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Find SNPs near a set of genomic features (SNPs, genes or regions)'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import bisect

from   operator                        import itemgetter

from   glu.lib.fileutils               import table_reader,table_writer

from   glu.modules.genedb              import open_genedb
from   glu.modules.genedb.find_regions import resolve_features
from   glu.modules.genedb.queries      import query_snps_by_location


HEADER = ['SNP_NAME','CHRMOSOME','LOCATION','STRAND','DISTANCE','DISTANCE_RANK',
          'REGION_START','REGION_END','FEATURE_NAME','FEATURE_STRAND','FEATURE_START','FEATURE_END','FEATURE_TYPE']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-g', '--genedb',   dest='genedb', metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_option('-u', '--upbases',   dest='upbases',   default=20000, type='int',  metavar='N',
                    help='upstream margin in bases (default=20000)')
  parser.add_option('-d', '--downbases', dest='downbases', default=10000, type='int',  metavar='N',
                    help='downstream margin in bases (default=10000)')
  parser.add_option('-U', '--upsnps',    dest='upsnps',                   type='int',  metavar='N',
                    help='maximum number of upstream SNPs (default=0 for no limit)')
  parser.add_option('-D', '--downsnps',  dest='downsnps',                 type='int',  metavar='N',
                    help='maximum number of downstream SNPs (default=0 for no limit)')
  parser.add_option('-o', '--output',  dest='output', default='-', metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  return parser


def feature_margin(start,end,strand,mup,mdown):
  if strand == '+':
    return start-mup,end+mdown
  elif strand == '-':
    return start-mdown,end+mup
  else:
    raise ValueError('Unknown feature orientation')


def process_results(results,start,end,strand,nup,ndown):
  loci   = sorted(results,key=itemgetter(1,2,0))
  locs   = [ loc[2] for loc in loci ]
  m1,m2  = (nup,ndown) if strand=='+' else (ndown,nup)

  istart = bisect.bisect_left(locs,start)
  iend   = bisect.bisect_right(locs,end)

  if m2 and len(loci)-iend>m2:
    loci = loci[:iend+m2]

  if m1 and istart > m1:
    loci   = loci[istart-m1:]
    iend  -= istart-m1
    istart = m1

  def _calc_rank_dist():
    for i,loc in enumerate(loci):
      pos = loc[2]

      rank = distance = 0
      if pos < start:
        distance = pos-start
        rank     = i-istart
      elif pos > end:
        distance = pos-end
        rank     = i-iend+1

      if strand == '-':
        rank = -rank

      yield list(loc)+[distance,rank]

  return _calc_rank_dist()


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    return

  con = open_genedb(options.genedb)
  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(HEADER)

  for infile in args:
    features = table_reader(infile,want_header=True,hyphen=sys.stdin)
    features = resolve_features(con,features,options)

    for name,chr,strand,start,end,mup,mdown,nup,ndown,featuretype in features:
      if featuretype == 'UNKNOWN':
        continue
      start = int(start)
      end   = int(end or start)
      nup   = nup or 0
      ndown = ndown or 0
      chrStart,chrEnd = feature_margin(start,end,strand,int(mup or 0),int(mdown or 0))

      results = query_snps_by_location(con,chr,chrStart,chrEnd)
      results = process_results(results,start,end,strand,int(nup or 0),int(ndown or 0))

      for result in results:
        result += [chrStart,chrEnd,name,strand,start,end,featuretype]
        out.writerow(result)


if __name__ == '__main__':
  main()
