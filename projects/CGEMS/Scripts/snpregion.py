__version__ = '0.1'

import os
import csv
import sys

from   itertools          import islice
from   operator           import itemgetter

from   biozilla.utils     import autofile
from   biozilla.genodata  import load_list,load_genomatrix
from   matrixsplit        import MatrixFileCache,split_fullname


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-d', '--maxdist',      dest='maxdist',      metavar='N', type='int',
                    help='Maximum distance between the reference snp and the adjacent snps')
  parser.add_option('-l', '--limit',        dest='limit',        metavar='N', type='int',
                    help='Limit the number of reference SNPs to extract')
  parser.add_option('-m', '--maxsnps',      dest='maxsnps',      metavar='N', type='int',
                    help='Maximum number of adjacent snps on either side of the reference snp')
  parser.add_option('-i', '--include',      dest='include',      metavar='FILE',
                    help='The file with the set of snps to be considered')
  parser.add_option('-r', '--reference',    dest='reference',    metavar='FILE',
                    help='The file with the set of reference snps')
  parser.add_option('-M', '--mapinfo',      dest='mapinfo',      metavar='FILE',
                    help='The snp manifest file containing chromsome and location for each snp')

  return parser


def load_snp_mapinfo(filename,rowlimit=None):
  r = csv.reader(autofile(filename),dialect='excel-tab')

  if rowlimit:
    r=islice(r,0,rowlimit)

  for row in r:
    yield row[0],row[1],int(row[2])


def filter_by_inclusion(mapinfo,include):
   for info in mapinfo:
     if info[0] in include:
       yield info


def find_adjacent_snps(refsnp,mapinfo,snpind,maxdist,maxsnps):
  if refsnp not in snpind:
    print >> sys.stderr,'   .. cannot find snp %s.  Skipping...' % refsnp
    return None

  i = snpind[refsnp]
  chr,loc = mapinfo[i][1:3]
  start,end = max(0,i-maxsnps),min(maxsnps+i+1,len(mapinfo))

  adj = []
  for xsnp,xchr,xloc in islice(mapinfo,start,end):
    if chr==xchr and abs(xloc-loc) < maxdist:
      adj.append(xsnp)

  return adj


def gen_adjacent_snps(refsnps,mapinfo,snpind,maxdist,maxsnps):
  for snp in refsnps:
    adjacent = find_adjacent_snps(snp,mapinfo,snpind,maxdist,maxsnps)
    if adjacent:
      yield snp,adjacent


def multiplexer(matrix,snpdict):
  header = matrix.next()
  header = ['key'] + header
  for rowkey,row in matrix:
    for key in snpdict.get(rowkey,[]):
      yield (key,),header,rowkey,row


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) !=1:
    parser.print_help(sys.stderr)
    return

  include = None
  if options.include:
    print >> sys.stderr, 'Loading include list'
    include = set(load_list(options.include,1))

  print >> sys.stderr, 'Loading map file'
  mapinfo = load_snp_mapinfo(options.mapinfo)

  if include is not None:
    mapinfo = filter_by_inclusion(mapinfo,include)

  mapinfo = sorted(mapinfo, key=itemgetter(1,2,0))
  snpind = dict((info[0],i) for i,info in enumerate(mapinfo))

  refsnps = load_list(options.reference)

  if options.limit:
    refsnps = islice(refsnps, options.limit)

  print >> sys.stderr, 'Finding adjacent regions'
  adjacent = gen_adjacent_snps(refsnps,mapinfo,snpind,options.maxdist,options.maxsnps)

  snpdict = {}
  for snp,adjsnps in adjacent:
    for adjsnp in adjsnps:
      snpdict.setdefault(adjsnp,[]).append(snp)

  print >> sys.stderr, 'Building region output'

  format,matrix = load_genomatrix(args[0],'ldat')

  filecache = MatrixFileCache('','ldat.gz')
  for m in multiplexer(matrix,snpdict):
    filecache.emit(*m)


if __name__=='__main__':
  main()
