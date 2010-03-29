# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Filter SAM/BAM files by read length and tile overlap'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import csv

from   collections import deque,defaultdict
from   operator    import attrgetter
from   itertools   import groupby

import pysam

from   glu.lib.fileutils import autofile


GOOD,TOOSHORT,LOWOVERLAP = range(3)


def read_tiles(filename):
  bed = csv.reader(autofile(filename),dialect='excel-tab')
  track = bed.next()

  tiles = defaultdict(list)
  for contig,start,end in bed:
    tiles[contig].append( (int(start),int(end)) )

  for contig in tiles:
    tiles[contig].sort()

  return tiles


def depth_filter(reads,options):
  for read in reads:
    if read.rlen<options.minreadlen:
      yield TOOSHORT,read
    else:
      yield GOOD,read


def tile_filter(allreads,references,tiles,options):
  for contig,reads in groupby(allreads, attrgetter('rname')):
    rname  = references[contig]
    ctiles = deque(tiles[rname])

    for read in reads:
      # Fast-path: fail short reads
      if read.rlen<options.minreadlen:
        yield TOOSHORT,read
        continue

      read_start = read.pos
      read_end   = read_start+read.rlen

      #print 'Read: %s:%s (%d - %d)' % (read.rname,read.qname,read_start,read_end)

      # Remove tiles that end prior to the start of the current read
      while ctiles and ctiles[0][1]<read_start:
        #print '  Tile: (%d - %d) dropped' % (ctiles[0][0],ctiles[0][1])
        ctiles.popleft()

      overlap_len = 0
      for tile_start,tile_end in ctiles:
        if tile_start>read_end:
          break

        overlap_start,overlap_end = max(tile_start,read_start),min(tile_end,read_end)
        overlap_len += overlap_end-overlap_start
        #print '  Tile: (%d - %d) overlaps %d bases' % (tile_start,tile_end,overlap_end-overlap_start)

      #print '  %d bases overlap' % overlap_len
      if overlap_len<options.minoverlap:
        yield LOWOVERLAP,read
      else:
        yield GOOD,read


class FilterStats(object):
  def __init__(self):
    self.stats = [0]*3

  def __call__(self, reads):
    stats = self.stats
    for status,read in reads:
      stats[status] += 1
      yield status,read

  def __getitem__(self, i):
    return self.stats[i]


def action_color(reads):
  for status,read in reads:
    read.is_reverse = status==GOOD
    yield read


def action_filter(reads):
  for status,read in reads:
    if status==GOOD:
      yield read


def sink_file(filename,template,reads):
  outbam = pysam.Samfile(filename, 'wb', template=template)
  for read in reads:
    outbam.write(read)
  outbam.close()


def sink_null(reads):
  for read in reads:
    pass


def percent(a,b):
  return a/b*100 if b else 0


def option_parser():
  import optparse

  usage = 'usage: %prog [options] in.bam tiles.bed out.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--minreadlen', dest='minreadlen', metavar='N', type='int', default=85,
                    help='Minimum read length filter (default=85)')

  parser.add_option('--tiles', dest='tiles', metavar='BED',
                    help='Single track BED file containing all tiled intervals')
  parser.add_option('--minoverlap', dest='minoverlap', metavar='N', type='int', default=20,
                    help='Minimum read overlap with any tile (default=20)')

  parser.add_option('--action', dest='action', metavar='X', default='filter',
                    help='Action to perform on failing reads (filter or color, default=filter)')

  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output BAM file')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  options.action = options.action.lower()

  if options.action not in ('filter','color'):
    raise ValueError('Invalid filter action selected')

  inbam  = pysam.Samfile(args[0], 'rb')

  if options.tiles:
    tiles = read_tiles(options.tiles)
    reads = tile_filter(inbam.fetch(),inbam.references,tiles,options)
  else:
    reads = depth_filter(inbam.fetch(),options)

  stats = FilterStats()
  reads = stats(reads)

  if options.action=='filter':
    reads = action_filter(reads)
  else:
    reads = action_color(reads)

  if options.output:
    sink_file(options.output, inbam, reads)
  else:
    sink_null(reads)

  inbam.close()

  total = sum(stats)

  print 'Read summary:'
  print '         Good: %7d (%5.2f%%)' % (stats[GOOD],      percent(stats[GOOD],      total))
  print '    Too short: %7d (%5.2f%%)' % (stats[TOOSHORT],  percent(stats[TOOSHORT],  total))
  print '  Low overlap: %7d (%5.2f%%)' % (stats[LOWOVERLAP],percent(stats[LOWOVERLAP],total))


if __name__=='__main__':
  if 1:
    main()
  else:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
