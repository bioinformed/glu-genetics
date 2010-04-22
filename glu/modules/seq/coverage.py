# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute summary statistics on SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys
import csv

from   operator    import attrgetter, itemgetter
from   itertools   import groupby
from   heapq       import heappush, heappop

import numpy as np
import pysam

from   glu.lib.fileutils            import autofile, guess_format
from   glu.modules.seq.targetfilter import read_targets


def percent(a,b):
  return a/b*100 if b else 0


def percent3(a,b):
  return a,b,percent(a,b)


def _interval_pileup_region(contig_regions):
  '''
  Generator that produces intervals with coordinates [start, end) and a
  sequence of interval names for regions with a constant number of
  overlapping intervals based on an input sequence of (potentially)
  intervals and names provided as (start, end, name).
  '''
  # Current position and window of regions
  pos    = 0
  window = []

  for read_start,read_end,name in contig_regions:
    # Yield interval prior to current read_start
    while window and window[0][0]<read_start:
      end = window[0][0]

      yield pos,end,window

      while window and window[0][0]==end:
        heappop(window)

      pos = end

    # Yield gap between last region and and current read start
    if pos!=read_start:
      yield pos,read_start,window
      pos = read_start

    # Add current region to the window
    heappush(window, (read_end,name) )

  # Drain the window once all regions have been added
  while window:
    end = window[0][0]

    # Yield current interval and reset the current location
    yield pos,end,window
    pos = end

    # Drain the window at the current endpoint
    while window and window[0][0]==end:
      heappop(window)


def interval_pileup(regions):
  '''
  Generator that produces a sequence of contig name, contig length, and a
  seriies of intervals with coordinates [start, end) and a sequence of
  interval names for regions with a constant number of overlapping intervals
  based on an input sequence of (potentially) intervals and names provided
  as (start, end, name).
  '''
  for contig_name,contig_len,contig_regions in regions:
    yield contig_name,contig_len,_interval_pileup_region(contig_regions)


def pileup_stats(regions,targets,options):
  '''
  Compute coverage depth statistics that may overlap with a series of
  disjoint target regions.
  '''
  maxcoverage  = options.maxcoverage
  coverage     = {}
  all_coverage = np.zeros( (2,maxcoverage+1), dtype=int )

  # Process each contig and match it with targets
  for contig_name,contig_len,contig_pileup in regions:
    ctargets = list(reversed(targets[contig_name]))

    assert contig_name not in coverage
    contig_coverage = coverage[contig_name] = np.zeros( (2,maxcoverage+1), dtype=int )

    #print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (contig_name,len(ctargets))

    # Process each pile-up interval within the contig
    start = end = 0
    for start,end,regions in contig_pileup:
      depth = len(regions)

      if depth > maxcoverage:
        depth = maxcoverage

      # Process all targets that overlap the current region
      while ctargets and ctargets[-1][0]<end:
        target_start,target_end = ctargets[-1]

        if start<target_start:
          contig_coverage[0,depth] += target_start-start
          start = target_start

        bound = min(target_end,end)
        contig_coverage[1,depth] += bound-start
        start = bound

        if target_end>end:
          break

        ctargets.pop()

      # Process the remaining gap between the last target and the end of region
      contig_coverage[0,depth] += end-start

    # Account for all targets beyond the last region
    for target_start,target_end in reversed(ctargets):
      contig_coverage[0,0] += max(0,target_start-end)
      contig_coverage[1,0] += target_end-max(end,target_start)
      end = target_end

    # Add the gap between the last region/target and the end of the contig
    contig_coverage[0,0] += max(0,contig_len-end)

    # Accumulate the current contig results
    all_coverage         += contig_coverage

  target_len = all_coverage.sum(axis=1)

  format = '%3d' + ' %10d %10d %6.2f%%'*3
  for depth in range(maxcoverage+1):
    print format % ((depth,) + percent3(all_coverage[1,depth],target_len[1])
                             + percent3(all_coverage[0,depth],target_len[0])
                             + percent3(all_coverage[:,depth].sum(),target_len.sum()))


def load_bam(filename,options):
  inbam = pysam.Samfile(filename, 'rb')

  try:
    contig_names = inbam.references
    contig_lens  = inbam.lengths
    reads        = inbam.fetch(region=options.region or None)

    for contig_tid,contig_reads in groupby(reads, attrgetter('rname')):
      contig_name  = contig_names[contig_tid]
      contig_len   = contig_lens[contig_tid]
      contig_reads = ( (read.pos,read.aend,read) for read in contig_reads )

      yield contig_name,contig_len,contig_reads

  finally:
    inbam.close()


def load_bed(filename,options):
  bed = csv.reader(autofile(filename),dialect='excel-tab')

  def regions(bed):
    for row in bed:
      n = len(row)
      if n<3 or row[0].startswith('track ') or row[0].startswith('#'):
        continue

      contig,start,end = row[:3]
      name = row[4] if n>4 else ''

      yield contig,int(start),int(end),name

  for contig_name,contig_regions in groupby(regions(bed), itemgetter(0)):
    contig_regions = ( (start,end,name) for contig,start,end,name in contig_regions )
    yield contig_name,0,contig_regions


def load_regions(filename,options):
  format = guess_format(filename, ['bam','bed','sam'])
  if format in ('sam','bam'):
    return load_bam(filename,options)
  elif format=='bed':
    return load_bed(filename,options)
  else:
    raise ValueError('Unknown format for input %s' % filename)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] in.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--targets', dest='targets', metavar='BED',
                    help='Single track BED file containing all targetd intervals')
  parser.add_option('--region', dest='region', metavar='REGION',
                    help='Region over which to compute as "", "contig", or "contig:start-stop".  '
                         'Default="" (all aligned reads)')
  parser.add_option('--maxcoverage', dest='maxcoverage', metavar='N', type='int', default=100,
                    help='Maximum coverage depth to track')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  if options.maxcoverage<1:
    options.maxcoverage = 1

  regions = load_regions(args[0],options)
  targets = read_targets(options.targets)
  pileup  = interval_pileup(regions)

  pileup_stats(pileup,targets,options)


if __name__=='__main__':
  main()
