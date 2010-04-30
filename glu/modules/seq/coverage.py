# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute summary statistics on SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys
import csv

from   collections import defaultdict
from   operator    import attrgetter, itemgetter
from   itertools   import groupby
from   heapq       import heappush, heappop

import numpy as np
import pysam

from   glu.lib.fileutils            import autofile, guess_format, table_writer
from   glu.modules.seq.targetfilter import read_targets


def percent(a,b):
  return a/b*100 if b else 0


def percent3(a,b):
  return [a,b,'%.2f' % percent(a,b)]


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
  maxcoverage     = options.maxcoverage[0]
  coverage_width  = options.maxcoverage[1]
  max_track       = maxcoverage//coverage_width
  targetout       = options.targetout
  coverage        = {}
  all_coverage    = np.zeros( (2,max_track+1), dtype=int )
  target_coverage = defaultdict(lambda: np.zeros(max_track+1, dtype=int )) if targetout else None

  # Process each contig and match it with targets
  for contig_name,contig_len,contig_pileup in regions:
    ctargets = list(reversed(targets[contig_name]))

    assert contig_name not in coverage
    contig_coverage = coverage[contig_name] = np.zeros( (2,max_track+1), dtype=int )

    #print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (contig_name,len(ctargets))

    # Process each pile-up interval within the contig
    start = end = 0
    for start,end,regions in contig_pileup:
      depth = len(regions)//coverage_width

      if depth > max_track:
        depth = max_track

      # Process all targets that overlap the current region
      while ctargets and ctargets[-1][0]<end:
        target_start,target_end,target_name = ctargets[-1]

        if start<target_start:
          contig_coverage[0,depth] += target_start-start
          start = target_start

        bound = min(target_end,end)
        contig_coverage[1,depth] += bound-start
        if targetout:
          target_coverage[target_name][depth] += bound-start
        start = bound

        if target_end>end:
          break

        ctargets.pop()

      # Process the remaining gap between the last target and the end of region
      contig_coverage[0,depth] += end-start

    # Account for all targets beyond the last region
    for target_start,target_end,target_name in reversed(ctargets):
      contig_coverage[0,0] += max(0,target_start-end)
      contig_coverage[1,0] += target_end-max(end,target_start)
      if targetout:
        target_coverage[target_name][0] += target_end-max(end,target_start)

      end = target_end

    # Add the gap between the last region/target and the end of the contig
    contig_coverage[0,0] += max(0,contig_len-end)

    # Accumulate the current contig results
    all_coverage         += contig_coverage

  return all_coverage,target_coverage


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
      name = row[3] if n>3 else ''

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


def parse_maxcoverage(options):
  try:
    maxcov = map(int,options.maxcoverage.split(':'))
    if len(maxcov)==1:
      maxcov.append(1)
    if len(maxcov)!=2:
      raise ValueError
    options.maxcoverage = tuple( max(1,m) for m in maxcov )
  except ValueError:
    raise ValueError('Invalid --maxcoverage parameter (%s)' % options.maxcoverage)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] in.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--targets', dest='targets', metavar='BED',
                    help='Single track BED file containing all targetd intervals')
  parser.add_option('--region', dest='region', metavar='REGION',
                    help='Region over which to compute as "", "contig", or "contig:start-stop".  '
                         'Default="" (all aligned reads)')
  parser.add_option('--maxcoverage', dest='maxcoverage', metavar='N:M', type='str', default='100:5',
                    help='Maximum coverage depth N to track in intervals of width M')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Overall coverage statistics')
  parser.add_option('-O', '--targetout', dest='targetout', metavar='FILE',
                    help='Per-target coverage statistics')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  parse_maxcoverage(options)

  regions = load_regions(args[0],options)
  targets = read_targets(options.targets)
  pileup  = interval_pileup(regions)

  all_coverage,target_coverage = pileup_stats(pileup,targets,options)

  maxcoverage     = options.maxcoverage[0]
  coverage_width  = options.maxcoverage[1]
  max_track       = maxcoverage//coverage_width

  if options.output:
    out = table_writer(options.output,hyphen=sys.stdout)

    out.writerow(['depth','on-target bases',  'target length',     'percent on-target',
                          'off-target bases', 'off-target length', 'percent off-target',
                          'bases',            'total length',      'percent'])

    target_len = all_coverage.sum(axis=1)
    genome_len = target_len.sum()

    for i in xrange(max_track+1):
      depth = i*coverage_width
      if i==max_track:
        depth = '>=%d' % depth
      elif coverage_width!=1:
        depth = '%d-%d' % (depth,depth+coverage_width-1)

      out.writerow([depth] + percent3(all_coverage[1,i],target_len[1])
                           + percent3(all_coverage[0,i],target_len[0])
                           + percent3(all_coverage[:,i].sum(),genome_len) )

  if options.targetout and target_coverage is not None:
    out = table_writer(options.targetout)

    header = ['target']
    if coverage_width==1:
      header += range(0,max_track)
    else:
      header += [ '%d-%d' % (i*coverage_width,(i+1)*coverage_width-1) for i in range(0,max_track)]
    header += [ '>=%d' % (max_track*coverage_width) ]

    out.writerow(header)

    for target_name in sorted(target_coverage):
      out.writerow([target_name]+target_coverage[target_name].tolist())


if __name__=='__main__':
  main()
