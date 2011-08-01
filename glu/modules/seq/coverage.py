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

from   glu.lib.progressbar    import progress_loop
from   glu.lib.fileutils      import autofile, guess_format, parse_augmented_name, table_writer

from   glu.lib.seqlib.bed     import read_features
from   glu.lib.seqlib.filter  import alignment_filter_options, filter_alignments


def percent(a,b):
  return a/b*100 if b else 0


def percent3(a,b):
  return [a,b,'%.2f' % percent(a,b)]


def _interval_pileup_contig(contig_regions):
  '''
  Generator that produces intervals with coordinates [start, end) and a
  sequence of interval names for regions with a constant number of
  overlapping intervals based on an input sequence of (potentially)
  intervals and names provided as (start, end, name).
  '''
  # Current position and window of regions
  pos    = 0
  window = []

  for align_start,align_end,name in contig_regions:
    if align_start is None or align_end is None:
      continue

    # Yield interval prior to current align_start
    while window and window[0][0]<align_start:
      end = window[0][0]

      yield pos,end,window

      while window and window[0][0]==end:
        heappop(window)

      pos = end

    # Yield gap between last region and and current align start
    if pos!=align_start:
      yield pos,align_start,window
      pos = align_start

    # Add current region to the window
    heappush(window, (align_end,name) )

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
  series of intervals with coordinates [start, end) and a sequence of
  interval names for regions with a constant number of overlapping intervals
  based on an input sequence of (potentially) intervals and names provided
  as (start, end, name).
  '''
  for contig_name,contig_len,contig_regions in regions:
    yield contig_name,contig_len,_interval_pileup_contig(contig_regions)


def _target_pileup_contig(contig_pileup,contig_len,targets):
  ctargets = list(reversed(targets))

  #print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (contig_name,len(ctargets))

  start = end = 0
  for start,end,regions in contig_pileup:
    depth = len(regions)

    # Process all targets that overlap the current region
    while ctargets and ctargets[-1][0]<end:
      target_start,target_end,target_name = ctargets[-1]

      if start<target_start:
        yield start,target_start,depth,None
        start = target_start

      bound = min(target_end,end)

      if start<bound:
        yield start,bound,depth,target_name

      start = bound

      if target_end>end:
        break

      ctargets.pop()

    # Process the remaining gap between the last target and the end of region
    if start<end:
      yield start,end,depth,None

  # Account for all targets beyond the last region
  for target_start,target_end,target_name in reversed(ctargets):
    if end<target_start:
      yield end,target_start,0,None

    bound = max(target_start,end)
    if bound<target_end:
      yield bound,target_end,0,target_name

    end = target_end

  # Add the gap between the last region/target and the end of the contig
  if contig_len is not None and end<contig_len:
    yield end,contig_len,0,None


def target_pileup(contig_intervals,targets):
  '''
  Compute coverage depth statistics that may overlap with a series of
  disjoint target regions.
  '''
  # Process each contig and match it with targets
  for contig_name,contig_len,contig_pileup in contig_intervals:
    yield contig_name,_target_pileup_contig(contig_pileup,contig_len,targets[contig_name])


def pileup_stats(target_intervals,options):
  '''
  Compute coverage depth statistics that may overlap with a series of
  disjoint target regions.
  '''
  max_coverage    = options.maxcoverage[0]
  coverage_width  = options.maxcoverage[1]
  max_track       = max_coverage//coverage_width
  targetout       = options.targetout
  intervalout     = options.intervalout
  coverage        = {}
  all_coverage    = np.zeros( (2,max_track+1), dtype=int )
  target_coverage = defaultdict(lambda: np.zeros(max_track+1, dtype=int )) if targetout else None

  # Process each contig and match it with targets
  for contig_name,contig_intervals in target_intervals:
    contig_coverage = coverage[contig_name] = np.zeros( (2,max_track+1), dtype=int )

    for start,end,depth,target_name in contig_intervals:
      size  = end-start
      depth = depth//coverage_width

      if depth > max_track:
        depth = max_track

      if target_name is not None:
        contig_coverage[1,depth] += size
        if targetout:
          target_coverage[contig_name,target_name][depth] += size
      else:
        contig_coverage[0,depth] += size

    # Accumulate the current contig results
    all_coverage += contig_coverage

  return all_coverage,target_coverage


def load_bam(filename,options):
  inbam = pysam.Samfile(filename, 'rb')

  try:
    contig_names = inbam.references
    contig_lens  = inbam.lengths
    aligns       = inbam.fetch(region=options.region or None)
    aligns       = progress_loop(aligns, label='Loading BAM file: ', units='alignments')
    aligns       = filter_alignments(aligns, options.includealign, options.excludealign)

    for tid,contig_aligns in groupby(aligns, attrgetter('tid')):
      contig_name   = contig_names[tid]
      contig_len    = contig_lens[tid]
      contig_aligns = ( (align.pos,align.aend,align) for align in contig_aligns )

      yield contig_name,contig_len,contig_aligns

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

  for contig_name,contig_regions in groupby(sorted(regions(bed)), itemgetter(0)):
    contig_regions = ( (start,end,name) for contig,start,end,name in contig_regions )
    yield contig_name,0,contig_regions


def load_regions(filename,options):
  format = guess_format(filename, ['bam','bed','sam'])
  name   = parse_augmented_name(filename,{})
  if format in ('sam','bam'):
    return load_bam(name,options)
  elif format=='bed':
    return load_bed(name,options)
  else:
    raise ValueError('Unknown format for input %s' % filename)


def output_summary_coverage(all_coverage, options):
  maxcoverage     = options.maxcoverage[0]
  coverage_width  = options.maxcoverage[1]
  max_track       = maxcoverage//coverage_width

  out = table_writer(options.output,hyphen=sys.stdout)

  out.writerow(['depth','on-target bases',  'target length',     'percent on-target',
                        'off-target bases', 'off-target length', 'percent off-target',
                        'bases',            'total length',      'percent'])

  target_len = all_coverage.sum(axis=1)
  genome_len = target_len.sum()

  results = []
  on_target = off_target = 0

  for i in reversed(xrange(max_track+1)):
    depth = '>=%d' % (i*coverage_width)

    on_target  += all_coverage[1,i]
    off_target += all_coverage[0,i]

    results.append([depth] + percent3(on_target,            target_len[1])
                           + percent3(off_target,           target_len[0])
                           + percent3(on_target+off_target, genome_len) )

  results.reverse()
  out.writerows(results)


def output_target_coverage(targets, target_coverage, options):
  maxcoverage     = options.maxcoverage[0]
  coverage_width  = options.maxcoverage[1]
  max_track       = maxcoverage//coverage_width

  target_dict = {}
  for target_contig,ctargets in targets.iteritems():
    for target_start,target_end,target_name in ctargets:
      target_dict[target_contig,target_name] = target_start,target_end

  out = table_writer(options.targetout)

  header = ['target','contig','target_start','target_end']
  if coverage_width==1:
    header += range(0,max_track)
  else:
    header += [ '%d-%d' % (i*coverage_width,(i+1)*coverage_width-1) for i in range(0,max_track)]
  header += [ '>=%d' % (max_track*coverage_width) ]

  out.writerow(header)

  for target_contig,target_name in sorted(target_coverage):
    target_start,target_end = target_dict[target_contig,target_name]
    coverage = target_coverage[target_contig,target_name].tolist()
    info = [target_name,target_contig,target_start,target_end]
    out.writerow(info+coverage)


def output_intervals(target_intervals, options):

  out = table_writer(options.intervalout)
  out.writerow(['contig','start','end','depth','target'])

  def _output_intervals_contig(target_intervals):
    name = [contig_name]
    for row in target_intervals:
      values = list(row)

      target_name = values[-1]
      if target_name is not None:
        values[-1] = '<on target>'
      elif not target_name:
        values[-1] = ''

      out.writerow(name+values)
      yield row

  for contig_name,target_intervals in target_intervals:
    yield contig_name,_output_intervals_contig(target_intervals)


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

  alignment_filter_options(parser)

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
  parser.add_option('--intervalout', dest='intervalout', metavar='FILE',
                    help='Output coverage by position interval')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  parse_maxcoverage(options)

  regions   = load_regions(args[0],options)
  targets   = read_features(options.targets)

  contig_intervals = interval_pileup(regions)
  target_intervals = target_pileup(contig_intervals, targets)

  if options.intervalout:
    target_intervals = output_intervals(target_intervals, options)

  all_coverage,target_coverage = pileup_stats(target_intervals, options)

  if options.output:
    output_summary_coverage(all_coverage, options)

  if options.targetout and target_coverage is not None:
    output_target_coverage(targets, target_coverage, options)


if __name__=='__main__':
  main()
