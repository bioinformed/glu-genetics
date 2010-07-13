# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Filter SAM/BAM files'
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


GOOD,UNALIGNED,TOOSHORT,LOWOVERLAP = range(4)


def merge_targets(targets):
  targets.sort()

  current_start,current_end,current_name = None,None,set()

  for target_start,target_end,target_name in targets:
    if current_start is None:
      current_start,current_end,current_name = target_start,target_end,set([target_name])
    elif current_end<target_start:
      yield current_start,current_end,','.join(sorted(current_name))
      current_start,current_end,current_name = target_start,target_end,set([target_name])
    else:
      current_name.add(target_name)
      current_end = max(current_end,target_end)

  if current_start is not None:
    yield current_start,current_end,','.join(sorted(current_name))


def read_targets(filename):
  targets = defaultdict(list)

  if not filename:
    return targets

  bed = csv.reader(autofile(filename),dialect='excel-tab')

  i = 1
  for row in bed:
    n = len(row)
    if n<3 or row[0].startswith('track ') or row[0].startswith('#'):
      continue

    contig,start,end = row[:3]

    if n>3:
      name = row[3]
    else:
      name = 'target_%06d' % i
      i += 1

    targets[contig].append( (int(start),int(end),name) )

  for contig in targets:
    targets[contig] = list(merge_targets(targets[contig]))

  return targets


def depth_filter(reads,options):
  for read in reads:
    if read.rlen<options.minreadlen:
      yield TOOSHORT,read
    else:
      yield GOOD,read


def target_filter(allreads,references,targets,options):
  contigs = set()

  for contig,reads in groupby(allreads, attrgetter('rname')):
    rname = references[contig] if contig>= 0 else 'unaligned'

    ctargets = deque(targets[rname])
    print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (rname,len(ctargets))

    assert rname not in contigs, 'Duplicate contig %s seen' % rname
    contigs.add(rname)

    if rname=='unaligned':
      for read in reads:
        yield UNALIGNED,read
      continue

    last = -1,-1

    for read in reads:
      read_start = read.pos
      read_end   = read.aend

      here = read_start,read_start
      assert last<=here, 'Out of order read (%s>%s)' % (last,here)
      last = here

      # Fast-path: fail short reads
      if read.rlen<options.minreadlen:
        yield TOOSHORT,read
        continue

      # Remove targets that end prior to the start of the current read
      while ctargets and ctargets[0][1]<read_start:
        ctargets.popleft()

      overlap_len = 0
      for target_start,target_end,target_name in ctargets:
        if target_start>read_end:
          break
        if read_start>=target_end:
          continue
        overlap_start,overlap_end = max(target_start,read_start),min(target_end,read_end)
        overlap_len += overlap_end-overlap_start

      if overlap_len<options.minoverlap:
        yield LOWOVERLAP,read
      else:
        yield GOOD,read


class FilterStats(object):
  def __init__(self):
    self.stats = [0]*4

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

  try:
    for read in reads:
      outbam.write(read)
  finally:
    outbam.close()


def sink_null(reads):
  for read in reads:
    pass


def percent(a,b):
  return a/b*100 if b else 0


def option_parser():
  import optparse

  usage = 'usage: %prog [options] in.bam'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--minreadlen', dest='minreadlen', metavar='N', type='int', default=85,
                    help='Minimum read length filter (default=85)')

  parser.add_option('--targets', dest='targets', metavar='BED',
                    help='Single track BED file containing all targeted intervals')
  parser.add_option('--minoverlap', dest='minoverlap', metavar='N', type='int', default=20,
                    help='Minimum read overlap with any target (default=20)')

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
  reads  = inbam.fetch(until_eof=True)

  if options.targets:
    targets = read_targets(options.targets)
    reads = target_filter(reads,inbam.references,targets,options)
  else:
    reads = depth_filter(reads,options)

  stats = FilterStats()
  reads = stats(reads)

  if options.action=='filter':
    reads = action_filter(reads)
  else:
    reads = action_color(reads)

  complete = False

  try:
    if options.output:
      sink_file(options.output, inbam, reads)
    else:
      sink_null(reads)

    complete = True

  except KeyboardInterrupt:
    pass

  finally:
    inbam.close()

  total = sum(stats)

  if complete:
    print 'Read summary:'
  else:
    print 'Partial read summary (execution interrupted):'

  print '         Good: %10d (%5.2f%%)' % (stats[GOOD],      percent(stats[GOOD],      total))
  print '    Unaligned: %10d (%5.2f%%)' % (stats[UNALIGNED], percent(stats[UNALIGNED], total))
  print '    Too short: %10d (%5.2f%%)' % (stats[TOOSHORT],  percent(stats[TOOSHORT],  total))
  print '  Low overlap: %10d (%5.2f%%)' % (stats[LOWOVERLAP],percent(stats[LOWOVERLAP],total))
  print

  if not complete:
    raise KeyboardInterrupt


if __name__=='__main__':
  main()
