# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Filter SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import csv

from   collections           import deque, defaultdict
from   operator              import attrgetter
from   itertools             import groupby

import pysam

from   glu.lib.utils         import consume

from   glu.lib.seqlib.bed    import read_features
from   glu.lib.seqlib.filter import alignment_filter_options, filter_alignments


PASS,UNALIGNED,TOOSHORT,LOWOVERLAP = range(4)


def depth_filter(aligns,options):
  for align in aligns:
    if align.rlen<options.minreadlen:
      yield TOOSHORT,align
    else:
      yield PASS,align


def target_filter(allaligns,references,targets,options):
  contigs = set()

  for contig,aligns in groupby(allaligns, attrgetter('rname')):
    rname = references[contig] if contig>= 0 else 'unaligned'

    ctargets = deque(targets[rname])
    print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (rname,len(ctargets))

    assert rname not in contigs, 'Duplicate contig %s seen' % rname
    contigs.add(rname)

    if rname=='unaligned':
      for align in aligns:
        yield UNALIGNED,align
      continue

    last = -1,-1

    for align in aligns:
      align_start = align.pos
      align_end   = align.aend

      here = align_start,align_start
      assert last<=here, 'Out of order align (%s>%s)' % (last,here)
      last = here

      # Fast-path: fail short aligns
      if align.rlen<options.minreadlen:
        yield TOOSHORT,align
        continue

      # Remove targets that end prior to the start of the current align
      while ctargets and ctargets[0][1]<align_start:
        ctargets.popleft()

      overlap_len = 0
      for target_start,target_end,target_name in ctargets:
        if target_start>align_end:
          break
        if align_start>=target_end:
          continue
        overlap_start,overlap_end = max(target_start,align_start),min(target_end,align_end)
        overlap_len += overlap_end-overlap_start

      if overlap_len<options.minoverlap:
        yield LOWOVERLAP,align
      else:
        yield PASS,align


class FilterStats(object):
  def __init__(self):
    self.stats = [0]*4

  def __call__(self, aligns):
    stats = self.stats
    for status,align in aligns:
      stats[status] += 1
      yield status,align

  def __getitem__(self, i):
    return self.stats[i]


def action_fail(aligns):
  for status,align in aligns:
    align.is_qcfail = status!=PASS
    yield align


def action_drop(aligns):
  for status,align in aligns:
    if status==PASS:
      yield align


def sink_file(filename,template,aligns):
  flags  = 'wb' if filename.endswith('.bam') else 'wh'
  outbam = pysam.Samfile(filename, flags, template=template)

  try:
    for align in aligns:
      outbam.write(align)
  finally:
    outbam.close()


def sink_null(aligns):
  consume(aligns)


def percent(a,b):
  return a/b*100 if b else 0


def option_parser():
  import optparse

  usage = 'usage: %prog [options] in.bam'
  parser = optparse.OptionParser(usage=usage)

  alignment_filter_options(parser)

  parser.add_option('--minreadlen', dest='minreadlen', metavar='N', type='int', default=0,
                    help='Minimum read length filter')
  parser.add_option('--targets', dest='targets', metavar='BED',
                    help='Single track BED file containing all targeted intervals')
  parser.add_option('--minoverlap', dest='minoverlap', metavar='N', type='int', default=1,
                    help='Minimum alignment overlap with any target (default=1)')

  parser.add_option('--action', dest='action', metavar='X', default='fail',
                    help='Action to perform on failing alignments (drop or fail, default=fail)')

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

  if options.action not in ('drop','fail'):
    raise ValueError('Invalid filter action selected')

  flags  = 'rb' if args[0].endswith('.bam') else 'r'
  inbam  = pysam.Samfile(args[0], flags)
  aligns = inbam.fetch(until_eof=True)
  aligns = filter_alignments(aligns, options.includealign, options.excludealign)

  if options.targets:
    targets = read_features(options.targets)
    aligns  = target_filter(aligns,inbam.references,targets,options)
  else:
    aligns  = depth_filter(aligns,options)

  stats  = FilterStats()
  aligns = stats(aligns)

  if options.action=='drop':
    aligns = action_drop(aligns)
  else:
    aligns = action_fail(aligns)

  complete = False

  try:
    if options.output:
      sink_file(options.output, inbam, aligns)
    else:
      sink_null(aligns)

    complete = True

  except KeyboardInterrupt:
    pass

  finally:
    inbam.close()

  total = sum(stats)

  if complete:
    print 'Alignment summary:'
  else:
    print 'Partial alignment summary (execution interrupted):'

  print '         Pass: %10d (%5.2f%%)' % (stats[PASS],      percent(stats[PASS],      total))
  print '    Unaligned: %10d (%5.2f%%)' % (stats[UNALIGNED], percent(stats[UNALIGNED], total))
  print '    Too short: %10d (%5.2f%%)' % (stats[TOOSHORT],  percent(stats[TOOSHORT],  total))
  print '  Low overlap: %10d (%5.2f%%)' % (stats[LOWOVERLAP],percent(stats[LOWOVERLAP],total))
  print

  if not complete:
    raise KeyboardInterrupt


if __name__=='__main__':
  main()
