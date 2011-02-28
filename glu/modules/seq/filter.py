# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Filter SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   collections           import deque
from   operator              import attrgetter, itemgetter
from   itertools             import groupby, imap

import pysam

from   glu.lib.utils         import consume
from   glu.lib.fileutils     import autofile, hyphen
from   glu.lib.progressbar   import progress_loop

from   glu.lib.seqlib.bed    import read_features
from   glu.lib.seqlib.filter import alignment_filter_options, filter_alignments


PASS,UNALIGNED,TOOSHORT,LOWOVERLAP = range(4)


def simple_filter(aligns,stats,options):
  minreadlen = options.minreadlen

  fail = options.action=='fail'

  reads_UNALIGNED  = 0
  bases_UNALIGNED  = 0
  reads_TOOSHORT   = 0
  bases_TOOSHORT   = 0
  reads_PASS       = 0
  bases_PASS       = 0

  try:
    for align in aligns:
      rlen = align.rlen

      if align.tid<0:
        reads_UNALIGNED += 1
        bases_UNALIGNED += rlen

        if fail:
          align.is_qcfail = 1
          yield align

      elif rlen<minreadlen:
        reads_TOOSHORT += 1
        bases_TOOSHORT += rlen

        if fail:
          align.is_qcfail = 1
          yield align

      else:
        reads_PASS += 1
        bases_PASS += rlen
        yield align

  finally:
    stats.reads[UNALIGNED]  += reads_UNALIGNED
    stats.bases[UNALIGNED]  += bases_UNALIGNED
    stats.reads[TOOSHORT]   += reads_TOOSHORT
    stats.bases[TOOSHORT]   += bases_TOOSHORT
    stats.reads[PASS]       += reads_PASS
    stats.bases[PASS]       += bases_PASS


def target_filter(aligns,references,targets,stats,options):
  contigs = set()

  minoverlap = options.minoverlap
  minreadlen = options.minreadlen
  nulltarget = sys.maxint,sys.maxint

  fail = options.action=='fail'

  reads_UNALIGNED  = 0
  bases_UNALIGNED  = 0
  reads_TOOSHORT   = 0
  bases_TOOSHORT   = 0
  reads_LOWOVERLAP = 0
  bases_LOWOVERLAP = 0
  reads_PASS       = 0
  bases_PASS       = 0

  try:
    for tid,contig_aligns in groupby(aligns, attrgetter('tid')):
      rname = references[tid] if tid>=0 else 'unaligned'

      ctargets = deque(targets.get(rname,[]))

      print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (rname,len(ctargets))

      assert rname not in contigs, 'Duplicate contig %s seen' % rname
      contigs.add(rname)

      if rname=='unaligned':
        for align in contig_aligns:
          reads_UNALIGNED += 1
          bases_UNALIGNED += align.rlen

          if fail:
            align.is_qcfail = True
            yield align

        continue

      next_start,next_end = ctargets[0][:2] if ctargets else nulltarget

      for align in contig_aligns:
        rlen = align.rlen

        # Fast-path 1: fail short aligns
        if rlen<minreadlen:
          reads_TOOSHORT += 1
          bases_TOOSHORT += rlen

          if fail:
            align.is_qcfail = True
            yield align

          continue

        align_end = align.aend

        # Fast-path 2: No overlap with next target (if any)
        if align_end<next_start:
          reads_LOWOVERLAP += 1
          bases_LOWOVERLAP += rlen

          if fail:
            align.is_qcfail = True
            yield align

          continue

        align_start = align.pos

        # Slow-path: Pop targets that end prior to the start of the current
        #            alignment.  This sets up the precondition no additional
        #            targets exist or that the next target begins at or before
        #            the currrent alignment start.
        while ctargets and next_end<align_start:
          ctargets.popleft()
          next_start,next_end = ctargets[0][:2] if ctargets else nulltarget

        # Compute overlap by iterating over targets until
        #   1: The current target starts after the end of the current
        #      alignment, which terminates the search.
        #   2: Otherwise, compute the amount overlap, which must be
        #      non-negative, and add it to the total overlap count.

        overlap_len = 0
        for target_start,target_end,target_name in ctargets:
          if target_start>align_end:
            break

          overlap_len += min(target_end,align_end)-max(target_start,align_start)

          # Determine if the degree of overlap is sufficient
          if overlap_len>=minoverlap:
            reads_PASS += 1
            bases_PASS += rlen
            yield align
            break

        else:
          reads_LOWOVERLAP += 1
          bases_LOWOVERLAP += rlen

          if fail:
            align.is_qcfail = True
            yield align

  finally:
    stats.reads[UNALIGNED]  += reads_UNALIGNED
    stats.bases[UNALIGNED]  += bases_UNALIGNED
    stats.reads[TOOSHORT]   += reads_TOOSHORT
    stats.bases[TOOSHORT]   += bases_TOOSHORT
    stats.reads[LOWOVERLAP] += reads_LOWOVERLAP
    stats.bases[LOWOVERLAP] += bases_LOWOVERLAP
    stats.reads[PASS]       += reads_PASS
    stats.bases[PASS]       += bases_PASS


def target_filter_generic(aligns,references,targets,stats,options):
  from   glu.lib.seqlib.intervaltree  import IntervalTree

  contigs = set()
  minoverlap = options.minoverlap

  fail  = options.action=='fail'
  reads = stats.reads
  bases = stats.bases

  for tid,contig_aligns in groupby(aligns, attrgetter('tid')):
    rname = references[tid] if tid>=0 else 'unaligned'

    contig_targets = targets.get(rname,[])

    targettree = IntervalTree()
    for target_start,target_end,target_name in contig_targets:
      targettree.insert(target_start,target_end)

    print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (rname,len(contig_targets))

    assert rname not in contigs, 'Duplicate contig %s seen' % rname
    contigs.add(rname)

    if rname=='unaligned':
      for align in contig_aligns:
        reads[UNALIGNED] += 1
        bases[UNALIGNED] += align.rlen

        if fail:
          align.is_qcfail = True
          yield align

      continue

    for align in contig_aligns:
      rlen = align.rlen

      if rlen<options.minreadlen:
        reads[TOOSHORT] += 1
        bases[TOOSHORT] += rlen

        if fail:
          align.is_qcfail = True
          yield align

        continue

      align_start = align.pos
      align_end   = align.aend

      overlap_len = 0
      for target in targettree.find(align_start, align_end):
        overlap_len += min(target.end,align_end)-max(target.start,align_start)

      if overlap_len<minoverlap:
        reads[LOWOVERLAP] += 1
        bases[LOWOVERLAP] += rlen

        if fail:
          align.is_qcfail = True
          yield align
      else:
        reads[PASS] += 1
        bases[PASS] += rlen
        yield align


class FilterStats(object):
  def __init__(self):
    self.reads = [0]*4
    self.bases = [0]*4


def sink_file(filename,header,references,lengths,aligns):
  flags  = 'wb' if filename.endswith('.bam') else 'wh'

  outbam = pysam.Samfile(filename, flags, header=header,
                                          referencenames=references,
                                          referencelengths=lengths)

  try:
    for align in aligns:
      outbam.write(align)

  finally:
    outbam.close()


def sink_null(aligns):
  consume(aligns)


def percent(a,b):
  return a/b*100 if b else 0


# UNUSED: Merging functionality was added to Pysam directly.  This code is
#         preserved for reference only.
def merge_bams(samfiles):
  from   copy                  import deepcopy
  from   glu.lib.imerge        import imerge
  from   glu.lib.utils         import chunk_iter

  header     = None
  references = None
  lengths    = None

  aligns     = []

  for samfile in samfiles:
    if references is None:
      references = samfile.references
    elif references != samfile.references:
      raise ValueError('Merging currently requires identical reference contigs')

    if lengths is None:
      lengths = samfile.lengths
    elif lengths != samfile.lengths:
      raise ValueError('Merging currently requires identical reference contigs lengths')

    if header is None:
      header = deepcopy(samfile.header)
    elif 'RG' in samfile.header:
      # FIXME: Assumes RGs are unique
      header.setdefault('RG',[]).extend(samfile.header['RG'])

    # Apply a read-ahead of 5000 alignments to minimize disk seeking at the
    # expense of a small amount of front-loaded latency
    aligns.append( chunk_iter(samfile.fetch(until_eof=True), 5000) )

  if len(aligns)>1:
    aligns = [ (((a.tid if a.tid>=0 else 99999)<<34|a.pos,a) for a     in align)
                                                             for align in aligns ]
    aligns = imap(itemgetter(1),imerge(aligns))
  else:
    aligns = aligns[0]

  return samfiles,header,references,lengths,aligns


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
  parser.add_option('-O', '--sumout', dest='sumout', metavar='FILE', default='-',
                    help='Summary output file')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    sys.exit(2)

  options.action = options.action.lower()

  if options.action not in ('drop','fail'):
    raise ValueError('Invalid filter action selected')

  samfiles   = []
  for filename in args:
    flags   = 'rb' if filename.endswith('.bam') else 'r'
    samfiles.append(pysam.Samfile(filename, flags))

  if 0: # Merging functionality is now part of our internal pysam tree
    samfiles,header,references,lengths,aligns = merge_bams(samfiles)
  else:
    merger     = pysam.SamMergeFiles(samfiles)
    header     = merger.header
    references = merger.references
    lengths    = merger.lengths
    aligns     = iter(merger)

  aligns = progress_loop(aligns, label='Loading BAM file(s): ', units='alignments')
  aligns = filter_alignments(aligns, options.includealign, options.excludealign)

  stats  = FilterStats()

  if options.targets:
    targets = read_features(options.targets)
    aligns  = target_filter(aligns,references,targets,stats,options)
  else:
    aligns  = simple_filter(aligns,stats,options)

  complete = False

  try:
    if options.output:
      sink_file(options.output,header,references,lengths,aligns)
    else:
      sink_null(aligns)

    complete = True

  except KeyboardInterrupt:
    pass

  finally:
    for samfile in samfiles:
      samfile.close()

  bases = stats.bases
  reads = stats.reads

  total_reads = sum(reads)
  total_bases = sum(bases)

  sumout = autofile(hyphen(options.sumout, sys.stdout),'w')

  if complete:
    sumout.write('Alignment summary:\n')
  else:
    sumout.write('Partial alignment summary (execution interrupted):\n')

  def stat_line(status):
    return '%12d (%6.2f%%) %17d (%6.2f%%)' % (reads[status], percent(reads[status], total_reads),
                                              bases[status], percent(bases[status], total_bases))

  sumout.write('           STATUS       READS          %            BASES          %\n')
  for status,name in enumerate(['Pass','Unaligned','Too short','Low overlap']):
    sumout.write('  %15s: %s\n' % (name, stat_line(status)))

  if not complete:
    raise KeyboardInterrupt


if __name__=='__main__':
  main()
