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
from   itertools             import groupby, imap, izip

import pysam

from   glu.lib.utils         import consume
from   glu.lib.fileutils     import autofile, hyphen, list_reader, table_reader, table_writer
from   glu.lib.progressbar   import progress_loop

from   glu.lib.seqlib.bed    import read_features
from   glu.lib.seqlib.filter import alignment_filter_options, filter_alignments


def get_read_group(align):
  tags = align.tags or []

  names,values = zip(*tags)
  idx = names.index('RG')
  if idx>=0:
    return values[idx]
  else:
    return None


def target_overlap(aligns,references,reference,targets,options):
  contigs     = set()
  minoverlap  = options.minoverlap
  minreadlen  = options.minreadlen
  minfoverlap = options.minfoverlap
  nulltarget  = sys.maxint,sys.maxint

  for tid,contig_aligns in groupby(aligns, attrgetter('tid')):
    rname = references[tid] if tid>=0 else 'unaligned'

    ctargets = deque(targets.get(rname,[]))

    print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (rname,len(ctargets))

    assert rname not in contigs, 'Duplicate contig %s seen' % rname
    contigs.add(rname)

    if rname=='unaligned':
      for align in contig_aligns:
        pass
      continue

    next_start,next_end = ctargets[0][:2] if ctargets else nulltarget

    for align in contig_aligns:
      rlen      = align.rlen
      align_end = align.aend

      # Fast-path 2: No overlap with next target (if any)
      if align_end<next_start:
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

      overlap     = []
      for target_start,target_end,target_name in ctargets:
        if target_start>align_end:
          break

        offby1       = target_start - align_start
        offby2       = target_end   - align_end
        offby        = min(abs(offby1),abs(offby2))

        overlap_len = min(target_end,align_end)-max(target_start,align_start)
        target_len  = target_end - target_start
        foverlap    = overlap_len/target_len

        # Determine if the degree of overlap is sufficient
        if foverlap<minfoverlap:
          continue

        if overlap_len<minoverlap:
          continue

        ref = reference.fetch(rname,align.pos,align.aend).upper()

        #print overlap_len,foverlap,target_name,get_read_group(align)
        overlap.append( (foverlap,target_name) )

        if 0:
          print align.qname

          a1,a2 = cigar_alignment(ref,align.seq,align.cigar)

          # Flip alignments if read was originally on the reverse strand
          if align.is_reverse:
            a1 = a1[::-1]
            a2 = a2[::-1]

          #a1,a2 = cigar_alignment(ref,align.query,cigar)
          print '   ref: %s' % a1
          print '  read: %s' % a2

      if overlap:
        overlap.sort()
        yield overlap[-1][1]


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('bamfile', nargs='+', help='Input BAM file(s)')

  alignment_filter_options(parser)

  parser.add_argument('--minreadlen', metavar='N', type=int, default=0,
                    help='Minimum read length filter')
  parser.add_argument('--targets', metavar='BED',
                    help='Single track BED file containing all targeted intervals')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--minoverlap', metavar='N', type=int, default=1,
                    help='Minimum alignment overlap with any target (default=1)')
  parser.add_argument('--minfoverlap', metavar='P', type=float, default=0,
                    help='Minimum fraction of target overlapping with alignment (default=0)')
  parser.add_argument('--setreadgroup', type=str, metavar='RGNAME',
                    help='Set all reads to specified read group name')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output coverage file')
  parser.add_argument('-O', '--sumout', metavar='FILE', default='-',
                    help='Summary output file')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  reference = pysam.Fastafile(options.reference)
  targets   = read_features(options.targets,merge=False)
  target_names = []

  for contig in targets:
    for target in targets[contig]:
      target_names.append(target[2])

  target_names.sort()

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['FILENAME']+target_names)

  samfiles = []
  for filename in options.bamfile:
    flags   = 'rb' if filename.endswith('.bam') else 'r'
    samfile = pysam.Samfile(filename, flags)

    header     = samfile.header
    references = samfile.references
    lengths    = samfile.lengths
    aligns     = iter(samfile)

    aligns = progress_loop(aligns, label='Loading BAM file(s): ', units='alignments')
    aligns = filter_alignments(aligns, options.includealign, options.excludealign)


    target_stats = {}
    for target_name in target_names:
      target_stats[target_name] = 0

    for target_name in target_overlap(aligns,references,reference,targets,options):
      target_stats[target_name] += 1

    target_stats = [ target_stats.get(target_name,0) for target_name in target_names ]
    out.writerow([filename]+target_stats)


if __name__=='__main__':
  main()
