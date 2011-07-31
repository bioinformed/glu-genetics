# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Remove or mark duplicate reads in  SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import random

from   operator                     import attrgetter, itemgetter
from   itertools                    import groupby

import pysam

from   glu.lib.utils                import gcdisabled
from   glu.lib.union_find           import union_find
from   glu.lib.progressbar          import progress_loop

from   glu.modules.seq.filter       import sink_file


def read_start_generic(align):
  return align.aend if align.is_reverse else align.pos


def read_start_454(align):
  seq = align.query

  if not align.is_reverse:
    align_start = align.pos

    if seq:
      offset = len(seq)-len(seq.lstrip(seq[0]))-1
      align_start += offset

    return align_start

  else:
    align_end = align.aend

    if seq:
      offset = len(seq)-len(seq.rstrip(seq[-1]))-1
      align_end -= offset

    return align_end


def read_groups(aligns,platform):
  platform = platform.lower()

  if platform=='454':
    get_read_start = read_start_454
  elif platform in ('','neutral'):
    get_read_start = read_start_generic
  else:
    raise ValueError('Unknown platform specified: %s' % platform)

  for tid,contig_aligns in groupby(aligns,attrgetter('tid')):
    if tid==-1:
      continue

    fwd = []
    rev = []

    for align in contig_aligns:
      read_start = get_read_start(align)
      if not align.is_reverse:
        fwd.append( (read_start,align) )
      else:
        rev.append( (read_start,align) )

    fwd.sort()
    rev.sort()

    for pos,group in groupby(fwd,itemgetter(0)):
      yield [ align for pos,align in group ]

    for pos,group in groupby(rev,itemgetter(0)):
      yield [ align for pos,align in group ]


def handle_duplicates(alignments, action, duplicates):
  action = action.lower()

  if action=='keep':
    for align in alignments:
      if align.qname in duplicates:
        align.is_duplicate = True
      yield align

  elif action=='drop':
    for align in alignments:
      if align.qname not in duplicates:
        yield align

  else:
    raise ValueError('Invalid action: %s' % action)


def pick_duplicates(groups, method):
  duplicates = set()
  method     = method.lower()

  if method=='random':
    def _pick(group):
      i = random.randint(0,len(group)-1)
      group[0],group[i] = group[i],group[0]
      return group

  elif method=='best':
    def _pick(group):
      group.sort(key=lambda a: (-a.qlen,a.qname) )
      return group

  else:
    raise ValueError('Invalid primary alignment action: %s' % method)

  for group in groups:
    if len(group)>1:
      group = _pick(list(group))
      group[0].is_duplicate = False
      for g in group[1:]:
        duplicates.add(g.qname)

  return duplicates


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('bamfile', help='Input BAM file')

  parser.add_argument('--platform', metavar='NAME', default='neutral',
                    help='Sequencing platform: neutral, 454 (default=neutral).')
  parser.add_argument('--action', metavar='ACTION', default='keep',
                    help='Action to perform for duplicate reads: keep, drop.  Default=keep')
  parser.add_argument('--pick', metavar='METHOD', default='best',
                    help='Method of selecting primary alignment when keeping duplicate reads: best, random.  Default=best')
  parser.add_argument('-o', '--output', metavar='FILE',
                    help='Output BAM file')

  return parser


def main():
  parser   = option_parser()
  options  = parser.parse_args()

  inbam    = pysam.Samfile(options.bamfile,'rb')
  aligns   = inbam.fetch()
  aligns   = progress_loop(aligns, label='Loading BAM file: ', units='alignments')
  groups   = read_groups(aligns,options.platform)

  seen   = set()
  group_count = 0
  total_len   = 0

  with gcdisabled():
    uf = union_find()
    for group in groups:
      for a in group:
        if a.qname not in seen:
          seen.add(a.qname)
          total_len += a.rlen

      if len(group)>1:
        #names = [a.qname for a in group]
        uf.union(*group)
        group_count += 1

    #print len(group),[ (a.qname,a.aend if a.is_reverse else a.pos) for a in group ]

  groups      = [ g for g in uf.sets() if len(g)>1 ]
  dup_count   = sum(len(g) for g in groups) - len(groups)
  total_count = len(seen)

  if options.output!='-':
    print
    print 'Statistics for %s:'   % options.bamfile
    print '  Total  Mbps: %0.2f' % (total_len/1000000)
    print '  Total Reads: %8d'   % total_count
    print '    Dup Reads: %8d'   % dup_count
    print '  %% Dup Reads: %0.2f' % (dup_count/total_count*100)
    print

  if options.output:
    duplicates = pick_duplicates(groups, options.pick)

    inbam    = pysam.Samfile(options.bamfile,'rb')
    aligns   = inbam.fetch(until_eof=True)
    aligns   = progress_loop(aligns, label='Saving BAM file: ', units='alignments')
    aligns   = handle_duplicates(aligns, options.action, duplicates)

    sink_file(options.output, inbam, aligns)


if __name__=='__main__':
  main()
