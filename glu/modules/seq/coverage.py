# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute summary statistics on SAM/BAM files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   operator    import attrgetter
from   itertools   import groupby

import numpy as np
import pysam

from   glu.modules.seq.targetfilter import read_targets


def percent(a,b):
  return a/b*100 if b else 0


def percent3(a,b):
  return a,b,percent(a,b)


def pileup_stats(samfile,targets,options):
  maxcoverage  = options.maxcoverage
  contig_names = samfile.references
  contig_lens  = samfile.lengths
  coverage     = {}
  all_coverage = np.zeros( (2,maxcoverage+1), dtype=int )

  pileup       = samfile.pileup(region=options.region or None)

  for tid,contig_pileup in groupby(pileup, attrgetter('tid')):
    if tid<0:
      continue

    contig_name = contig_names[tid]
    contig_len  = contig_lens[tid]
    ctargets    = list(reversed(targets[contig_name]))

    assert contig_name not in coverage
    contig_coverage = coverage[contig_name] = np.zeros( (2,maxcoverage+1), dtype=int )

    print >> sys.stderr, '[INFO] Processing contig=%s targets=%d' % (contig_name,len(ctargets))

    last = -1
    target_missed = 0
    for location in contig_pileup:
      pos   = location.pos
      depth = location.n if location.n<maxcoverage else maxcoverage

      while ctargets and ctargets[-1][1]<=pos:
        target_start,target_end = ctargets[-1]
        target_missed += target_end - max(last+1,target_start)
        ctargets.pop()

      if ctargets and ctargets[-1][0]<=pos<ctargets[-1][1]:
        target_start,target_end = ctargets[-1]
        target_missed += pos - max(last+1,target_start)
        target_covered = 1
      else:
        target_covered = 0

      last = pos

      contig_coverage[target_covered,depth] += 1

      if 0:
        quality = 0
        bases   = []
        quals   = []
        for read in location.pileups:
          if not read.is_del:
            qpos = read.qpos
            base = read.alignment.seq[qpos]
            qual = read.alignment.qual[qpos]
            quality += ord(qual)-33
          bases.append(base)
          quals.append(qual)

        bases = ''.join(bases)
        quals = ''.join(quals)

        quality /= len(bases)
        quality  = 10**(-quality/10)

      #print '%s:%-10d depth=%d, targets=%d' % (contig_name,location.pos,location.n,targets_covered)

    for target_start,target_end in ctargets:
      target_missed += target_end - max(last+1,target_start)

    contig_coverage[1,0] += target_missed
    contig_coverage[0,0]  = max(0,contig_len-contig_coverage.sum())
    all_coverage         += contig_coverage

  target_len = all_coverage.sum(axis=1)

  format = '%3d' + ' %10d %10d %6.2f%%'*3
  for depth in range(maxcoverage+1):
    print format % ((depth,) + percent3(all_coverage[1,depth],target_len[1])
                             + percent3(all_coverage[0,depth],target_len[0])
                             + percent3(all_coverage[:,depth].sum(),target_len.sum()))


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

  inbam   = pysam.Samfile(args[0], 'rb')
  targets = read_targets(options.targets)

  try:
    pileup_stats(inbam,targets,options)

  finally:
    inbam.close()


if __name__=='__main__':
  main()
