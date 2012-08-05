# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert a VCF file to an ldat file'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

import numpy as np

from   itertools          import izip

from   glu.lib.fileutils  import table_writer
from   glu.lib.seqlib.vcf import VCFReader


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input VCF variant file')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser   = option_parser()
  options  = parser.parse_args()

  vcf      = VCFReader(options.variants,hyphen=sys.stdin)
  out      = table_writer(options.output,hyphen=sys.stdout)

  samples  = vcf.samples
  n        = len(samples)

  genomaps = {}

  out.writerow(['ldat']+samples)

  seen     = set()
  non_autosomes = set(['chrX','chrY','chrM','X','Y','M','Mt','MT'])

  for v in vcf:
    if v.chrom in non_autosomes or not v.names or len(v.var)!=1:
      continue

    a,b = ab = v.ref,v.var[0]

    # Count SNPs only for now
    if len(a)!=1 or len(b)!=1:
      continue

    try:
      name = next(name for name in v.names if name.startswith('rs'))
    except StopIteration:
      continue

    if name in seen:
      continue

    seen.add(name)

    genomap = genomaps.get(ab)

    if genomap is None:
      genomap = genomaps[ab] = {'0/0':a+a,'0/1':a+b,'1/0':a+b,'1/1':b+b,
                                '0|0':a+a,'0|1':a+b,'1|0':a+b,'1|1':b+b}

    genos = [ genomap.get(g[0],'  ') for g in v.genos ]

    if genos.count('  ')>n/2:
      continue

    out.writerow( [name] + genos )


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
      stats = pstats.Stats(prof,stream=sys.stderr)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
