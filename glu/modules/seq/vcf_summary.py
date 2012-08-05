# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Produce summary statistics from an annotated VCF file'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

import numpy as np

from   itertools          import izip

from   glu.lib.fileutils  import table_writer
from   glu.lib.seqlib.vcf import VCFReader


Ti_alleles = set( [('A','G'), ('G','A'), ('C','T'), ('T','C')] )


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

  ex_ti_tv = np.zeros( (2,2,n), dtype=int )
  ge_ti_tv = np.zeros( (2,2,n), dtype=int )

  ex_novel = np.zeros(n, dtype=int)
  ge_novel = np.zeros(n, dtype=int)

  ex_known = np.zeros(n, dtype=int)
  ge_known = np.zeros(n, dtype=int)

  # Bloody double negatives
  not_variant = set(['.','0','0/0','./0','0/.','./.','0|.','.|0'])

  bad = set(['SegDup','Repeat'])

  for v in vcf:
    #if bad.intersection(v.filter):
    #  continue

    #if 'SegDup' in v.filter or 'Repeat' in v.filter:
    #  continue

    #if 'SegDup' not in v.filter and 'Repeat' not in v.filter:
    #  continue

    #if float(v.qual)<30:
    #  continue

    # Count SNPs only for now
    if len(v.ref)!=1 or len(v.var)!=1 or len(v.var[0])!=1:
      continue

    #missing = np.fromiter( ('.'  in g             for g in v.genos), dtype=bool, count=n )
    #if missing.sum()>n/2:
    #  continue

    variant = np.fromiter( (g[0] not in not_variant for g in v.genos), dtype=bool, count=n )

    if not variant.sum():
      continue

    infomap = {}
    for inf in v.info:
      if '=' in inf:
        key,value = inf.split('=',1)
      else:
        key,value = inf,''

      infomap[key] = value

    if not v.names and int(infomap.get('REFVAR_OUTGROUP_COUNT',0))>80:
      continue

    if v.names:
      ge_known += variant
    else:
      ge_novel += variant

    alleles = v.ref,v.var[0]

    ge_ti_tv[bool(v.names),bool(alleles in Ti_alleles)] += variant

    if 'CDS' in infomap.get('GENE_LOCATION',[]):
      if v.names:
        ex_known += variant
      else:
        ex_novel += variant

      ex_ti_tv[bool(v.names)][bool(alleles in Ti_alleles)] += variant

  ge_ti = ge_ti_tv[0][1] + ge_ti_tv[1][1]
  ge_tv = ge_ti_tv[0][0] + ge_ti_tv[1][0]
  ex_ti = ex_ti_tv[0][1] + ex_ti_tv[1][1]
  ex_tv = ex_ti_tv[0][0] + ex_ti_tv[1][0]

  out.writerow(['SAMPLE','GENOME_KNOWN_VARIANTS','GENOME_NOVEL_VARIANTS','GENOME_NOVEL_FRACTION',
                         'GENOME_TiTv_ALL', 'GENOME_TiTv_KNOWN', 'GENOME_TiTv_NOVEL',
                         'EXOME_KNOWN_VARIANTS','EXOME_NOVEL_VARIANTS','EXOME_NOVEL_FRACTION',
                         'EXOME_TiTv_ALL', 'EXOME_TiTv_KNOWN', 'EXOME_TiTv_NOVEL'])

  rows = izip(samples,ge_known,ge_novel,ge_novel/(ge_novel+ge_known),
                      ge_ti/ge_tv,
                      ge_ti_tv[1][1]/ge_ti_tv[1][0],
                      ge_ti_tv[0][1]/ge_ti_tv[0][0],
                      ex_known,
                      ex_novel,
                      ex_novel/(ex_novel+ex_known),
                      ex_ti/ex_tv,
                      ex_ti_tv[1][1]/ex_ti_tv[1][0],
                      ex_ti_tv[0][1]/ex_ti_tv[0][0])

  out.writerows(rows)


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
