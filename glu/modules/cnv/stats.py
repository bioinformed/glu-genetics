# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Generate statistics on LRR and BAF for original and renormalized values'
__copyright__ = 'Copyright (c) 2012, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

from   math                 import modf
from   itertools            import repeat

import numpy       as np
import scipy.stats as stats

from   glu.lib.progressbar  import progress_loop
from   glu.lib.fileutils    import table_writer
from   glu.modules.cnv.gdat import GDATFile, BatchTableWriter, parallel_gdat_iter, \
                                   create_gdat_qn, baf_deviation


def skewtest(s):
  if len(s)<10:
    return 1.0
  return stats.skewtest(s)[0]


def normalize_chromosome(chrom):
  chrom = chrom.strip().upper()
  if not chrom:
    return None
  if chrom.startswith('CHR'):
    chrom = chrom[3:]
  if chrom=='MT':
    chrom = 'M'
  return chrom


def parse_chromosome_region(region):
  if ':' not in region:
    chrom      = normalize_chromosome(region)
    start      = 0
    stop       = sys.maxint
  else:
    chrom,rest = region.split(':',1)
    start,stop = rest.split('-')
    chrom      = normalize_chromosome(chrom)
    start      = int(start.replace(',',''))-1
    stop       = int(stop.replace(',',''))

  return chrom,start,stop


def parse_chromosome_regions(regions):
  for region in regions:
    region = region.strip()
    if region:
      yield parse_chromosome_region(region)


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdat',                         help='Input GDAT file')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                      help='Output results (default is "-" for standard out)')
  parser.add_argument('--chrmask',    action='append',
                                      help='Show only chrom:start-stop region for events on chrom.  May be specified multiple times.')
  parser.add_argument('--chromosomes',metavar='CHRS', default='',
                                      help='Comma separated list of chromosomes to plot (blank for all)')
  parser.add_argument('-P', '--progress', action='store_true',
                         help='Show analysis progress bar, if possible')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gdat      = GDATFile(options.gdat,'r+')
  name      = os.path.splitext(os.path.basename(options.gdat))[0]
  n         = gdat.sample_count
  s         = gdat.snp_count
  samples   = gdat.samples
  BAF_orig  = gdat['BAF']
  LRR_orig  = gdat['LRR']
  genotypes = gdat['Genotype']

  if 'LRR_QN' in gdat:
    BAF_norm  = gdat['BAF_QN']
    LRR_norm  = gdat['LRR_QN']
  else:
    BAF_empty    = np.empty( (BAF_orig.shape[1],), dtype=float )
    LRR_empty    = np.empty( (LRR_orig.shape[1],), dtype=float )
    BAF_empty[:] = np.nan
    LRR_empty[:] = np.nan
    BAF_norm     = repeat(BAF_empty)
    LRR_norm     = repeat(LRR_empty)

  print 'samples=%d, SNPs=%d' % (n,s)

  np.set_printoptions(linewidth=120,precision=2,suppress=True)
  np.seterr(all='ignore')

  valid_snps,mask = s,None

  if options.chromosomes or options.chrmask:
    chrom_index = gdat.chromosome_index

    mask = np.zeros(gdat.snp_count, dtype=bool)

    regions = (options.chrmask or []) + (options.chromosomes or '').split(',')
    regions = parse_chromosome_regions(regions)

    for chrom,start,stop in regions:
      pos,indices          = chrom_index[chrom]
      region_mask          = (pos>=start)&(pos<stop)
      region_indices       = indices[region_mask]
      mask[region_indices] = True

    valid_snps = mask.sum()
    print 'Found %d valid probes' % valid_snps

  data = parallel_gdat_iter(LRR_orig,LRR_norm,BAF_orig,BAF_norm,genotypes)

  if options.progress:
    data = progress_loop(data, length=n, units='samples', label='ASSAY: ')

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['ASSAY_ID','GDAT',
                'COUNT_ORIG','LRR_ORIG_MEAN','LRR_ORIG_SD','BAF_ORIG_SKEW',
                             'BDEV_ORIG_MEAN',  'BDEV_ORIG_SD',
                             'BAF_ORIG_AA_MEAN','BAF_ORIG_AA_SD',
                             'BAF_ORIG_AB_MEAN','BAF_ORIG_AB_SD',
                             'BAF_ORIG_BB_MEAN','BAF_ORIG_BB_SD',
                             'MISSING_RATE_ORIG',
                'COUNT_NORM','LRR_NORM_MEAN','LRR_NORM_SD','BAF_NORM_SKEW',
                             'BDEV_NORM_MEAN',  'BDEV_NORM_SD',
                             'BAF_NORM_AA_MEAN','BAF_NORM_AA_SD',
                             'BAF_NORM_AB_MEAN','BAF_NORM_AB_SD',
                             'BAF_NORM_BB_MEAN','BAF_NORM_BB_SD',
                             'MISSING_RATE_NORM'])

  for i,(lrr_orig,lrr_norm,baf_orig,baf_norm,genos) in enumerate(data):
    if not options.progress and options.output!='-':
      print '  Sample %5d / %d' % (i+1,n)

    if mask is not None:
      lrr_orig     = lrr_orig[mask]
      lrr_norm     = lrr_norm[mask]
      baf_orig     = baf_orig[mask]
      baf_norm     = baf_norm[mask]
      genos        = genos[mask]

    assay          = samples[i]
    AA             = genos=='AA'
    AB             = genos=='AB'
    BB             = genos=='BB'
    NN             = genos=='  '

    bdev_orig      = baf_deviation(baf_orig,genos)
    bdev_norm      = baf_deviation(baf_norm,genos)

    valid          = np.isfinite(baf_orig)&np.isfinite(lrr_orig)
    count_orig     = valid.sum()
    baf_orig_skew  = skewtest(baf_orig[valid&AB])
    lrr_orig_mu    = lrr_orig[valid   ].mean()
    lrr_orig_sd    = lrr_orig[valid   ].std(ddof=1)
    bdev_orig_mu   = bdev_orig[valid  ].mean()
    bdev_orig_sd   = bdev_orig[valid  ].std(ddof=1)
    baf_orig_aa_mu = baf_orig[valid&AA].mean()
    baf_orig_aa_sd = baf_orig[valid&AA].std(ddof=1)
    baf_orig_ab_mu = baf_orig[valid&AB].mean()
    baf_orig_ab_sd = baf_orig[valid&AB].std(ddof=1)
    baf_orig_bb_mu = baf_orig[valid&BB].mean()
    baf_orig_bb_sd = baf_orig[valid&BB].std(ddof=1)
    missing_orig   = 1-valid.sum()/valid_snps
    #missing_orig  = (valid&NN).sum()/valid_snps

    valid          = np.isfinite(baf_norm)&np.isfinite(lrr_norm)
    count_norm     = valid.sum()
    baf_norm_skew  = skewtest(baf_norm[valid&AB])
    lrr_norm_mu    = lrr_norm[valid   ].mean()
    lrr_norm_sd    = lrr_norm[valid   ].std(ddof=1)
    bdev_norm_mu   = bdev_norm[valid  ].mean()
    bdev_norm_sd   = bdev_norm[valid  ].std(ddof=1)
    baf_norm_aa_mu = baf_norm[valid&AA].mean()
    baf_norm_aa_sd = baf_norm[valid&AA].std(ddof=1)
    baf_norm_ab_mu = baf_norm[valid&AB].mean()
    baf_norm_ab_sd = baf_norm[valid&AB].std(ddof=1)
    baf_norm_bb_mu = baf_norm[valid&BB].mean()
    baf_norm_bb_sd = baf_norm[valid&BB].std(ddof=1)
    missing_norm   = 1-valid.sum()/valid_snps
    #missing_norm  = (valid&NN).sum()/valid_snps

    out.writerow([assay,name,
                  count_orig,lrr_orig_mu,lrr_orig_sd,baf_orig_skew,
                             bdev_orig_mu,bdev_orig_sd,
                             baf_orig_aa_mu,baf_orig_aa_sd,
                             baf_orig_ab_mu,baf_orig_ab_sd,
                             baf_orig_bb_mu,baf_orig_bb_sd,missing_orig,
                  count_norm,lrr_norm_mu,lrr_norm_sd,baf_norm_skew,
                             bdev_norm_mu,bdev_norm_sd,
                             baf_norm_aa_mu,baf_norm_aa_sd,
                             baf_norm_ab_mu,baf_norm_ab_sd,
                             baf_norm_bb_mu,baf_norm_bb_sd,missing_norm,
                 ])

  gdat.close()


if __name__ == '__main__':
  main()
