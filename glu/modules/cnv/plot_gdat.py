 # -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Plot intensity and allelic ratio data from a single GDAT file'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

import numpy as np
import scipy
import scipy.stats

from   itertools            import groupby
from   operator             import itemgetter,attrgetter

from   collections          import namedtuple

from   glu.lib.fileutils    import table_reader, table_writer

from   glu.modules.cnv.gdat import GDATFile, get_gcmodel, gc_correct
from   glu.modules.cnv.plot import plot_chromosome


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdat',         help='Input GDAT file')

  parser.add_argument('--signal',     metavar='TYPE', default='norm',
                                      help='LRR/BAF signal processing: raw or norm (default)')
  parser.add_argument('--chromosomes',metavar='CHRS', default='',
                                      help='Comma separated list of chromosomes to plot (blank for all)')
  parser.add_argument('--gcmodel',    metavar='GCM', help='GC model file to use')
  parser.add_argument('--gcmodeldir', metavar='DIR', help='Directory containing GC model files for automatic selection based on manifest name')
  parser.add_argument('--outdir',     metavar='DIR', help='Plot output directory', default='.')
  parser.add_argument('--template',   metavar='VAL', default='{GDAT}_{ASSAY}_chr{CHROMOSOME}.png',
                                      help='Plot name template (default={GDAT}_{ASSAY}_chr{CHROMOSOME}.png)')
  parser.add_argument('--title',      metavar='VAL', default='Intensity plot of {GDAT}:{ASSAY} chr{CHROMOSOME}',
                                      help="Plot name template (default='Intensity plot of {GDAT}:{ASSAY} chr{CHROMOSOME}')")

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)
  gdat      = GDATFile(options.gdat)

  gdatname  = '_'.join(os.path.basename(gdat.filename).replace('.gdat','').split('_')[:2])
  manifest  = gdat.attrs['ManifestName'].replace('.bpm','')

  print 'Loading manifest for %s...' % manifest
  chrom_index = gdat.chromosome_index

  if gccorrect:
    print 'Loading GC/CpG model for %s...' % manifest
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    gcdesign,gcmask = get_gcmodel(filename,chrom_index)

  chromosomes = [ c.strip().upper() for c in options.chromosomes.split(',')  if c.strip() ]

  if not chromosomes:
    chromosomes = sorted(c for c in chrom_index if c)

  for offset,assay in enumerate(gdat.samples):
    print 'GDAT:',gdatname,assay
    assay_id,genos,lrr,baf  = gdat.cnv_data(offset, raw=options.signal.lower()=='raw')
    mask     = np.isfinite(lrr)&(lrr>=-2)&(lrr<=2)
    lrr     -= lrr[mask].mean()

    if gccorrect:
      lrr_adj= gc_correct(lrr, gcdesign, gcmask)

    for chrom in chromosomes:
      chrom = chrom.upper()
      if chrom.startswith('CHR'):
        chrom = chrom[3:]
      if chrom.upper()=='MT':
        chrom = 'M'

      pos,index    = chrom_index[chrom]
      chrom_baf    = baf[index]
      chrom_lrr    = None

      if gccorrect:
        chrom_lrr = lrr_adj[index]
        mask = np.isfinite(chrom_lrr)&np.isfinite(chrom_baf)&(chrom_lrr>=-2)&(chrom_lrr<=2)

        if mask.sum()<100:
          chrom_lrr = None

      if chrom_lrr is None:
        chrom_lrr = lrr[index]
        mask      = np.isfinite(chrom_lrr)&np.isfinite(chrom_baf)&(chrom_lrr>=-2)&(chrom_lrr<=2)

      print '  chr%-3s LRR=%6.3f +/ %.3f, BAF=%.2f +/ %.2f' %  \
                  (chrom,chrom_lrr[mask].mean(),1.96*chrom_lrr[mask].std(),
                         chrom_baf[mask].mean(),1.96*chrom_baf[mask].std())

      variables = dict(GDAT=gdatname,ASSAY=assay,CHROMOSOME=chrom)

      plotname = options.template.format(**variables)
      plotname = '%s/%s' % (options.outdir,plotname)
      title    = options.title.format(**variables)

      plot_chromosome(plotname, pos[mask], chrom_lrr[mask], chrom_baf[mask],
                      title=title)


if __name__=='__main__':
  main()
