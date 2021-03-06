# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Plot intensity and allelic ratio data from a single GDAT file'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

import numpy as np

from   glu.lib.fileutils         import table_reader

from   glu.modules.cnv.gdat      import GDATFile, get_gcmodel, gc_correct
from   glu.modules.cnv.plot      import plot_chromosome

from   glu.lib.genolib.transform import GenoTransform


def load_extravars(filename):
  extravars = {}

  if not filename:
    return extravars,extravars

  data   = table_reader(filename)
  header = next(data)

  if 'GDAT' not in header:
    raise ValueError('extravars must have a "GDAT" column')

  if 'ASSAY' not in header:
    raise ValueError('extravars must have an "ASSAY" column')

  for row in data:
    vars  = dict(zip(header,row))
    gdat  = vars.get('GDAT')
    assay = vars.get('ASSAY')

    if not gdat or not assay:
      continue

    key = gdat,assay

    if key in extravars:
      raise KeyError('duplicate extra variable information for %s' % key)

    extravars[key] = vars

  extradummy = dict( (k,'') for k in header )

  return extravars,extradummy


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdat',         help='Input GDAT file')

  parser.add_argument('--signal',     metavar='TYPE', default='norm', choices=['norm','raw'],
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
  parser.add_argument('--extravars',  metavar='FILE', 
                                      help='Table of extra variables on each sample to add information to the '
                                      'title and filename template. Table must have GDAT and ASSAY columns.')
  parser.add_argument('--includesamples', metavar='FILE', action='append',
                                      help='List of samples to include, all others will be skipped')
  parser.add_argument('--excludesamples', metavar='FILE', action='append',
                                      help='List of samples to exclude, only samples not present will be kept')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)
  gdat      = GDATFile(options.gdat)

  gdatraw   = os.path.basename(gdat.filename).replace('.gdat','')
  gdatname  = '_'.join(gdatraw.split('_')[:2])
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

  transform = GenoTransform.from_object(options)

  if transform is not None:
    include  = transform.samples.include
    exclude  = transform.samples.exclude
  else:
    include  = exclude = None

  extravars,extradummy = load_extravars(options.extravars)

  for offset,assay in enumerate(gdat.samples):
    if include is not None and assay not in include:
      continue

    if exclude is not None and assay in exclude:
      continue

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
      chrom_genos  = genos[index]

      if gccorrect:
        chrom_lrr = lrr_adj[index]
        mask = np.isfinite(chrom_lrr)&np.isfinite(chrom_baf)&(chrom_lrr>=-2)&(chrom_lrr<=2)

        if mask.sum()<100:
          chrom_lrr = None

      if chrom_lrr is None:
        chrom_lrr = lrr[index]
        mask      = np.isfinite(chrom_lrr)&np.isfinite(chrom_baf)&(chrom_lrr>=-2)&(chrom_lrr<=2)

      hets    = (chrom_genos=='AB')&mask
      het_baf = chrom_baf[hets]

      print '  chr%-3s LRR=%6.3f +- %.3f, BAF=%.2f +- %.2f BAF_AB=%.2f +- %.2f' %  \
                  (chrom,chrom_lrr[mask].mean(),1.96*chrom_lrr[mask].std(),
                         chrom_baf[mask].mean(),1.96*chrom_baf[mask].std(),
                         het_baf.mean(),1.96*het_baf.std())

      extra     = extravars.get( (gdatraw,assay), {} )
      if not extra:
        extra   = extravars.get( (gdatname,assay), {} )
      defvars   = dict(GDAT=gdatname,ASSAY=assay,CHROMOSOME=chrom)

      variables = extradummy.copy()
      variables.update(extra)
      variables.update(defvars)

      plotname = options.template.format(**variables)
      plotname = '%s/%s' % (options.outdir,plotname)
      title    = options.title.format(**variables)

      plot_chromosome(plotname, pos[mask], chrom_lrr[mask], chrom_baf[mask],
                      title=title)


if __name__=='__main__':
  main()
