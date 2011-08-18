# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Plot intensity and allelic ratio data from indexed GDAT files'
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

from   glu.lib.fileutils    import table_reader, table_writer, cook_table, table_options

from   glu.modules.cnv.gdat import GDATIndex, get_gcmodel, gc_correct
from   glu.modules.cnv.plot import plot_chromosome


def get_assay_normal_mask(lrr,chrom_indices,events,options):
  normal_mask = np.ones_like(lrr,dtype=bool)

  if options.chromid and options.segstart and options.segstop:
    for chrom,chrom_events in groupby(events,key=attrgetter(options.chromid)):
      if chrom.startswith('chr'):
        chrom = chrom[3:]
      if chrom.upper()=='MT':
        chrom = 'M'

      pos,index    = chrom_indices[chrom]

      for i,event in enumerate(events):
        start      = int(getattr(event,options.segstart))
        stop       = int(getattr(event,options.segstop ))
        event_mask = (pos>=start)&(pos<stop)

        normal_mask[index] &= ~event_mask

  return normal_mask


def get_chrom_normal_mask(pos,events,options):
  normal_mask = np.ones_like(pos,dtype=bool)

  if options.chromid and options.segstart and options.segstop:
    for i,event in enumerate(events):
      start        = int(getattr(event,options.segstart))
      stop         = int(getattr(event,options.segstop ))
      event_mask   = (pos>=start)&(pos<=stop)
      normal_mask &= ~event_mask

  return normal_mask


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('events',       help='Tablular or delimited file of events')
  parser.add_argument('index',        help='Index file created by cnv.index_gdats')

  parser.add_argument('--gcmodel',    metavar='GCM', help='GC model file to use')
  parser.add_argument('--gcmodeldir', metavar='DIR', help='GC models directory')
  parser.add_argument('--outdir',     metavar='DIR', help='Plot output directory', default='.')
  parser.add_argument('--template',   metavar='VAL', default='{GDAT}_{ASSAY}_chr{CHROMOSOME}.png',
                                      help='Plot name template (default={GDAT}_{ASSAY}_chr{CHROMOSOME}.png)')
  parser.add_argument('--title',      metavar='VAL', default='Intensity plot of {GDAT}:{ASSAY} chr{CHROMOSOME}',
                                      help="Plot name template (default='Intensity plot of {GDAT}:{ASSAY} chr{CHROMOSOME}')")
  parser.add_argument('--assayid',    metavar='COLUMN',  default='ASSAY_ID',
                                      help='Assay ID column name (default=ASSAY_ID)')
  parser.add_argument('--chromid',    metavar='COLUMN',  default='CHROMOSOME',
                                      help='Chromosome column name (default=CHROMOSOME)')
  parser.add_argument('--segstart',   metavar='COLUMN', default='SEGSTART',
                                      help='Segment start index column name (default=SEGSTART)')
  parser.add_argument('--segstop',    metavar='COLUMN', default='SEGSTOP',
                                      help='Segment stop index column name (default=SEGSTOP)')
  table_options(parser)

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)
  indexfile = GDATIndex(options.index)

  gcmodels  = {}
  chip_indices = {}

  events    = table_reader(options.events)
  events    = cook_table(events,options)
  header    = next(events)
  Event     = namedtuple('Event', header)._make
  events    = [ Event(e) for e in events ]

  events.sort(key=attrgetter(options.assayid,options.chromid))

  for assay,assay_events in groupby(events,key=attrgetter(options.assayid)):
    assay_events = list(assay_events)

    locations  = indexfile.get(assay)

    if not locations:
      print 'ASSAY %s not found' % (assay)
      continue

    for gdat,offset in locations:
      gdatname      = '_'.join(os.path.basename(gdat.filename).replace('.gdat','').split('_')[:2])
      manifest      = gdat.attrs['ManifestName'].replace('.bpm','')

      chrom_indices = chip_indices.get(manifest)

      if chrom_indices is None:
        print 'Loading mapping information for %s...' % manifest
        chrom_indices = chip_indices[manifest] = gdat.chromosome_index

      print 'GDAT:',gdatname,assay

      assay_id,genos,lrr,baf = gdat.cnv_data(offset)
      normal_mask  = get_assay_normal_mask(lrr,chrom_indices,assay_events,options)
      mask         = normal_mask&(lrr>=-2)&(lrr<=2)
      lrr         -= lrr[mask].mean()

      if gccorrect:
        gcmodel  = gcmodels.get(manifest)
        if gcmodel is None:
          print 'Loading GC/CpG model for %s...' % manifest
          filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
          gcmodels[manifest] = gcmodel = get_gcmodel(filename,chrom_indices)

        gcdesign,gcmask = gcmodel
        lrr_adj         = gc_correct(lrr, gcdesign, gcmask&normal_mask,minval=-3,maxval=3)

      for chrom,chrom_events in groupby(assay_events,key=attrgetter(options.chromid)):
        if chrom.startswith('chr'):
          chrom = chrom[3:]

        chrom_events = list(chrom_events)
        pos,index    = chrom_indices[chrom]

        chrom_genos  = genos[index]
        chrom_baf    = baf[index]
        chrom_lrr    = None

        normal       = get_chrom_normal_mask(pos, chrom_events, options)

        if gccorrect:
          chrom_lrr  = lrr_adj[index]
          valid      = (chrom_lrr>=-2)&(chrom_lrr<=2)

          if valid.sum()<100:
            chrom_lrr = None

        if chrom_lrr is None:
          chrom_lrr  = lrr[index]
          valid      = (chrom_lrr>=-2)&(chrom_lrr<=2)

        normal       =   normal &valid
        abnormal     = (~normal)&valid

        print '  chr%-3s Normal: LRR=%6.3f += %.3f, BAF=%.2f +- %.2f' %  \
                    (chrom,chrom_lrr[normal].mean(),1.96*chrom_lrr[normal].std(),
                           chrom_baf[normal].mean(),1.96*chrom_baf[normal].std())
        print '       Abnormal: LRR=%6.3f += %.3f, BAF=%.2f +- %.2f' %  \
                          (chrom_lrr[abnormal].mean(),1.96*chrom_lrr[abnormal].std(),
                           chrom_baf[abnormal].mean(),1.96*chrom_baf[abnormal].std())

        variables = chrom_events[0]._asdict()
        variables.update(GDAT=gdatname,ASSAY=assay,CHROMOSOME=chrom)

        plotname = options.template.format(**variables)
        plotname = '%s/%s' % (options.outdir,plotname)
        title    = options.title.format(**variables)

        plot_chromosome(plotname, pos[valid], chrom_lrr[valid], chrom_baf[valid], genos=chrom_genos[valid],
                        title=title, events=chrom_events, startattr=options.segstart, stopattr=options.segstop)


if __name__=='__main__':
  main()
