# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Plot intensity and allelic ratio data from indexed GDAT files'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os

import numpy as np

from   itertools            import groupby
from   operator             import attrgetter

from   collections          import namedtuple

from   glu.lib.fileutils    import table_reader, cook_table, table_options
from   glu.lib.recordtype   import recordtype

from   glu.modules.cnv.gdat import GDATIndex, get_gcmodel, gc_correct
from   glu.modules.cnv.plot import plot_chromosome


Event = recordtype('Event', 'chrom start stop mask')


def fit_gmm1(components,x):
  from scikits.learn.mixture   import GMM, DPGMM, VBGMM

  x = np.asanyarray(x).reshape(-1,1)
  #gmm_baf = VBGMM(components,verbose=True,min_covar=0.01)
  gmm_baf = GMM(components)
  gmm_baf.fit(x)

  mu    = np.array(gmm_baf.means,  dtype=float).reshape(-1)
  order = mu.argsort()
  mu    = mu[order]
  sd    = np.array(gmm_baf.covars,  dtype=float).reshape(-1)[order]**0.5
  ws    = np.array(gmm_baf.weights, dtype=float).reshape(-1)[order]
  logL  = gmm_baf.score(x)

  return logL,mu,sd,ws


def fit_gmm2(k,x):
  from   glu.modules.cnv.gmm     import gaussian_mixture

  w  = np.ones(k)/k
  mu = np.array([ i/(k-1) for i in range(k) ], dtype=float)
  sd = np.array([0.01]+[0.05]*(k-2)+[0.01], dtype=float)

  logL,w,mu,sd = gaussian_mixture(k,x, w=w, mu=mu, sd=sd, min_w=0.05, min_sd=0.005)

  return logL,mu,sd,w


def fit_gmm3(k,x,wnoise=0.02):
  import mixture

  k -= 2

  dataset = mixture.DataSet()
  dataset.fromArray(x)

  hom_a   = mixture.NormalDistribution(0,0.05)
  hom_b   = mixture.NormalDistribution(1,0.05)

  hom_a.min_sigma = 0.005
  hom_b.min_sigma = 0.005

  noise   = mixture.UniformDistribution(0,1)

  hets    = [ mixture.NormalDistribution( (i+1)/(k+1), 0.05) for i in range(k) ]
  for het in hets:
    het.min_sigma = 0.03

  dists   = [hom_a]+hets+[hom_b]
  homw    = 0.25 if k else (0.50-wnoise/2)
  weights = [homw]+[(0.5-wnoise)/k if k else 0]*k+[homw]

  if wnoise>0:
    dists.append(noise)
    weights.append(wnoise)

  mix = mixture.MixtureModel(len(dists), weights, dists)

  post,logL = mix.EM(dataset, 100, delta=0.001, silent=True)

  w     = np.array(mix.pi,                                                       dtype=float)
  mu    = np.array([hom_a.mu   ] + [het.mu    for het in hets] + [hom_b.mu    ], dtype=float)
  sd    = np.array([hom_a.sigma] + [het.sigma for het in hets] + [hom_b.sigma ], dtype=float)

  order = mu.argsort()
  mu    = mu[order]
  sd    = sd[order]
  w     =  w[order]

  if 0:
    print '  %d band component mixture: logL=%f' % (k+2,logL)
    print '      AIC: %.2f' % (4*(k+2)-2*logL)
    print '    MEANS:',mu
    print '       SD:',sd
    print '  WEIGHTS:',w
    print

  return logL,mu,sd,w


def fit_gmm_baf(baf,max_baf=4):
  #mask = (baf>0)&(baf<1)

  mask = np.isfinite(baf)

  if mask.sum()<100:
    return None

  baf  = baf[mask]

  import time

  t0        = time.time()
  fits      = [ fit_gmm3(c,baf) for c in range(2,max_baf+1) ]
  AIC       = np.array([ (20*c-2*f[0]) for c,f in enumerate(fits,2) ], dtype=float)
  best      = AIC.argmin()

  logL,mu,sd,ws = fits[best]
  c             = best+2

  if 0:
    print '        TIME: %.2f s' % (time.time()-t0)
    print '  COMPONENTS:',c
    print '        logL:',logL
    print '         AIC:',2*(c-logL)
    print '       MEANS:',mu
    print '          SD:',sd
    print '     WEIGHTS:',ws
    print

  return c,mu,sd,ws


def norm_chromosome(chrom):
  if chrom.startswith('chr'):
    chrom = chrom[3:]
  if chrom.upper()=='MT':
    chrom = 'M'
  return chrom


def get_assay_normal_mask(lrr,chrom_indices,events,options):
  normal_mask = np.ones_like(lrr,dtype=bool)

  if options.chromid and options.segstart and options.segstop:
    for chrom,chrom_events in groupby(events,key=attrgetter(options.chromid)):
      chrom        = norm_chromosome(chrom)
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


def build_chromosome_masks(pos, baf, lrr, mask=None, events=None, chromid=None, segstart=None, segstop=None):
  if mask is not None:
    valid_mask  = np.isfinite(lrr)&mask
  else:
    valid_mask  = mask

  normal_mask   =  valid_mask.copy()
  abnormal_mask = np.zeros_like(normal_mask,dtype=bool)

  eventrecs     = []

  for i,event in enumerate(events or []):
    chrom          = getattr(event,chromid)
    chrom          = norm_chromosome(chrom)
    start          = int(getattr(event,segstart))
    stop           = int(getattr(event,segstop ))

    event_mask     = pos >= start-1
    event_mask    &= pos <= stop
    event_mask    &= valid_mask

    normal_mask   &= ~event_mask
    abnormal_mask |=  event_mask

    eventrecs.append( Event(chrom,start,stop,event_mask) )

  return normal_mask,abnormal_mask,eventrecs


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
  parser.add_argument('--segstart',   metavar='COLUMN', default='SEG_START',
                                      help='Segment start index column name (default=SEG_START)')
  parser.add_argument('--segstop',    metavar='COLUMN', default='SEG_STOP',
                                      help='Segment stop index column name (default=SEG_STOP)')
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
      print 'ASSAY %s not found' % assay
      continue

    for gdat,offset in locations:
      gdatname      = '_'.join(os.path.splitext(os.path.basename(gdat.filename))[0].split('_')[:2])
      manifest      = gdat.attrs['ManifestName'].replace('.bpm','')

      chrom_indices = chip_indices.get(manifest)

      if chrom_indices is None:
        print 'Loading mapping information for %s...' % manifest
        chrom_indices = chip_indices[manifest] = gdat.chromosome_index

      print 'GDAT:',gdatname,assay

      assay_id,genos,lrr,baf = gdat.cnv_data(offset)
      normal_mask  = get_assay_normal_mask(lrr,chrom_indices,assay_events,options)
      mask         = normal_mask&np.isfinite(lrr)&(lrr>=-2)&(lrr<=2)

      if not mask.sum():
        print '  ASSAY %s DOES NOT CONTAIN VALID LRR/BAF data' % assay
        continue

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
        chrom        = norm_chromosome(chrom)
        chrom_events = list(chrom_events)
        pos,index    = chrom_indices[chrom]

        chrom_genos  = genos[index]
        chrom_baf    = baf[index]
        chrom_lrr    = None

        if gccorrect:
          chrom_lrr  = lrr_adj[index]
          valid      = (chrom_lrr>=-2)&(chrom_lrr<=2)

          if valid.sum()<100:
            chrom_lrr = None

        if chrom_lrr is None:
          chrom_lrr  = lrr[index]
          valid      = (chrom_lrr>=-2)&(chrom_lrr<=2)

        (normal_mask,
         abnormal_mask,
         eventrecs)  = build_chromosome_masks(pos, chrom_baf, chrom_lrr, mask=valid, events=chrom_events,
                                              chromid=options.chromid, segstart=options.segstart,
                                              segstop=options.segstop)

        lrr_normal   = chrom_lrr[normal_mask]
        baf_normal   = chrom_baf[normal_mask]

        lrr_abnormal = chrom_lrr[abnormal_mask]
        baf_abnormal = chrom_baf[abnormal_mask]

        print '  chr%-3s Normal: probes=%5d LRR=%6.3f += %.3f, BAF=%.2f +- %.2f' %  \
                    (chrom,len(lrr_normal),
                           lrr_normal.mean(),1.96*lrr_normal.std(),
                           baf_normal.mean(),1.96*baf_normal.std())

        print '       Abnormal: probes=%5d LRR=%6.3f += %.3f, BAF=%.2f +- %.2f' %  \
                          (len(lrr_abnormal),
                           lrr_abnormal.mean(),1.96*lrr_abnormal.std(),
                           baf_abnormal.mean(),1.96*baf_abnormal.std())

        baf_gmm_normal   = fit_gmm_baf(baf_normal,  max_baf=3)
        baf_gmm_abnormal = fit_gmm_baf(baf_abnormal,max_baf=4)

        variables = chrom_events[0]._asdict()
        variables.update(GDAT=gdatname,ASSAY=assay,CHROMOSOME=chrom)

        plotname = options.template.format(**variables)
        plotname = '%s/%s' % (options.outdir,plotname)
        title    = options.title.format(**variables)

        plot_chromosome(plotname, pos, chrom_lrr, chrom_baf, normal_mask, abnormal_mask,
                        baf_gmm_normal, baf_gmm_abnormal,
                        genos=chrom_genos, title=title, events=eventrecs)


if __name__=='__main__':
  main()
