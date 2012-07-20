# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Plot intensity and allelic ratio data from indexed GDAT files'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os

import numpy as np
import scipy.stats

from   math                      import sqrt, fabs
from   itertools                 import groupby
from   operator                  import attrgetter

from   collections               import namedtuple

from   glu.lib.fileutils         import table_reader, table_writer, cook_table, table_options
from   glu.lib.recordtype        import recordtype

from   glu.modules.cnv.gdat      import GDATIndex, get_gcmodel, gc_correct
from   glu.modules.cnv.plot      import plot_chromosome
from   glu.modules.cnv.normalize import quantile

from   glu.lib.genolib.transform import GenoTransform


Event = recordtype('Event', 'chrom start stop fields lrr baf baf_model state pmosaic')


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
    print '  %d band component mixture: logL=%f' % (k,logL)
    print '      AIC: %.2f' % (4*(k+2)-2*logL)
    print '    MEANS:',mu
    print '       SD:',sd
    print '  WEIGHTS:',w
    print

  return logL,mu,sd,w


def fit_gmm4(k,x,wnoise=0.02):
  import mixture

  dataset = mixture.DataSet()
  dataset.fromArray(x)

  if k==0:
    wnoise = 1

  noise   = mixture.UniformDistribution(0,1)
  hets    = [ mixture.NormalDistribution( (i+1)/(k+1), 0.05) for i in range(k) ]
  for het in hets:
    het.min_sigma = 0.03

  dists   = hets[:]
  weights = [(1-wnoise)/k if k else 0]*k

  if wnoise>0:
    dists.append(noise)
    weights.append(wnoise)

  mix = mixture.MixtureModel(len(dists), weights, dists)

  try:
    post,logL = mix.EM(dataset, 100, delta=0.001, silent=True)
  except mixture.ConvergenceFailureEM:
    return None

  w     = np.array(mix.pi,                      dtype=float)
  mu    = np.array([het.mu    for het in hets], dtype=float)
  sd    = np.array([het.sigma for het in hets], dtype=float)

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


def fit_gmm_baf(baf,max_baf=2):
  #mask = (baf>0)&(baf<1)

  mask = np.isfinite(baf)

  if mask.sum()<100:
    return None

  baf  = baf[mask]

  import time

  t0        = time.time()
  fits      = [ fit_gmm4(c,baf) for c in range(0,max_baf+1) ]

  if None in fits:
    return None

  AIC       = np.array([ (20*c-2*f[0]) for c,f in enumerate(fits) ], dtype=float)
  best      = AIC.argmin()

  logL,mu,sd,ws = fits[best]
  c             = best #+2

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


def norm_position(p):
  if p is None:
    return None
  p = p.replace(',','')
  try:
    return int(p)
  except ValueError:
    return None


def build_normal_mask(s,chrom_indices,events,options):
  normal_mask = np.ones(s,dtype=bool)

  if options.chromid and options.segstart and options.segstop:
    for chrom,chrom_events in groupby(events,key=attrgetter(options.chromid)):
      chrom        = norm_chromosome(chrom)
      pos,index    = chrom_indices[chrom]

      for i,event in enumerate(events):
        start      = norm_position(getattr(event,options.segstart,None))
        stop       = norm_position(getattr(event,options.segstop, None))

        if start is None or stop is None:
          normal_mask[index] = False
        else:
          event_mask = (pos>=start)&(pos<stop)
          normal_mask[index] &= ~event_mask

  elif options.chromid:
    for chrom,chrom_events in groupby(events,key=attrgetter(options.chromid)):
      chrom              = norm_chromosome(chrom)
      pos,index          = chrom_indices[chrom]
      normal_mask[index] = False

  return normal_mask


def build_events(lrr, baf, chrom_indices, valid_mask, events, options):
  eventrecs       = []

  if options.chromid :
    for i,event in enumerate(events or []):
      chrom         = getattr(event,options.chromid)
      chrom         = norm_chromosome(chrom)
      pos,indices   = chrom_indices[chrom]
      event_mask    = None
      start,stop    = None,None

      if options.segstart and options.segstop:
        start       = norm_position(getattr(event,options.segstart,None))
        stop        = norm_position(getattr(event,options.segstop, None))

        if start is not None and stop is not None:
          event_mask  = valid_mask[indices]
          event_mask &= pos >= start-1
          event_mask &= pos <= stop

      if event_mask is not None and event_mask.sum()>10:
        event_lrr   = lrr[indices][event_mask]
        event_baf   = baf[indices][event_mask]
      else:
        event_lrr   = None
        event_baf   = None

      eventrecs.append( Event(chrom,start,stop,event,event_lrr,event_baf,None,None,None) )

  return eventrecs


def valid_autosome_mask(s,chrom_indices):
  missing      = [np.array([], dtype=int)]*2
  chr6_pos,\
  chr6_indices = chrom_indices['6']
  hla          = chr6_indices[(chr6_pos>26000000)&(chr6_pos<34000000)]
  chrX         = chrom_indices.get('X', missing)[1]
  chrY         = chrom_indices.get('Y', missing)[1]
  chrXY        = chrom_indices.get('XY',missing)[1]
  chrM         = chrom_indices.get('M', missing)[1]

  valid        = np.ones(s, dtype=bool)
  valid[hla ]  = False
  valid[chrX]  = False
  valid[chrY]  = False
  valid[chrXY] = False
  valid[chrM]  = False

  return valid


def mirror_baf_distribution(baf,genos,mask,trim=0.20,threshold=0.0001,maxsize=100000):
  aa_mask       = genos=='AA'
  aa_mask      &= mask
  aa_mask      &= baf>threshold

  bb_mask       = genos=='BB'
  bb_mask      &= mask
  bb_mask      &= baf<(1-threshold)

  aa_baf        = baf[aa_mask]
  bb_baf        = baf[bb_mask]

  aa_baf.sort()
  bb_baf.sort()

  aa_limit      = int(len(aa_baf)*(1-trim))
  bb_limit      = int(len(bb_baf)*   trim )

  print 'AA: Extracting %d elements, threshold value=%f' % (aa_limit,aa_baf[aa_limit])
  print 'BB: Extracting %d elements, threshold value=%f' % (bb_limit,bb_baf[bb_limit])

  aa_baf        =  aa_baf[:aa_limit]
  bb_baf        =  bb_baf[bb_limit:]

  aa_baf_max    = aa_baf.max()
  bb_baf_min    = bb_baf.min()

  aa_baf        =   -aa_baf
  bb_baf        =  2-bb_baf

  np.random.shuffle(aa_baf)
  np.random.shuffle(bb_baf)

  if maxsize:
    aa_baf      = aa_baf[:maxsize]
    bb_baf      = bb_baf[:maxsize]

  return aa_baf,aa_baf_max,bb_baf,bb_baf_min


def augment_event_baf(baf,aa_baf,aa_baf_max,bb_baf,bb_baf_min):
  baf_mask = (baf>0)&(baf<1)
  aa_count = ((baf<aa_baf_max)&baf_mask).sum()
  bb_count = ((baf>bb_baf_min)&baf_mask).sum()
  baf_aug = np.concatenate( (baf[baf_mask],aa_baf[:aa_count],bb_baf[:bb_count]) )
  print '  Adding BAF values: AA<%f = %d, BB>%f = %d.  Size %d -> %d' % (aa_baf_max,aa_count,bb_baf_min,bb_count,len(baf),len(baf_aug))

  return baf_aug


def strip_event_baf(baf,aa_baf_max,bb_baf_min):
  baf_mask = (baf>aa_baf_max)&(baf<bb_baf_min)
  baf_aug  = baf[baf_mask]
  return baf_aug


def point_distance(a, b, c):
  t  = b[0]-a[0], b[1]-a[1]          # Vector ab
  dd = sqrt(t[0]**2+t[1]**2)         # Length of ab
  t  = t[0]/dd, t[1]/dd              # unit vector of ab
  n  = -t[1], t[0]                   # normal unit vector to ab
  ac = c[0]-a[0], c[1]-a[1]          # vector ac
  return fabs(ac[0]*n[0]+ac[1]*n[1]) # Projection of ac to n (the minimum distance)


def infer_mosaic_state(event,sd_threshold=0.25):
  if event.lrr is None:
    return

  lrr_mean  = event.lrr.mean()
  lrr_sd    = event.lrr.std(ddof=1)

  if event.baf_model and event.baf_model[0]==2:
    c,mu,sd,ws = event.baf_model
    d   = mu[1]-mu[0]

    p_neutral = d
    p_gain    = 2*d/(1-d)
    p_loss    = 2*d/(1+d)

    zero         = (0.0, 0.00)
    pure_neutral = (1.0, 0.00)
    pure_gain    = (1.0, 0.2682)
    pure_loss    = (1.0,-0.3427)

    d_neutral = point_distance(zero, pure_neutral, (p_neutral,lrr_mean))
    d_gain    = point_distance(zero, pure_gain,    (p_gain,   lrr_mean))
    d_loss    = point_distance(zero, pure_loss,    (p_loss,   lrr_mean))

    d_min     = min(d_neutral,d_gain,d_loss)

    if d_min==d_neutral:
      state = 'NEUTRAL'
      p     = p_neutral
    elif d_min==d_loss:
      state = 'LOSS'
      p     = p_loss
    else:
      state = 'GAIN'
      p     = p_gain

    event.state   = state
    event.pmosaic = max(0,min(1.,p))
  else:
    snr = abs(lrr_mean)/lrr_sd
    if snr>sd_threshold:
      if event.lrr.mean()<0:
        event.state = 'LOSS'
      else:
        event.state = 'GAIN'
    else:
      event.state   = 'NEUTRAL'


def update_event_fields(event):
  fields          = event.fields
  fields.Probes   = len(event.lrr)
  fields.STATE    = event.state
  fields.PMosaic  = '%.5f' % event.pmosaic if event.pmosaic is not None else ''
  fields.LRR_mean = '%.3f' % event.lrr.mean()
  fields.LRR_sd   = '%.3f' % event.lrr.std(ddof=1)

  if not event.baf_model:
    fields.BAF_BANDS = ''
    fields.BAF_means = ''
    fields.BAF_sds   = ''
  else:
    c,mu,sd,ws = event.baf_model

    if not c:
      fields.BAF_BANDS = 0
      fields.BAF_means = ''
      fields.BAF_sds   = ''
    else:
      fields.BAF_BANDS = c
      fields.BAF_means = ','.join('%.3f' % m for m in mu)
      fields.BAF_sds   = ','.join('%.3f' % s for s in sd)


def t_test(mu1,sd1,n1,mu2,sd2,n2):
  ss1   = sd1**2 / n1
  ss2   = sd2**2 / n2
  ssp   = ss1 + ss2
  t     = (mu1 - mu2) / ssp**0.5
  df    = ssp**2 // (ss1**2/(n1-1) + ss2**2/(n2-1))
  p     = 2*scipy.stats.distributions.t.cdf(-abs(t),df)

  return t,df,p


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('events',       help='Tablular or delimited file of events')
  parser.add_argument('index',        help='Index file created by cnv.index_gdats')

  parser.add_argument('--gcmodel',    metavar='GCM', help='GC model file to use')
  parser.add_argument('--gcmodeldir', metavar='DIR', help='GC models directory')
  parser.add_argument('--outdir',     metavar='DIR', help='Plot output directory', default='.')
  parser.add_argument('--outevents',  metavar='FILE', help='Output annotated events')
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
  parser.add_argument('--includesamples', metavar='FILE', action='append',
                                      help='List of samples to include, all others will be skipped')
  parser.add_argument('--excludesamples', metavar='FILE', action='append',
                                      help='List of samples to exclude, only samples not present will be kept')

  table_options(parser)

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)
  indexfile = GDATIndex(options.index)

  gcmodels  = {}
  chip_indices = {}

  add_fields = ['Probes','STATE','PMosaic','LRR_mean','LRR_sd','BAF_BANDS','BAF_means','BAF_sds']

  events    = table_reader(options.events)
  events    = cook_table(events,options)

  header    = next(events)
  extra     = [ h for h in add_fields if h not in header ]
  header   += extra
  blanks    = ['']*len(extra)

  Event     = recordtype('Event', header)._make
  events    = [ Event(e+blanks) for e in events ]

  events.sort(key=attrgetter(options.assayid,options.chromid))

  out       = table_writer(options.outevents) if options.outevents else None

  if out:
    out.writerow(header)

  transform = GenoTransform.from_object(options)

  if transform is not None:
    include  = transform.samples.include
    exclude  = transform.samples.exclude
  else:
    include  = exclude = None

  for assay,assay_events in groupby(events,key=attrgetter(options.assayid)):
    if include is not None and assay not in include:
      continue

    if exclude is not None and assay in exclude:
      continue

    assay_events = list(assay_events)

    locations  = list(indexfile.get(assay))

    if not locations:
      print 'ASSAY %s not found' % assay
      continue

    for gdat,offset in locations:
      gdatname        = '_'.join(os.path.splitext(os.path.basename(gdat.filename))[0].split('_')[:2])
      manifest        = gdat.attrs['ManifestName'].replace('.bpm','')

      chrom_indices,\
      autosome_mask   = chip_indices.get(manifest, (None,None))

      if chrom_indices is None:
        print 'Loading mapping information for %s...' % manifest
        chrom_indices = gdat.chromosome_index
        autosome_mask = valid_autosome_mask(gdat.snp_count,chrom_indices)

        chip_indices[manifest] = chrom_indices,autosome_mask

      print 'GDAT:ASSAY =',gdatname,assay

      assay_id,genos,lrr,baf = gdat.cnv_data(offset)
      normal_mask            = build_normal_mask(len(lrr),chrom_indices,assay_events,options)
      valid_mask             = (lrr>=-2)&(lrr<=2)
      valid_mask            &= autosome_mask

      base_mask              = normal_mask&valid_mask

      if not base_mask.sum():
        print '  ASSAY %s DOES NOT CONTAIN VALID LRR/BAF data' % assay
        continue

      lrr -= lrr[base_mask].mean()

      if gccorrect:
        gcmodel  = gcmodels.get(manifest)
        if gcmodel is None:
          print 'Loading GC/CpG model for %s...' % manifest
          filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
          gcmodels[manifest] = gcmodel = get_gcmodel(filename,chrom_indices)

        gcdesign,gcmask = gcmodel
        lrr             = gc_correct(lrr, gcdesign, gcmask&valid_mask)

      normal_lrr    = lrr[base_mask]
      normal_lrr_sd = normal_lrr.std(ddof=1)
      print '  ALL    Normal: probes=%7d LRR=%6.3f +- %.3f' % (len(normal_lrr),normal_lrr.mean(),
                                                                          1.96*normal_lrr_sd)

      aa_baf,aa_baf_max,\
      bb_baf,bb_baf_min = mirror_baf_distribution(baf,genos,base_mask)

      eventrecs = build_events(lrr, baf, chrom_indices, valid_mask, assay_events, options)
      eventrecs.sort(key=attrgetter('chrom'))

      for i,event in enumerate(eventrecs):
        if event.lrr is None:
          print '  EVENT %02d: Skipped... no valid probes.' % (i+1)
          continue

        event_baf = strip_event_baf(event.baf,aa_baf_max,bb_baf_min)
        #event_baf = augment_event_baf(event.baf,aa_baf,aa_baf_max,bb_baf,bb_baf_min)
        event.baf_model = fit_gmm_baf(event_baf,max_baf=2)

        infer_mosaic_state(event)
        update_event_fields(event)

        if out:
          out.writerow(event.fields)

        print '  EVENT %2d: chr%-3s probes=%d LRR=%6.3f +- %.3f STATE=%s %%Mosaic=%s' \
                    % (i+1, event.chrom,len(event.lrr),event.lrr.mean(),1.96*event.lrr.std(ddof=1),
                            event.state, ('%.2f' % (100*event.pmosaic) if event.pmosaic is not None else '?'))

      print

      for chrom,chrom_events in groupby(eventrecs,key=attrgetter('chrom')):
        chrom        = norm_chromosome(chrom)
        chrom_events = list(chrom_events)
        pos,index    = chrom_indices[chrom]

        #chrom_valid  = valid_mask[index]
        #chrom_normal = normal_mask[index]
        chrom_genos  = genos[index]
        chrom_baf    = baf[index]
        chrom_lrr    = lrr[index]

        states    = ','.join('%s:%s%%' % (e.state,int(100*e.pmosaic) if e.pmosaic else '?') for e in chrom_events)
        variables = chrom_events[0].fields._asdict()
        variables.update(GDAT=gdatname,ASSAY=assay,CHROMOSOME=chrom,STATES=states)

        plotname = options.template.format(**variables)
        plotname = '%s/%s' % (options.outdir,plotname)
        title    = options.title.format(**variables)

        plot_chromosome(plotname, pos, chrom_lrr, chrom_baf, genos=chrom_genos,
                        title=title, events=chrom_events)


if __name__=='__main__':
  main()
