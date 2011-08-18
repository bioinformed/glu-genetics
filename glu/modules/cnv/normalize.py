# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Perform quantile normalization and re-estimate LRR and BAF'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   math                 import modf

import numpy as np

from   glu.lib.fileutils    import table_reader

from   glu.lib.glm          import Linear

from   glu.modules.cnv.gdat import GDATFile, gdat_encode_f, gdat_decode, get_gcmodel, gc_correct


BAF_TYPE   = np.int16
BAF_SCALE  =  10000
BAF_NAN    = np.iinfo(BAF_TYPE).min
BAF_MIN    = np.iinfo(BAF_TYPE).min+1
BAF_MAX    = np.iinfo(BAF_TYPE).max

LRR_TYPE   = np.int32
LRR_SCALE  = 100000
LRR_NAN    = np.iinfo(LRR_TYPE).min
LRR_MIN    = np.iinfo(LRR_TYPE).min+1
LRR_MAX    = np.iinfo(LRR_TYPE).max


def boxcox(x,l=0,d=0,base=np.e):
  '''
  Shifted Box-Cox Power Transformation

  l    : lambda         (default=0)
  d    : delta          (defailt=0)
  base : logarithm base (default=e)

  See: http://en.wikipedia.org/wiki/Power_transform
       http://www.jstor.org/pss/2984418
  '''
  assert l>=0

  x = x.astype(float)-d

  if l==0:
    y = np.log(x)
  else:
    y = (x**l-1)/l

  if base!=np.e:
    y /= np.log(b)

  return y


def gmpt(x,l=0,d=0,base=np.e):
  '''
  General Modulus Power Transformation

  l    : lambda         (default=0)
  d    : delta          (defailt=0)
  base : logarithm base (default=e)

  See: http://en.wikipedia.org/wiki/Power_transform
       http://www.informaworld.com/smpp/content~db=all~content=a780040516

  '''
  assert l>=0

  x = np.asarray(x, dtype=float)

  if l==0:
    y = np.sign(x)*np.log(np.abs(x)+1)
  else:
    y = np.sign(x)*((np.abs(x)+1)**l-1)/l

  if base!=np.e:
    y /= np.log(base)

  return y


def quantile_normalize(x,y):
  mask       = np.isfinite(x)&np.isfinite(y)
  xx         = x[mask]
  yy         = y[mask]

  xorder     = xx.argsort()
  yorder     = yy.argsort()

  values     = (xx[xorder]+yy[yorder])/2.0

  xx[xorder] = values
  yy[yorder] = values

  x[mask]    = xx
  y[mask]    = yy

  return x,y


def outlier_mask(data):
  mask    = np.isfinite(data)
  valid   = data[mask]
  if len(valid)<20:
    return mask|(~mask)
  valid.sort()
  n       = valid.shape[0]
  vmin    = np.min(valid[ 5],valid[int(n*0.01)])
  vmax    = np.max(valid[-5],valid[int(n*0.99)])
  outlier = (~mask)|(data<vmin)|(data>vmax)
  return outlier


def trimmed_mean(data, p):
  mask    = np.isfinite(data)
  valid   = data[mask]
  valid.sort()
  n       = valid.shape[0]
  lower   = int(p*n)
  upper   = n-lower
  return valid[lower:upper].mean(axis=0)


def compute_lrr_baf(t,r,r_AA,r_AB,r_BB,t_AA,t_AB,t_BB):
  # Set assays with no observed hets to use average of two
  # homozygotes for r and theta
  mask    = np.isfinite(t_AA)&np.isfinite(t_BB)&(~np.isfinite(t_AB))
  t_AB[mask] = (t_AA[mask]+t_BB[mask])/2.0
  r_AB[mask] = (r_AA[mask]+r_BB[mask])/2.0

  t0      = (t<=t_AA)
  t1      = (t_AA< t)&(t<t_AB)
  t2      = (t_AB<=t)&(t<t_BB)
  t3      = (t>=t_BB)

  rex1    = (r_AB-r_AA)*(t-t_AA)/(t_AB-t_AA)+r_AA
  rex2    = (r_BB-r_AB)*(t-t_AB)/(t_BB-t_AB)+r_AB

  baf1    =    (t-t_AA)/(t_AB-t_AA)  / 2
  baf2    = (1+(t-t_AB)/(t_BB-t_AB)) / 2

  rex     = np.empty_like(r)
  rex.fill(np.nan)

  rex[t1] = rex1[t1]
  rex[t2] = rex2[t2]
  rex[t0] = r_AA[t0]
  rex[t3] = r_BB[t3]

  if 0:
    print '  t_AA',t_AA[:10]
    print '  t_AB',t_AB[:10]
    print '  t_BB',t_BB[:10]
    print '  r_AA',r_AA[:10]
    print '  r_AB',r_AB[:10]
    print '  r_BB',r_BB[:10]
    print '     t',t[:10]
    print '     r',r[:10]
    print '   rex',rex[:10]
    print '   lrr',lrr[:10]

  rex    += 1e-6
  lrr     = np.log2(r/rex)

  baf     = np.empty_like(t)
  baf.fill(np.nan)

  baf[t1] = baf1[t1]
  baf[t2] = baf2[t2]
  baf[t0] = 0
  baf[t3] = 1

  return lrr,baf


def compute_qn(x,y,threshold=1.5):
  xmax    = threshold*x
  ymax    = threshold*y

  xqn,yqn = quantile_normalize(x,y)

  x       = np.min([xqn,xmax],axis=0)
  y       = np.min([yqn,ymax],axis=0)

  return x,y


def update_centers(r,t,g,r_g,t_g,n_g,genos,mask):
  mask       = (genos==g)&mask
  n_g       +=   mask
  r_g[mask] += r[mask]
  t_g[mask] += t[mask]


def finalize_centers(r_g,t_g,n_g):
  r_g       /= n_g
  t_g       /= n_g


def create_gdat_qn(gdatobject,s,n):
  comp         = dict(compression='gzip',compression_opts=5)
  chunks       = (1,s)
  shape        = (n,s)
  shuffle      = False

  gdat         = gdatobject.gdat
  BAF_QN       = gdat.require_dataset('BAF_QN', shape, BAF_TYPE,
                                      maxshape=shape,chunks=chunks,shuffle=shuffle,
                                      fillvalue=BAF_NAN,**comp)
  BAF_QN.attrs['SCALE'] = BAF_SCALE
  BAF_QN.attrs['NAN']   = BAF_NAN
  BAF_QN.attrs['MIN']   = BAF_MIN
  BAF_QN.attrs['MAX']   = BAF_MAX

  LRR_QN       = gdat.require_dataset('LRR_QN', shape, LRR_TYPE,
                                      maxshape=shape,chunks=chunks,shuffle=shuffle,
                                      fillvalue=LRR_NAN,**comp)

  LRR_QN.attrs['SCALE'] = LRR_SCALE
  LRR_QN.attrs['NAN']   = LRR_NAN
  LRR_QN.attrs['MIN']   = LRR_MIN
  LRR_QN.attrs['MAX']   = LRR_MAX


def regress_center(p, design, dmask, genos, p_AA, p_AB, p_BB, thin=None):
  AA           =  genos=='AA'
  AB           =  genos=='AB'
  BB           =  genos=='BB'

  valid        = (genos!='  ')
  valid       &= np.isfinite(p)
  valid       &= dmask

  design[:,1]  = p

  expected     = np.empty_like(p)
  expected[AA] = p_AA[AA]
  expected[AB] = p_AB[AB]
  expected[BB] = p_BB[BB]

  valid       &= np.isfinite(expected)

  expected     = expected[valid].reshape(-1,1)
  design       = design[valid]

  if thin:
    expected   = expected[::thin]
    design     = design[::thin]

  lm = Linear(expected, design)

  lm.fit()

  ss_t  = np.var(lm.y,ddof=1)
  r2    = 1 - lm.ss/ss_t
  beta  = lm.beta.reshape(-1)

  print '  Regress to center: r2=%.2f beta=%s' % (r2,beta)

  return beta


def adjust_center(p,beta,design,dmask):
  p_adj        = np.empty_like(p)
  p_adj.fill(np.nan)

  valid        = np.isfinite(p)
  valid       &= dmask
  design[:,1]  = p
  p_adj[valid] = np.dot(design[valid],beta.reshape(-1,1))

  return p_adj


def quantile(data,k):
  mask = np.isfinite(data)
  data = data[mask]

  data.sort()

  n    = len(data)-1
  f,w  = modf(n*k)
  w    = int(w)

  if f < 1e-10:
    q = data[w]
  else:
    q = (1-f)*data[w] + f*data[w+1]

  return q


def gdat_decoder(table):
  if 'SCALE' not in table.attrs:
    return rows

  def _decoder(table):
    scale = table.attrs['SCALE']
    nan   = table.attrs['NAN']
    for row in table:
      yield gdat_decode(row, scale, nan)

  return _decoder(table)


def block_iter(table,chunksize):
  start = 0
  last  = len(table)

  while start<last:
    stop = min(last,start+chunksize)
    chunk = table[start:stop]
    for row in chunk:
      yield row
    start = stop


def parallel_gdat_iter(*tables):
  iters = [ gdat_decoder(table,block_iter(table,table.chunks[0])) for table in tables ]
  return izip(*iters)


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdat',                          help='Input GDAT file')

  parser.add_argument('--gcmodel',     metavar='GCM',  help='GC model file to use')
  parser.add_argument('--gcmodeldir',  metavar='DIR',  help='GC models directory')
  parser.add_argument('--maxmissing', metavar='RATE', default=0.02, type=float,
                          help='Maximum missing rate for assays to estimate cluster centers (default=0.02)')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)

  gdat      = GDATFile(options.gdat,'r+')
  n         = gdat.sample_count
  s         = gdat.snp_count
  qnorm     = True

  if gccorrect:
    manifest = gdat.attrs['ManifestName'].replace('.bpm','')
    print 'Loading GC/CpG model for %s...' % manifest
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    design,dmask = get_gcmodel(filename, extra_terms=1)
  else:
    design  = np.ones( (s,2), dtype=float )
    dmask   = design[0]==1

  X         = gdat['X']
  Y         = gdat['Y']
  GC        = gdat['GC']
  LRR       = gdat['LRR']
  BAF       = gdat['BAF']
  genotypes = gdat['Genotype']

  print 'SNPs=%d, samples=%d' % (s,n)

  np.set_printoptions(linewidth=120,precision=2,suppress=True)
  np.seterr(all='ignore')

  print 'PASS 1: Estimate intensity distribution...'

  n_AA        = np.zeros(s, dtype=int  )
  x_AA        = np.zeros(s, dtype=float)
  y_AA        = np.zeros(s, dtype=float)
  n_AB        = np.zeros(s, dtype=int  )
  x_AB        = np.zeros(s, dtype=float)
  y_AB        = np.zeros(s, dtype=float)
  n_BB        = np.zeros(s, dtype=int  )
  x_BB        = np.zeros(s, dtype=float)
  y_BB        = np.zeros(s, dtype=float)

  skipped     = 0

  for i,(x,y,genos) in enumerate(parallel_gdat_iter(X,Y,genotypes)):
    print 'Sample %5d / %d' % (i+1,n)

    missing   = (genos=='  ').sum() / s

    if missing > options.maxmissing:
      #print '  ... skipping due to missing rate %.2f%% > %.2f%%' % (missing*100,options.maxmissing*100)
      skipped += 1
      continue

    mask      = np.isfinite(x)&np.isfinite(y)

    update_centers(x,y,'AA',x_AA,y_AA,n_AA,genos,mask)
    update_centers(x,y,'AB',x_AB,y_AB,n_AB,genos,mask)
    update_centers(x,y,'BB',x_BB,y_BB,n_BB,genos,mask)

  finalize_centers(x_AA,y_AA,n_AA)
  finalize_centers(x_AB,y_AB,n_AB)
  finalize_centers(x_BB,y_BB,n_BB)

  if skipped:
    print '  Skipped %d samples (%.2f%% ) for missing rate > %.2f%%' % (skipped,skipped/n*100,options.maxmissing*100)

  print 'PASS 2: Quantile normalization...'

  beta_xs     = []
  beta_ys     = []

  n_AA        = np.zeros(s, dtype=int  )
  r_AA        = np.zeros(s, dtype=float)
  t_AA        = np.zeros(s, dtype=float)
  n_AB        = np.zeros(s, dtype=int  )
  r_AB        = np.zeros(s, dtype=float)
  t_AB        = np.zeros(s, dtype=float)
  n_BB        = np.zeros(s, dtype=int  )
  r_BB        = np.zeros(s, dtype=float)
  t_BB        = np.zeros(s, dtype=float)

  for i,(x,y,genos,gqual) in enumerate(parallel_gdat_iter(X,Y,genotypes,GC)):
    print 'Sample %5d / %d' % (i+1,n)

    min_qual  = quantile(gqual,0.5)

    qmask     = gqual>=min_qual
    qmask    &= dmask

    beta_x    = regress_center(x,design,dmask,genos,x_AA,x_AB,x_BB,thin=5)
    beta_y    = regress_center(y,design,dmask,genos,y_AA,y_AB,y_BB,thin=5)

    beta_xs.append(beta_x)
    beta_ys.append(beta_y)

    missing   = (genos=='  ').sum() / s
    if missing > options.maxmissing:
      continue

    x         = adjust_center(x,beta_x,design,dmask)
    y         = adjust_center(y,beta_y,design,dmask)

    if qnorm:
      x,y     = compute_qn(x,y)

    r         = x+y
    t         = (2/np.pi)*np.arctan2(y,x)

    mask      = np.isfinite(r)&np.isfinite(t)

    update_centers(r,t,'AA',r_AA,t_AA,n_AA,genos,mask)
    update_centers(r,t,'AB',r_AB,t_AB,n_AB,genos,mask)
    update_centers(r,t,'BB',r_BB,t_BB,n_BB,genos,mask)

  finalize_centers(r_AA,t_AA,n_AA)
  finalize_centers(r_AB,t_AB,n_AB)
  finalize_centers(r_BB,t_BB,n_BB)

  print 'PASS 3: Updating LRR and BAF...'

  bad         = t_AA>=t_AB
  bad        |= t_AB>=t_BB
  bad        |= t_AA>=t_BB

  t_AA[bad]   = np.nan
  t_AB[bad]   = np.nan
  t_BB[bad]   = np.nan

  create_gdat_qn(gdat,s,n)

  LRR_QN = gdat['LRR_QN']
  BAF_QN = gdat['BAF_QN']

  for i,(x,y,genos) in enumerate(parallel_gdat_iter(X,Y,genotypes)):
    print 'Sample %5d / %d' % (i+1,n)

    x         = adjust_center(x,beta_xs[i],design,dmask)
    y         = adjust_center(y,beta_ys[i],design,dmask)

    if qnorm:
      x,y     = compute_qn(x,y)

    r         = x+y
    t         = (2/np.pi)*np.arctan2(y,x)

    beta_r    = regress_center(r,design,dmask,genos,r_AA,r_AB,r_BB,thin=5)
    r         = adjust_center(r,beta_r,design,dmask)

    lrr,baf   = compute_lrr_baf(t,r,r_AA,r_AB,r_BB,t_AA,t_AB,t_BB)

    LRR_QN[i] = gdat_encode_f(lrr, LRR_SCALE, LRR_MIN, LRR_MAX, LRR_NAN)
    BAF_QN[i] = gdat_encode_f(baf, BAF_SCALE, BAF_MIN, BAF_MAX, BAF_NAN)


if __name__ == '__main__':
  main()
