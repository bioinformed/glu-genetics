# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Perform quantile normalization and re-estimate LRR and BAF'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import numpy    as np
import numpy.ma as ma

from   glu.lib.fileutils import table_reader

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


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdat',                          help='Input GDAT file')

  parser.add_argument('--gcmodel',     metavar='GCM',  help='GC model file to use')
  parser.add_argument('--gcmodeldir',  metavar='DIR',  help='GC models directory')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)

  gdat      = GDATFile(options.gdat,'r+')
  qnorm     = True

  if gccorrect:
    manifest = gdat.attrs['ManifestName'].replace('.bpm','')
    print 'Loading GC/CpG model for %s...' % manifest
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    gcdesign,gcmask = get_gcmodel(filename, gdat.chromosome_index)

  X         = gdat['X']
  Y         = gdat['Y']
  LRR       = gdat['LRR']
  BAF       = gdat['BAF']

  x_scale   = X.attrs['SCALE']
  x_nan     = X.attrs['NAN']

  y_scale   = Y.attrs['SCALE']
  y_nan     = Y.attrs['NAN']

  lrr_scale = LRR.attrs['SCALE']
  lrr_nan   = LRR.attrs['NAN']

  baf_scale = BAF.attrs['SCALE']
  baf_nan   = BAF.attrs['NAN']

  genos     = gdat['Genotype']

  n,s       = X.shape

  print 'SNPs=%d, samples=%d' % (s,n)

  np.set_printoptions(linewidth=120,precision=2,suppress=True)
  np.seterr(all='ignore')

  n_AA        = np.zeros(s, dtype=int  )
  r_AA        = np.zeros(s, dtype=float)
  t_AA        = np.zeros(s, dtype=float)
  n_AB        = np.zeros(s, dtype=int  )
  r_AB        = np.zeros(s, dtype=float)
  t_AB        = np.zeros(s, dtype=float)
  n_BB        = np.zeros(s, dtype=int  )
  r_BB        = np.zeros(s, dtype=float)
  t_BB        = np.zeros(s, dtype=float)

  for i in xrange(n):
    print 'Sample %5d / %d' % (i+1,n)

    x         = gdat_decode(X[i], x_scale, x_nan)
    y         = gdat_decode(Y[i], y_scale, y_nan)

    if qnorm:
      x,y     = compute_qn(x,y)

    r         = x+y
    t         = (2/np.pi)*np.arctan2(y,x)

    genosi    = genos[i]

    mask      = np.isfinite(r)&np.isfinite(t)

    update_centers(r,t,'AA',r_AA,t_AA,n_AA,genosi,mask)
    update_centers(r,t,'AB',r_AB,t_AB,n_AB,genosi,mask)
    update_centers(r,t,'BB',r_BB,t_BB,n_BB,genosi,mask)

  finalize_centers(r_AA,t_AA,n_AA)
  finalize_centers(r_AB,t_AB,n_AB)
  finalize_centers(r_BB,t_BB,n_BB)

  bad         = t_AA>=t_AB
  bad        |= t_AB>=t_BB
  bad        |= t_AA>=t_BB

  t_AA[bad]   = np.nan
  t_AB[bad]   = np.nan
  t_BB[bad]   = np.nan

  create_gdat_qn(gdat,s,n)

  LRR_QN = gdat['LRR_QN']
  BAF_QN = gdat['BAF_QN']

  for i in xrange(n):
    print 'Sample %5d / %d' % (i+1,n)

    x         = gdat_decode(X[i], x_scale, x_nan)
    y         = gdat_decode(Y[i], y_scale, y_nan)

    if qnorm:
      x,y     = compute_qn(x,y)

    r         = x+y
    t         = (2/np.pi)*np.arctan2(y,x)

    lrr,baf   = compute_lrr_baf(t,r,r_AA,r_AB,r_BB,t_AA,t_AB,t_BB)

    if gccorrect:
      lrr     = gc_correct(lrr, gcdesign, gcmask, minval=-2, maxval=2)

    LRR_QN[i] = gdat_encode_f(lrr, LRR_SCALE, LRR_MIN, LRR_MAX, LRR_NAN)
    BAF_QN[i] = gdat_encode_f(baf, BAF_SCALE, BAF_MIN, BAF_MAX, BAF_NAN)


if __name__ == '__main__':
  main()
