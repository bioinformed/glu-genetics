from __future__ import division

import numpy as np

import matplotlib.pyplot   as plt
import matplotlib          as mpl
import matplotlib.gridspec as gridspec

from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib.ticker       import NullFormatter

from   scipy.stats             import gaussian_kde
from   scikits.learn.mixture   import GMM, DPGMM, VBGMM

from   glu.modules.cnv.gmm     import gaussian_mixture


def fit_gmm(components,x):
  x = np.asanyarray(x).reshape(-1,1)
  #gmm_baf = VBGMM(components,verbose=True,min_covar=0.01)
  gmm_baf = GMM(components)
  gmm_baf.fit(x)

  means = np.array(gmm_baf.means,  dtype=float).reshape(-1)
  order = means.argsort()
  means = means[order]
  sds   = np.array(gmm_baf.covars,  dtype=float).reshape(-1)[order]**0.5
  ws    = np.array(gmm_baf.weights, dtype=float).reshape(-1)[order]
  logL  = gmm_baf.score(x)

  return logL,means,sds,ws


def fit_gmm2(k,x):
  w  = np.ones(k)/k
  mu = np.array([ i/(k-1) for i in range(k) ], dtype=float)
  sd = np.array([0.05]+[0.2]*(k-2)+[0.05], dtype=float)

  logL,w,mu,sd = gaussian_mixture(k,x, w=w, mu=mu, sd=sd,min_w=0.15)

  return logL,mu,sd,w


def plot_chromosome(filename, pos, lrr, baf, genos=None, title=None,
                              events=None, startattr=None, stopattr=None):

  np.set_printoptions(linewidth=120,precision=2,suppress=True)
  np.seterr(all='ignore')

  nullfmt  = NullFormatter()

  plt.clf()

  pos = pos/1000000.
  fig = plt.figure(1, figsize=(10.5,5.5))

  if title:
    fig.suptitle(title)

  gs = gridspec.GridSpec(2,5,width_ratios=[20,0.5,2,0.001,2],height_ratios=[1,19])

  main = plt.subplot(gs[1,0])
  main.axis( [pos.min(),pos.max(),-2,2] )
  main.tick_params(axis='x', which='major', direction='out', length=4, width=1, color='black')
  main.locator_params(axis='x', tight=True, nbins=15)

  sec = main.twinx()
  sec.axis( [pos.min(),pos.max(),0,1] )

  main.scatter(pos,  lrr, c='black', s=2, linewidths=0, zorder=4)

  if genos is None:
    sec.scatter(pos, baf, c='red',   s=2, linewidths=0, zorder=4)
  else:
    hets = genos=='AB'
    miss = genos=='  '
    homs = (~miss)&(~hets)

    sec.scatter(pos[homs], baf[homs], c='darkred', marker='o', s=2, linewidths=0, alpha=0.8, zorder=4)
    sec.scatter(pos[hets], baf[hets], c='red',     marker='o', s=2, linewidths=0, alpha=0.5, zorder=5)
    sec.scatter(pos[miss], baf[miss], c='red',     marker='o', s=2, linewidths=0, alpha=0.8, zorder=6)
    #sec.scatter(pos[hets], baf[hets], c='red',     marker='+', s=2, linewidths=1, alpha=0.8, zorder=5, edgecolor='red')
    #sec.scatter(pos[miss], baf[miss], c='white', marker='o', s=2, linewidths=1, alpha=0.8, zorder=6, edgecolor='red')

  main.set_xlabel('Chromosome Location (Mbps)')
  main.set_ylabel('Log Intensity Ratio (LRR)')
  sec.set_ylabel('B-Allele Frequency (BAF)')

  for tl in main.get_yticklabels():
    tl.set_color('black')

  for tl in sec.get_yticklabels():
    tl.set_color('red')

  if events:
    evsegs = plt.subplot(gs[0,0])
    evsegs.set_xlim(main.get_xlim())
    evsegs.set_ylim( (0,1) )
    evsegs.xaxis.set_major_formatter(nullfmt)
    evsegs.yaxis.set_major_formatter(nullfmt)

    event_mask = (lrr!=lrr)

    for i,event in enumerate(events):
      start = int(getattr(event,startattr))/1000000.
      stop  = int(getattr(event,stopattr ))/1000000.
      mid   = (start+stop)/2.0

      mask  = (pos>=start-1)&(pos<=stop)&np.isfinite(lrr)
      event_mask |= mask

      if mask.sum()>1:
        lrr_mean  = lrr[mask].mean()
        lrr_bmean = (lrr_mean+2)/4
        sec.plot([start,stop],[lrr_bmean,lrr_bmean],marker='o', color='blue', linewidth=2, zorder=100)
        sec.axvline(start, color='black', zorder=1, alpha=0.3, linewidth=2)
        sec.axvline(stop,  color='black', zorder=1, alpha=0.3, linewidth=2)

      evsegs.plot([start,stop],[0,0], marker='o', color='blue', linewidth=5)
      evsegs.annotate('%d' % i, xy=(mid,0), xycoords='data',
                                xytext=(1,0), horizontalalignment='center')

  binwidth = 0.02
  bins     = np.arange(-2, 2., binwidth)

  all_hist = plt.subplot(gs[1,2])
  all_hist.xaxis.set_major_formatter(nullfmt)
  all_hist.yaxis.set_major_formatter(nullfmt)
  all_hist.set_xlabel('Normal')

  if events:
    mask    &= event_mask
    all_lrr  = lrr[~event_mask]
    all_baf  = baf[~event_mask]
  else:
    all_lrr  = lrr
    all_baf  = baf

  mask     = (all_baf>0.01)&(all_baf<0.99)
  if mask.sum()>5:
    all_hbaf = (all_baf[mask]-0.5)*4
    all_hist.hist(all_lrr,  bins=bins, color='black', label='LRR',
                            normed=True, orientation='horizontal', histtype='step')
    all_hist.hist(all_hbaf, bins=bins, color='red', label='BAF',
                            normed=True, orientation='horizontal', histtype='step')

  #if genos is not None:
  #  bdev_homs  = np.min([baf,1-baf], axis=0)
  #  bdev_hets  = np.abs(baf-0.5)
  #  bdev_miss  = np.min([bdev_homs,bdev_hets],axis=0)
  #
  #  bdev       = np.empty_like(baf)
  #
  #  bdev[homs] = bdev_homs[homs]
  #  bdev[hets] = bdev_hets[hets]
  #  bdev[miss] = bdev_miss[miss]
  #
  #  #mask = (bdev>0.01)
  #  #all_hist.hist( (bdev[mask]-0.5)*4, bins=bins, color='yellow', label='Bdev',
  #  #                             normed=True, orientation='horizontal', histtype='step')

  all_hist.set_ylim(main.get_ylim())

  if events:
    seg_hist   = plt.subplot(gs[1,4])
    seg_hist.xaxis.set_major_formatter(nullfmt)
    seg_hist.yaxis.set_major_formatter(nullfmt)
    seg_hist.set_xlabel('Events')

    #bdev_mask  = (bdev>0.02)&event_mask

    #event_bdev = (bdev[bdev_mask]-0.25)*8

    #gmm_baf = DPGMM(4,alpha=0.00001)
    #gmm_baf = VBGMM(5)

    if 0:
      event_baf = baf[event_mask]
      fits      = [ fit_gmm2(c,event_baf) for c in range(2,5) ]
      AIC       = np.array([ 2*(c-f[0].sum()) for c,f in enumerate(fits,2) ], dtype=float)
      best      = AIC.argmin()

      logL,means,sds,ws = fits[best]
      c                 = best+2
      print '  COMPONENTS:',c
      print '        logL:',logL.sum()
      print '         AIC:',2*(c-logL.sum())
      print '       MEANS:',means
      print '          SD:',sds
      print '     WEIGHTS:',ws
      print

      colors = ['yellow','lightgreen']*5

      for i in range(c):
        left  = means[i]-sds[i]
        right = means[i]+sds[i]
        width = right-left
        peak  = mpl.patches.Rectangle( (0,left*4-2), 2, width*4, color=colors[i], zorder=1)
        peak.set_clip_box(seg_hist.bbox)
        peak.set_alpha(0.5)
        seg_hist.add_patch(peak)

    baf_mask   = (baf >0.01)&(baf<0.99)&event_mask

    if event_mask.sum()>5:
      event_lrr  = lrr[event_mask ]
      event_baf  = ( baf[baf_mask ]-0.50)*4

      seg_hist.hist(event_lrr,  bins=bins, color='blue', label='LRR',
                                normed=True, orientation='horizontal', histtype='step')
      seg_hist.hist(event_baf,  bins=bins, color='red', label='BAF',
                                normed=True, orientation='horizontal', histtype='step')
      #seg_hist.hist(event_bdev, bins=bins, color='green', label='Bdev',
      #                          normed=True, orientation='horizontal', histtype='step')

      seg_hist.set_ylim(main.get_ylim())

  plt.show()
  plt.savefig(filename)

