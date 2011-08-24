from __future__ import division

import numpy as np

import matplotlib.pyplot   as plt
import matplotlib          as mpl
import matplotlib.gridspec as gridspec

from   matplotlib.ticker       import NullFormatter

from   scipy.stats             import gaussian_kde

from   glu.lib.recordtype      import recordtype


def plot_baf_bars(seg_hist,c,mu,sd):
  colors = ['lightgreen']*c

  for i in range(c):
    left  = mu[i]-sd[i]
    right = mu[i]+sd[i]
    width = right-left
    peak  = mpl.patches.Rectangle( (0,left*4-2), 2, width*4, color=colors[i], zorder=1)
    peak.set_clip_box(seg_hist.bbox)
    peak.set_alpha(0.5)
    seg_hist.add_patch(peak)
    seg_hist.axhline(y=mu[i]*4-2, color='black', zorder=100)


def plot_hist(hist,bins,lrr,baf,gmm=None):
  nullfmt  = NullFormatter()
  hist.xaxis.set_major_formatter(nullfmt)
  hist.yaxis.set_major_formatter(nullfmt)

  if lrr is not None and len(lrr)>100:
    hist.hist(lrr,  bins=bins, color='black', label='LRR',
                    normed=True, orientation='horizontal', histtype='step')

  if baf is None or lrr is None:
    return

  if gmm is not None:
    c,mu,sd,ws = gmm
    plot_baf_bars(hist,c,mu,sd)

  mask = (baf>0.01)&(baf<0.99)
  if mask.sum()>100:
    hbaf  = (baf[mask]-0.50)*4
    hist.hist(hbaf, bins=bins, color='red', label='BAF',
                    normed=True, orientation='horizontal', histtype='step')


def plot_chromosome(filename, pos, lrr, baf, genos=None, title=None, events=None):
  valid           = np.isfinite(lrr)
  lrr             = lrr[valid]
  baf             = baf[valid]
  normal_lrr      = lrr
  normal_baf      = baf
  pos             = pos[valid]/1000000.
  genos           = genos[valid] if genos is not None else None

  if events:
    event_mask    = np.zeros_like(lrr,dtype=bool)

    for i,event in enumerate(events):
      start       = event.start / 1000000
      stop        = event.stop  / 1000000
      event_mask |= (pos>=start)&(pos<=stop)

    normal_lrr = normal_lrr[~event_mask]
    normal_baf = normal_baf[~event_mask]

  np.set_printoptions(linewidth=120,precision=2,suppress=True)
  np.seterr(all='ignore')

  nullfmt  = NullFormatter()

  plt.clf()

  fig = plt.figure(1, figsize=(10.5,5.5))

  if title:
    fig.suptitle(title)

  fig.subplots_adjust(left=0.0001, right=0.9999, wspace=0.0001)

  histograms = (len(normal_lrr)>0) + len(events or [])
  widths = [2.5,20,2.25,2]+[0.1,2]*(histograms-1)
  nsub   = 4+2*(histograms-1)
  gs = gridspec.GridSpec(2,nsub,width_ratios=widths,height_ratios=[1,19])

  main = plt.subplot(gs[1,1])
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

  main.set_xlabel('Chromosome Location (Mbps)')
  main.set_ylabel('Log Intensity Ratio (LRR)')
  sec.set_ylabel('B-Allele Frequency (BAF)')

  for tl in main.get_yticklabels():
    tl.set_color('black')

  for tl in sec.get_yticklabels():
    tl.set_color('red')

  binwidth = 0.02
  bins     = np.arange(-2, 2, binwidth)

  if events:
    evsegs = plt.subplot(gs[0,1])
    evsegs.set_xlim(main.get_xlim())
    evsegs.set_ylim( (0,1) )
    evsegs.xaxis.set_major_formatter(nullfmt)
    evsegs.yaxis.set_major_formatter(nullfmt)

    for i,event in enumerate(events):
      start = event.start / 1000000
      stop  = event.stop  / 1000000
      mid   = (event.start+event.stop)/2.0

      if event.lrr is not None and len(event.lrr):
        lrr_mean  = event.lrr.mean()
        lrr_bmean = (lrr_mean+2)/4

        sec.plot([start,stop],[lrr_bmean,lrr_bmean],
                 marker='o', color='blue', linewidth=2, zorder=100)

        sec.axvline(start, color='black', zorder=1, alpha=0.3, linewidth=2)
        sec.axvline(stop,  color='black', zorder=1, alpha=0.3, linewidth=2)

      evsegs.plot([start,stop],[0,0], marker='o', color='blue', linewidth=5)
      evsegs.annotate('%d' % i, xy=(mid,0), xycoords='data',
                                xytext=(1,0), horizontalalignment='center')


  h = 3
  if len(normal_lrr):
    norm_hist = plt.subplot(gs[1,h])
    h        += 2

    norm_hist.set_xlabel('Norm')
    norm_hist.set_ylim(main.get_ylim())
    plot_hist(norm_hist,bins,normal_lrr,normal_baf)

  for i,event in enumerate(events):
    seg_hist  = plt.subplot(gs[1,h])
    h        += 2

    seg_hist.set_ylim(main.get_ylim())
    seg_hist.set_xlabel('Seg%d' % (i+1))

    plot_hist(seg_hist,bins,event.lrr,event.baf,event.baf_model)

  plt.show()
  plt.savefig(filename)
