import numpy as np

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec

from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib.ticker       import NullFormatter


def plot_chromosome(filename, startattr, stopattr, pos, lrr, baf, title=None, events=None):
  nullfmt  = NullFormatter()

  plt.clf()

  pos = pos/1000000.
  fig = plt.figure(1, figsize=(10.5,5.5))

  if title:
    fig.suptitle(title)

  gs = gridspec.GridSpec(2,5,width_ratios=[20,0.5,2,0.001,2],height_ratios=[1,19])

  main = plt.subplot(gs[1,0])
  main.axis( [pos.min(),pos.max(),-2,2] )

  sec = main.twinx()
  sec.axis( [pos.min(),pos.max(),0,1] )

  main.scatter(pos, lrr, c='black', s=2, linewidths=0)
  sec.scatter(pos,  baf, c='red',   s=2, linewidths=0)

  main.set_xlabel('Chromosome Location (Mbps)')
  main.set_ylabel('Log Intensity Ratio (LRR)')
  sec.set_ylabel('B-Allele Frequency (BAF)')

  for tl in main.get_yticklabels():
    tl.set_color('black')

  for tl in sec.get_yticklabels():
    tl.set_color('red')

  binwidth = 0.02
  bins     = np.arange(-2, 2., binwidth)

  all_hist = plt.subplot(gs[1,2])
  all_hist.xaxis.set_major_formatter(nullfmt)
  all_hist.yaxis.set_major_formatter(nullfmt)

  mask = (baf>0.01)&(baf<0.99)
  hbaf = (baf[mask]-0.5)*4

  all_hist.hist(lrr,  bins=bins, color='black', label='LRR',
                      normed=True, orientation='horizontal', histtype='step')
  all_hist.hist(hbaf, bins=bins, color='red', label='BAF',
                      normed=True, orientation='horizontal', histtype='step')
  all_hist.set_ylim(main.get_ylim())

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

      #print '  SEG SIZE:',mask.sum()
      if mask.sum()>1:
        lrr_mean = lrr[mask].mean()
        #print '    SEG:',chromosome,start,stop,lrr_mean
        main.plot([start,stop],[lrr_mean,lrr_mean],marker='o', color='blue', linewidth=5)

      evsegs.plot([start,stop],[0,0], marker='o', color='blue', linewidth=5)
      evsegs.annotate('%d' % i, xy=(mid,0), xycoords='data',
                                xytext=(1,0), horizontalalignment='center')

    seg_hist = plt.subplot(gs[1,4])
    seg_hist.xaxis.set_major_formatter(nullfmt)
    seg_hist.yaxis.set_major_formatter(nullfmt)

    mask = (baf>0.01)&(baf<0.99)&event_mask

    seg_hist.hist(lrr[event_mask],  bins=bins, color='blue', label='LRR',
                                    normed=True, orientation='horizontal', histtype='step')
    seg_hist.hist((baf[mask]-0.5)*4,bins=bins, color='red', label='BAF',
                                    normed=True, orientation='horizontal', histtype='step')
    seg_hist.set_ylim(main.get_ylim())

  plt.show()
  plt.savefig(filename)
