# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Build GC/CpG model file (GCM) from an Illumina manifest file (BPM) and an indexed reference genome (FASTA)'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

from   itertools   import groupby
from   operator    import itemgetter

import h5py
import numpy as np
import pysam

from   glu.lib.illumina          import manifest_snps
from   glu.lib.seqlib.gc         import gc_window, cpg_window


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('manifest', help='Illumina BPM manifest file')

  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('-o', '--output', metavar='FILE', required=True,
                    help='Output model name')

  return parser


def create_gcmodel(filename, num_terms, num_snps):
  gcmodel      = h5py.File(filename, 'w')
  comp         = {}
  shuffle      = False
  vstr         = h5py.special_dtype(vlen=str)
  snptype      = [ ('name',vstr),('chromosome',vstr),('location',np.uint32) ]

  gcmodel.create_dataset('Terms', (num_terms,), vstr, maxshape=(None,),  **comp)
  gcmodel.create_dataset('GC',    (num_terms,num_snps), np.float32, fillvalue=np.nan, shuffle=shuffle,**comp)

  return gcmodel


def iter_pick_index(seq, indices, missing=None):
  indices = sorted(indices, reverse=True)

  if not indices:
    return

  last = indices[-1]

  for i,item in enumerate(seq):
    while i==last:
      yield item
      indices.pop()
      if not indices:
        return
      last = indices[-1]

  for i in indices:
    yield missing


def pick_gc_window(seq,winsize,indices):
  values = gc_window(seq,winsize)
  values = iter_pick_index(values,indices)
  return np.fromiter(values,count=len(indices),dtype=np.float32)


def pick_cpg_window(seq,winsize,indices):
  values = cpg_window(seq,winsize)
  values = list(iter_pick_index(values,indices))
  return np.array(values,dtype=np.float32)


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if 0:
    terms   = ['GC_1Mb',
               'CpG_1Mb',]
    windows = [50]
  elif 0:
    terms   = ['GC_1Mb', 'GC_500Kb', 'GC_250Kb', 'GC_100Kb', 'GC_50Kb', 'GC_25Kb', 'GC_10Kb', 'GC_5Kb', 'GC_1Kb', 'GC_500b', 'GC_100b', 'GC_50b',
               'CpG_1Mb','CpG_500Kb','CpG_250Kb','CpG_100Kb','CpG_50Kb','CpG_25Kb','CpG_10Kb','CpG_5Kb','CpG_1Kb','CpG_500b','CpG_100b','CpG_50b']
    windows = [1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000, 500, 100, 50]
  else:
    terms   = ['GC_1Mb', 'GC_250Kb', 'GC_100Kb', 'GC_50Kb', 'GC_10Kb', 'GC_1Kb', 'GC_500b', 'GC_100b', 'GC_50b',
               'CpG_1Mb','CpG_250Kb','CpG_100Kb','CpG_50Kb','CpG_10Kb','CpG_1Kb','CpG_500b','CpG_100b','CpG_50b']
    windows = [1000000, 250000, 100000, 50000, 10000, 1000, 500, 100, 50]

  reference = pysam.Fastafile(options.reference)

  sys.stderr.write('Reading Illumina manifest\n')
  allsnps = manifest_snps(options.manifest)
  allsnps = [ (row+(i,)) for i,row in enumerate(allsnps) ]
  allsnps.sort(key=itemgetter(1,2,0))

  m       = len(windows)
  n       = len(terms)
  s       = len(allsnps)

  assert n==m*2

  gcmodel               = create_gcmodel(options.output, n, s)
  gcmodel['Terms'][:]   = terms
  gcdata                = gcmodel['GC']

  attrs = gcmodel.attrs
  attrs['GLU_FORMAT']   = 'gcm'
  attrs['GLU_VERSION']  = 1
  attrs['ManifestPath'] = options.manifest
  attrs['ManifestName'] = os.path.basename(options.manifest)
  attrs['SNPCount']     = s
  attrs['TermCount']    = n

  for chrom,chrom_snps in groupby(allsnps, itemgetter(1)):
    try:
      seq = reference.fetch('chr'+chrom).upper()
    except IndexError:
      sys.stderr.write('Skipping region: %s\n' % chrom)
      continue

    seqlen = len(seq)
    sys.stderr.write('Processing region: chr%s, len=%d\n' % (chrom,seqlen))

    positions= []
    indices  = []

    skip     = 0
    for row in chrom_snps:
      pos    = row[2]-1
      index  = row[-1]
      assert 0<=index<s

      if 0<=pos<seqlen:
        positions.append(pos)
        indices.append(s-index-1)
      else:
        skip += 1

    if not indices:
      sys.stderr.write('  No valid SNPs found, skipping region: %s (%d SNPs)\n' % (chrom,skip))
      continue

    sys.stderr.write('  Generating GC/CpG windows for %d SNPs (skipped %d)\n' % (len(indices),skip))
    sys.stderr.write('  ... window size:')

    indices = np.array(indices, dtype=int)

    for i,winsize in enumerate(windows):
      sys.stderr.write(' %d' % winsize)

      win          = pick_cpg_window(seq,winsize,positions)

      gc           = gcdata[i]
      gc[indices]  = win[:,0]
      gcdata[i]    = gc

      cpg          = gcdata[i+m]
      cpg[indices] = win[:,1]
      gcdata[i+m]  = cpg

    filled = np.isfinite(gcdata[0]).sum()
    empty  = s-filled

    sys.stderr.write('\n')
    sys.stderr.write('  Wrote GC/CpG results (%d filled, %d empty)\n' % (filled,empty))

  sys.stderr.write('Closing GCM file\n')

  gcmodel.close()


if __name__=='__main__':
  main()
