# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert a GDAT file into series of GADA input files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import os
import sys

from   itertools                 import izip

import h5py
import numpy as np

from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib.transform import GenoTransform

from   glu.modules.cnv.gdat      import get_gcmodel, gc_correct



def decode(data, scale, nanval):
  nans       = data==nanval
  data       = data.astype(float)
  data[nans] = np.nan
  data      /= scale
  return data


def cnv_data(gdat,transform):
  include   = transform.samples.include
  exclude   = transform.samples.exclude
  rename    = transform.samples.rename

  snps      = gdat['SNPs'][:]
  samples   = gdat['Samples'][:]
  lrr       = gdat['LRR']
  baf       = gdat['BAF']

  lrr_scale = lrr.attrs['SCALE']
  lrr_nan   = lrr.attrs['NAN']

  baf_scale = baf.attrs['SCALE']
  baf_nan   = baf.attrs['NAN']

  for i,sample in enumerate(samples):
    if include is not None and sample not in include:
      continue

    if exclude is not None and sample in exclude:
      continue

    if rename is not None:
      sample = rename.get(sample,sample)

    if not sample:
      sys.stderr.write('Invalid null sample name... skipping.\n')
      continue

    sample_lrr = decode(lrr[i], lrr_scale, lrr_nan)
    sample_baf = decode(baf[i], baf_scale, baf_nan)

    yield sample,snps,sample_lrr,sample_baf


def split_fullname(filename):
  dirname  = os.path.dirname(filename)

  # Get filename
  filename = os.path.basename(filename)

  # Split filename into 1 or 2 parts up to the first '.'
  parts = filename.split('.',1)

  # Combine dirname and up to the first '.' of filename as prefix
  prefix = os.path.join(dirname,parts[0])

  sep = '_' if parts[0] else ''

  # Suffix the remainder of filename after the first '.'
  suffix = '' if len(parts) == 1 else parts[1]

  return prefix,sep,suffix


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdatfile', help='GLU GDAT input file')

  parser.add_argument('--gcmodel',    metavar='GCM', help='GC model file to use')
  parser.add_argument('--gcmodeldir', metavar='DIR', help='GC models directory')
  parser.add_argument('--includesamples', metavar='FILE', action='append',
                      help='List of samples to include, all others will be skipped')
  parser.add_argument('--excludesamples', metavar='FILE', action='append',
                    help='List of samples to exclude, only samples not present will be kept')
  parser.add_argument('--renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  parser.add_argument('-o', '--output', metavar='FILE', required=True,
                    help='Output genotype file name template')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  gccorrect = bool(options.gcmodel or options.gcmodeldir)
  prefix,sep,suffix = split_fullname(options.output)

  transform = GenoTransform.from_object(options)

  gdat      = h5py.File(options.gdatfile,'r')
  manifest  = gdat.attrs['ManifestName'].replace('.bpm','')
  data      = cnv_data(gdat,transform)

  if gccorrect:
    print 'Loading GC/CpG model for %s...' % manifest
    manifest  = gdat.attrs['ManifestName'].replace('.bpm','')
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    gcdesign,gcmask = get_gcmodel(filename)

  for sample,snps,lrr,baf in data:
    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chromosome','Position','Log R Ratio','B Allele Freq') )

    if gccorrect:
      lrr = gc_correct(lrr, gcdesign, gcmask)

    sample_data = izip(snps,lrr,baf)
    for (lname,chrom,location,alleles_forward),l,b in sample_data:
      out.writerow( (lname,chrom,location,l,b) )


if __name__ == '__main__':
  main()
