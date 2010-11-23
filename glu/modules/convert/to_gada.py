# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert an Illumina Locus by DNA report file into a GLU genotype file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

from   itertools                 import izip

import h5py
import numpy as np

from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib.transform import GenoTransform


def decode(data, scale, nanval):
  nans       = data==nanval
  data       = data.astype(float)
  data[nans] = np.nan
  data      /= scale
  return data


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
  import optparse

  usage = 'usage: %prog [options] gdat'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--includesamples', dest='includesamples', metavar='FILE', action='append',
                    help='List of samples to include, all others will be skipped')
  parser.add_option('--excludesamples', dest='excludesamples', metavar='FILE', action='append',
                    help='List of samples to exclude, only samples not present will be kept')
  parser.add_option('--renamesamples', dest='renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output genotype file name')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  if not options.output:
    sys.stderr.write('ERROR: An output filename template must be specified')
    sys.exit(2)

  prefix,sep,suffix = split_fullname(options.output)

  transform = GenoTransform.from_object(options)

  gdat      = h5py.File(args[0],'r')

  snps      = gdat['SNPs'][:]
  samples   = gdat['Samples'][:]
  genos     = gdat['Genotype']
  lrr       = gdat['LRR']
  baf       = gdat['BAF']

  lrr_scale = lrr.attrs['SCALE']
  lrr_nan   = lrr.attrs['NAN']

  baf_scale = baf.attrs['SCALE']
  baf_nan   = baf.attrs['NAN']

  include   = transform.samples.include
  exclude   = transform.samples.exclude
  rename    = transform.samples.rename

  genomap   = {'AA':'AA','AB':'AB','BB':'BB','  ':'NC'}

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

    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chr','Position','Log.R.Ratio','B.Allele.Freq','GType') )

    sample_lrr = decode(lrr[i], lrr_scale, lrr_nan)
    sample_baf = decode(baf[i], baf_scale, baf_nan)

    sample_data = izip(snps,sample_lrr,sample_baf,genos[i])
    for (lname,chrom,location,alleles_forward),l,b,geno in sample_data:
      out.writerow( (lname,chrom,location,l,b,genomap[geno]) )


if __name__ == '__main__':
  main()
