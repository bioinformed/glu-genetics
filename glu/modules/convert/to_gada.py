# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert a GDAT file into series of GADA input files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import os
import sys

from   itertools                 import izip

import numpy as np

from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib.transform import GenoTransform

from   glu.modules.cnv.gdat      import GDATFile, get_gcmodel, gc_correct


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
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('gdatfile', help='GLU GDAT input file')

  parser.add_argument('--signal',     metavar='TYPE', default='norm', choices=['norm','raw'],
                      help='LRR/BAF signal processing: raw or norm (default)')
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

  gdat      = GDATFile(options.gdatfile)
  data      = gdat.cnv_iter(transform, raw=options.signal.lower()=='raw')
  genomap   = {'AA':'AA','AB':'AB','BB':'BB','  ':'NC'}

  if gccorrect:
    manifest = gdat.attrs['ManifestName'].replace('.bpm','')
    print 'Loading GC/CpG model for %s...' % manifest
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    gcdesign,gcmask = get_gcmodel(filename, gdat.chromosome_index)

  for sample,snps,geno,lrr,baf in data:
    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chr','Position','Log.R.Ratio','B.Allele.Freq','GType') )

    if gccorrect:
      lrr = gc_correct(lrr, gcdesign, gcmask, minval=-2, maxval=2)

    sample_data = izip(snps,lrr,baf,geno)
    for (lname,chrom,location,alleles_forward),l,b,geno in sample_data:
      out.writerow( (lname,chrom,location,l,b,genomap[geno]) )


if __name__ == '__main__':
  main()
