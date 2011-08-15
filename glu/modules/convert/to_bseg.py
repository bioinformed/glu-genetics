# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert a GDAT file into series of BAFSegmentation input files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys

from   itertools                   import izip

import h5py

from   glu.lib.fileutils           import table_writer

from   glu.lib.genolib.transform   import GenoTransform

from   glu.modules.convert.to_gada import split_fullname, cnv_data
from   glu.modules.cnv.gdat        import get_gcmodel, gc_correct


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
  data      = cnv_data(gdat,transform)

  if gccorrect:
    print 'Loading GC/CpG model for %s...' % manifest
    manifest  = gdat.attrs['ManifestName'].replace('.bpm','')
    filename = options.gcmodel or '%s/%s.gcm' % (options.gcmodeldir,manifest)
    gcdesign,gcmask = get_gcmodel(filename)

  for sample,snps,geno,lrr,baf in data:
    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chr','Position','B Allele Frequency','Log R Ratio') )

    if gccorrect:
      lrr = gc_correct(lrr, gcdesign, gcmask)

    sample_data = izip(snps,lrr,baf)
    for (lname,chrom,location,alleles_forward),l,b in sample_data:
      out.writerow( (lname,chrom,location,b,l) )


if __name__ == '__main__':
  main()
