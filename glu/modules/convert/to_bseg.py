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
from   glu.modules.convert.to_gada import split_fullname, option_parser, cnv_data


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if not options.output:
    sys.stderr.write('ERROR: An output filename template must be specified')
    sys.exit(2)

  prefix,sep,suffix = split_fullname(options.output)

  transform = GenoTransform.from_object(options)
  gdat      = h5py.File(options.gdatfile,'r')
  data      = cnv_data(gdat,transform)

  for sample,snps,geno,lrr,baf in data:
    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chr','Position','B Allele Frequency','Log R Ratio') )

    sample_data = izip(snps,lrr,baf)
    for (lname,chrom,location,alleles_forward),l,b in sample_data:
      out.writerow( (lname,chrom,location,b,l) )


if __name__ == '__main__':
  main()
