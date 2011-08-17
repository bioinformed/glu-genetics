# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert a GDAT file into series of PennCNV input files'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys

from   itertools                   import izip

from   glu.lib.fileutils           import table_writer

from   glu.lib.genolib.transform   import GenoTransform
from   glu.modules.convert.to_gada import split_fullname, option_parser
from   glu.modules.cnv.gdat        import GDATFile


def main():
  parser  = option_parser()
  options = parser.parse_args()

  prefix,sep,suffix = split_fullname(options.output)

  transform = GenoTransform.from_object(options)
  gdat      = GDATFile(options.gdatfile)
  data      = gdat.cnv_iter(transform)
  genomap   = {'AA':'AA','AB':'AB','BB':'BB','  ':'NC'}

  for sample,snps,geno,lrr,baf in data:
    out = table_writer('%s%s%s.%s' % (prefix,sep,sample,suffix))
    out.writerow( ('Name','Chr','Position','GType','Log R Ratio','B Allele Freq') )

    sample_data = izip(snps,geno,lrr,baf)
    for (lname,chrom,location,alleles_forward),geno,l,b in sample_data:
      out.writerow( (lname,chrom,location,genomap[geno],l,b) )


if __name__ == '__main__':
  main()
