# -*- coding: utf-8 -*-
'''
File:          maf.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs various transformation on a genotype files

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   operator                  import itemgetter
from   itertools                 import izip

from   glu.lib.fileutils         import autofile,hyphen
from   glu.lib.genolib           import load_genostream
from   glu.lib.genolib.genoarray import count_alleles_from_genos


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f','--format',  dest='format', metavar='string',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of transformed data (default is "-" for standard out)')
  parser.add_option('-n', '--includesamples', dest='includesamples', metavar='FILE',
                    help='Include list for those samples to only use')
  parser.add_option('-u', '--includeloci', dest='includeloci', metavar='FILE',
                    help='Include list for those loci to only use')
  parser.add_option('-x', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='Exclude a list of samples')
  parser.add_option('-e', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='Exclude a list of loci')

  return parser


def compute_maf(model,genos):
  counts = count_alleles_from_genos(model,genos)
  n = sum(counts[1:])

  as   = model.alleles[1:]
  fs   = [ float(f)/n for f in counts[1:] ]
  afs  = sorted(izip(as,fs),key=itemgetter(1),reverse=True)
  miss = float(counts[0])/sum(counts)
  return miss,afs


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  loci = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                 genome=options.loci).as_ldat()

  loci = loci.transformed(include_loci=options.includeloci,
                          exclude_loci=options.excludeloci,
                          include_samples=options.includesamples,
                          exclude_samples=options.excludesamples)

  outfile = autofile(hyphen(options.output,sys.stdout),'w')
  out     = csv.writer(outfile,dialect='excel-tab')

  out.writerow(['Locus','%Missing','Major allele','frequency','...'])
  for (lname,genos),model in izip(loci,loci.models):
    missing,afs = compute_maf(model,genos)
    row = [lname, '%8.6f' % missing]
    for a,f in afs:
      row += [a,'%8.6f' % f]
    out.writerow(row)


if __name__ == '__main__':
  main()
