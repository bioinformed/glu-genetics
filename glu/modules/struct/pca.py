# -*- coding: utf-8 -*-
'''
File:          pca.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2008-03-06

Abstract:      Performs Principle Components Analysis on genotype data

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   itertools                 import izip

import numpy

from   glu.lib.fileutils         import hyphen,autofile

from   glu.lib.genolib           import load_genostream
from   glu.lib.genolib.genoarray import count_genotypes, count_alleles_from_genocounts, minor_allele_from_allelecounts


def pca(data,center=True):
  data = numpy.asarray(data)

  if center:
    data = data-data.mean(axis=0)

  u,s,vt = numpy.linalg.svd(data, full_matrices = 0)
  return s**2/len(data),vt


def do_pca(genos,output,nvectors,center=True):
  genos = genos.as_ldat()

  data = []
  for (locus,geno),model in izip(genos,genos.models):
    genocounts   = count_genotypes(model,geno)
    allelecounts = count_alleles_from_genocounts(model,genocounts)

    try:
      minor,freq = minor_allele_from_allelecounts(model,allelecounts)
    except ValueError:
      minor = None

    if minor is None:
      continue

    other = [ a for a,n in izip(model.alleles[1:],allelecounts[1:]) if a!=minor and n ]

    if len(other) != 1:
      continue

    other = other[0]

    # Form numeric encoding of informative data (0,.5,1)
    gmap = { model[other,other]:0.0, model[other,minor]:0.5, model[minor,minor]:1.0 }
    inf  = [ gmap[g] for g in geno if g ]

    # Mean center and normalize scale
    s,n  = sum(inf),len(inf)
    avg  = s/n
    p    = (s+0.5)/(n+1)
    norm = numpy.sqrt(p*(1-p))

    # Form genotype to score mapping
    gmap = { model[None,None]  :  0.0,
             model[other,other]: (0.0-avg)/norm,
             model[other,minor]: (0.5-avg)/norm,
             model[minor,minor]: (1.0-avg)/norm }

    data.append( numpy.array([ gmap[g] for g in geno ], dtype=float) )

  data = numpy.asarray(data, dtype=float)

  values,vectors = pca(data,center=center)

  out = csv.writer(autofile(hyphen(output,sys.stdout),'w'),dialect='excel-tab')
  out.writerow(values[:nvectors])
  for i,sample in enumerate(genos.samples):
    out.writerow( [sample]+['%7.4f' % v for v in vectors[:nvectors,i]] )


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output transformed data to FILE (default is "-" for standard out)')
  parser.add_option('-f','--format',  dest='format',
                    help='Input format of genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')

  parser.add_option('-n', '--includesamples', dest='includesamples', metavar='FILE',
                    help='Include list for those samples to only use in the transformation and output')
  parser.add_option('-u', '--includeloci', dest='includeloci', metavar='FILE',
                    help='Include list for those loci to only use in the transformation and output')
  parser.add_option('-x', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='Exclude a list of samples from the transformation and output')
  parser.add_option('-e', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='Exclude a list of loci from the transformation and output')

  parser.add_option('--vectors', dest='vectors', metavar='N', type='int', default=10,
                    help='Output the top N eigenvectors.  Default=10')
  parser.add_option('--center', dest='center', metavar='YN', default='Y',
                    help='Perform second stage mean centering.  Values: Y or N.  Default=Y. '
                         'Set to N to match Eigenstrat.')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                  genome=options.loci)

  genos = genos.transformed(include_loci=options.includeloci,
                            exclude_loci=options.excludeloci,
                            include_samples=options.includesamples,
                            exclude_samples=options.excludesamples)

  center = options.center.upper()=='Y'
  do_pca(genos,options.output,nvectors=options.vectors,center=center)


if __name__ == '__main__':
  main()
