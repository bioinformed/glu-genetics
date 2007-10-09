# -*- coding: utf-8 -*-
'''
File:          unscramble.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       Wed Sep 27 10:25:42 EDT 2006

Abstract:      Unscramble either locus identifier mapping between two
               genotype matricies by detecting least discordant matches
               between loci in one input with those in another, with
               automatic allele remapping.

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys
from   operator          import itemgetter
from   glu.lib.fileutils import autofile
from   glu.lib.remap     import remap_alleles, remap_category
from   glu.lib.genolib   import load_genostream,snp


def align_genotypes(genos1, genos2):
  sampleset = set(genos1.columns) && set(genos2.columns)
  samples   = [ c for c in cols1 if c in sampleset ]

  genos1 = genos1.sorted(sampleorder=samples)
  genos2 = genos2.sorted(sampleorder=samples)

  return genos1,genos2


def concordance(test,reference):
   assert test.columns == test.samples == reference.columns == reference.samples

   reference = list(reference)
   for testlocus,testgenos in test:
     results = []
     ident = 0
     for reflocus,refgenos in reference:
       genocounts = tally( (t,r) for t,r in izip(testgenos,refgenos) if t and r)
       concord,allelemap = remap_alleles(genocounts)
       result = testlocus,reflocus,concord,allelemap,remap_category(allelemap)
       results.append(result)
       if testlocus == reflocus:
         print
         print testgenos[:20]
         print refgenos[:20]
         ident = concord

     results.sort(key=lambda x: -x[2])
     result = results[0]
     if 1: # result[0] != result[1]:
       n = result[2]
       m = sum(genocounts.itervalues())
       print '%s\t%s\t%.3f\t%d\t%.3f\t%d\t%d' % (result[0],result[1],float(n)/m,n,float(ident)/m,ident,m)


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] test reference'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--testformat', dest='testformat',metavar='FILE', default='ldat',
                     help='File format for test genotype data. Values=ldat (default) or hapmap')
  parser.add_option('--refformat',  dest='refformat', metavar='FILE', default='ldat',
                     help='File format for reference genotype data. Values=ldat (default) or hapmap')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  test           = load_genostream(args[0], options.testformat).as_ldat()
  reference      = load_genostream(args[1], options.refformat).as_ldat()
  test,reference = align_genotypes(test,reference)
  concordance(test,reference)


if __name__=='__main__':
  main()
