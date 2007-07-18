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
from   itertools         import islice, chain
from   glu.lib.utils     import autofile, peekfirst
from   glu.lib.genoarray import snp_acgt
from   glu.lib.remap     import remap_alleles, remap_category
from   glu.lib.genodata  import load_genomatrixstream, reorder_genomatrix_columns


def align_genotypes(genos1, genos2):
  genos1 = list(genos1)
  genos2 = list(genos2)

  cols1,genos1 = peekfirst(genos1)
  cols2,genos2 = peekfirst(genos2)

  genos1 = list(genos1)
  genos2 = list(genos2)

  colset = set(cols1) & set(cols2)
  cols   = [ c for c in cols1 if c in colset ]

  genos1 = reorder_genomatrix_columns(genos1, cols)
  genos2 = reorder_genomatrix_columns(genos2, cols)

  genos1 = list(genos1)
  genos2 = list(genos2)

  return genos1,genos2


def concordance(test,reference):
   test        = iter(test)
   testsamples = test.next()

   reference   = iter(reference)
   refsamples  = reference.next()

   assert testsamples == refsamples

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
         print list(snp_acgt.genos_str(testgenos[:20]))
         print list(snp_acgt.genos_str(refgenos[:20]))
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

  test           = load_genomatrixstream(args[0], format=options.testformat)
  reference      = load_genomatrixstream(args[1], format=options.refformat)
  reference      = list(reference)
  test,reference = align_genotypes(test,reference)

  concordance(test,reference)


if __name__=='__main__':
  main()
