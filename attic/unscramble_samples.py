# -*- coding: utf-8 -*-
'''
File:          unscramble.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       Wed Sep 27 10:25:42 EDT 2006

Abstract:      Unscramble either locus identifier mapping between two
               genotype matricies by detecting least discordant matches
               between samples in one input with those in another.

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys

from   operator                  import itemgetter

from   glu.lib.fileutils         import table_writer
from   glu.lib.genolib           import load_genostream
from   glu.lib.genolib.genoarray import genoarray_concordance


def align_genotypes(genos1, genos2):
  locusset = set(genos1.columns) & set(genos2.columns)
  loci     = [ c for c in genos1.columns if c in locusset ]

  genos1 = genos1.transformed(orderloci=loci,includeloci=locusset)
  genos2 = genos2.transformed(orderloci=loci,includeloci=locusset,
                              recode_models=genos1.genome)

  return genos1,genos2


def concordance(genos1,genos2):
   assert genos1.columns == genos1.loci == genos2.columns == genos2.loci

   genos2 = genos2.materialize()
   results = []
   for name1,sample1 in genos1:
     for name2,sample2 in genos2:
       matches,comparisons = genoarray_concordance(sample1,sample2)

       if not comparisons:
         continue

       conc = float(matches)/comparisons
       if float(matches)/comparisons > 0.90:
         results.append( [name1,name2,matches,comparisons,conc ] )

   results.sort(key=itemgetter(-1,-2,-3),reverse=True)

   out = table_writer(sys.stdout)
   out.writerow( ['SAMPLE1','SAMPLE2','MATCHES','COMPARISONS','CONCORDANCE'] )
   for name1,name2,matches,comparisons,conc in results:
     out.writerow( [name1,name2,str(matches),str(comparisons),'%f' % (conc*100,)] )


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] set1 set2'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--format1', dest='format1',metavar='FILE',
                     help='File format for genotype data.')
  parser.add_option('-g', '--genorepr1', dest='genorepr1', metavar='REP',
                    help='Genotype representation')
  parser.add_option('--format2',  dest='format2', metavar='FILE',
                     help='File format for genotype data.')
  parser.add_option('-G', '--genorepr2', dest='genorepr2', metavar='REP',
                    help='Genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  genos1 = load_genostream(args[0], format=options.format1,
                                    genorepr=options.genorepr1,
                                    genome=options.loci).as_sdat()
  genos2 = load_genostream(args[1], format=options.format2,
                                    genorepr=options.genorepr2,
                                    genome=options.loci).as_sdat()

  genos1,genos2 = align_genotypes(genos1,genos2)
  concordance(genos1,genos2)


if __name__=='__main__':
  main()
