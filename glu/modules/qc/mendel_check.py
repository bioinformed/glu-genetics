# -*- coding: utf-8 -*-
'''
File:          mendel_check.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Abstract:      Detects non-Mendialian inheritance among parent offspring trios

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import csv

from operator  import itemgetter
from itertools import repeat


def percent(a,b):
  if not b:
    return 0.
  return float(a)/b*100


def genorepr(g):
  if not g:
    return None
  elif '/' in g:
    a1,a2 = g.split('/')
  elif len(g) == 2:
    a1,a2 = g
  else:
    raise ValueError, 'Invalid genotype: %s' % g

  if a1>a2:
    a1,a2 = a2,a1

  return intern(a1),intern(a2)


def load_genomatrix(genofile):
  rows = csv.reader(genofile,dialect='excel-tab')

  header = rows.next()[1:]

  yield header

  for row in rows:
    label = row[0]
    genos = map(genorepr,row[1:])
    yield label,genos


def parent_offspring_concordance(parent1, parent2, child, locusstats):
  # Must have informative child and at least one parent
  if not child or not (parent1 and parent2):
    return 0,0

  if not parent1:
    parent1 = repeat(None)
  if not parent2:
    parent2 = repeat(None)

  i = n = 0
  for locusstat,p1,p2,c in zip(locusstats,parent1,parent2,child):
    # Must have informative child and at least one parent genotype
    if not c or not (p1 and p2):
      continue

    locusstat[1] += 1

    # If ensure that p1 is not missing
    if not p2 and p1:
      p1,p2 = p2,p1

    # Check Parent1 -> Offspring case
    if p1 and not p2 and c[0] in p1 or c[1] in p1:
      i += 1
      locusstat[0] += 1

    # CHeck Parent1,Parent2 -> Offspring case
    elif (c[0] in p1 and c[1] in p2) or (c[1] in p1 and c[0] in p2):
      i += 1
      locusstat[0] += 1

    n += 1

  return i,n



def load_pedigree(pedfilename):
  pedfile = csv.reader(file(pedfilename),dialect='excel-tab')
  return [ row[:3] for row in pedfile if row[0] ]


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-p', '--pedigree', dest='pedigree', metavar='FILE',
                    help='A tab delimited pedigree file ')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output Mendelian concordance by sample')
  parser.add_option('-l', '--locout', dest='locout', metavar='FILE',
                    help='Output Mendelian concordance by locus')
  return parser


def main():
  import sys

  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  if not options.pedigree:
    print >> sys.stderr, 'A pedigree file must be specified (--pedigree)'
    return

  relationships = load_pedigree(options.pedigree)

  if options.output == '-':
    sampleout = sys.stdout
  else:
    sampleout = file(options.output, 'w')

  sampleout = csv.writer(sampleout,dialect='excel-tab')

  samples = load_genomatrix(file(args[0]))

  # Get locus names
  loci = samples.next()

  # Index samples
  samples = dict(samples)

  # Initialize statistics
  samplestats = []
  locusstats  = [ [0,0] for i in range(len(loci)) ]

  # Check all parent child relationships
  for child,parent1,parent2 in relationships:
    i,n = parent_offspring_concordance(samples.get(parent1), samples.get(parent2),
                                       samples.get(child), locusstats)
    if i!=n:
      samplestats.append( (child,parent1,parent2,i,n,percent(i,n)) )

  # Build and sort resulting statistics
  locusstats = sorted( (locus,i,n,percent(i,n)) for locus,(i,n) in zip(loci,locusstats) )
  samplestats.sort(key=itemgetter(5))

  # Produce output
  sampleout.writerow( ['CHILD','PARENT1','PARENT2','CONCORDANT','TOTAL','RATE'] )
  sampleout.writerows(samplestats)

  if options.locout:
    locout = csv.writer(file(options.locout,'w'),dialect='excel-tab')
    locout.writerow( ['LOCUS','CONCORDANT','TOTAL','RATE'] )
    locout.writerows(locusstats)


if __name__ == '__main__':
  main()
