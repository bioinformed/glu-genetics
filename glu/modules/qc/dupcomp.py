# -*- coding: utf-8 -*-
'''
File:          dupcomp.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       June 13, 2006

Abstract:      Reports the discordance between the supplied duplicates for
               each locus.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

## NOTE ## THIS MODULE IS SLATED TO BE DEPRECATED

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   textwrap              import fill

from   glu.lib.utils         import percent, pair_generator
from   glu.lib.fileutils     import autofile
from   glu.lib.union_find    import union_find


def duplicate_output_concordance(data, details=False):
  data = sorted( ((percent(c,n),l,c,n,d) for l,c,n,d in data))
  print
  print '   Locus            Conc / Total   %Conc       Discordant Individuals'
  print '   ------------     ------------   -------     ---------------------------'
  e = ' '*22
  for r,(p,l,c,n,diff) in enumerate(data):
    txt=''
    if diff and details:
      txt = ','.join(diff[0])
    print '   %-12s   %4d/%4d       %5.1f%%     %s' % (l,c,n,p,txt)
    if diff and details:
      for d in diff[1:]:
        txt = ','.join(d)
        print '   %67s' % txt
  print


def locus_concordance(dupsets, genos, locus_ids):
  for i,locus in enumerate(locus_ids):
    conc=comp=0
    diff = []

    for dset in dupsets:
      g = [ (ind,genos[ind][i+1]) for ind in dset if genos[ind][i+1] ]

      for (i1,g1),(i2,g2) in pair_generator(g):
        if g1==g2:
          conc+=1
        else:
          diff.append( (i1,i2) )

        comp+=1

    yield locus,conc,comp,diff


def load_samples(filename, format):
  samples = csv.reader(autofile(filename),dialect=format)
  for row in samples:
    yield tuple(intern(r) for r in row)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] dupfile genofile'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', default='excel',
                    help='Output of duplicate sets')
  return parser


def main():
  import sys,time

  parser = option_parser()
  options,args = parser.parse_args()

  if not len(args) == 2:
    parser.print_help()
    return

  print 'Duplicate Locus Concordance'
  print '  Duplicates:    ', args[0]
  print '  Genotypes:     ', args[1]
  print '  Timestamp:     ', time.asctime()

  expected_dupset = union_find()

  print >> sys.stderr, 'Loading duplicates...',
  for dupset in csv.reader(autofile(args[0]), options.format):
    for dup in dupset:
      expected_dupset.union(dupset[0],dup)
  print >> sys.stderr, 'Done.'

  print >> sys.stderr, 'Loading genotypes...',
  samples = load_samples(args[1], options.format)
  locus_ids = samples.next()[1:]

  genos={}
  for s in samples:
    genos[s[0]]=s[1:]

  dupsets = expected_dupset.sets()

  print >> sys.stderr, 'Done.'

  loc_conc = locus_concordance(dupsets, genos, locus_ids)
  loc_conc = list(loc_conc)

  print 'LOCUS CONCORDANCE FOR EXPECTED DUPLICATE INDIVIDUALS SUMMARY'
  duplicate_output_concordance(loc_conc)
  print 'LOCUS CONCORDANCE FOR EXPECTED DUPLICATE INDIVIDUALS'
  duplicate_output_concordance(loc_conc, details=True)


if __name__ == '__main__':
  main()
