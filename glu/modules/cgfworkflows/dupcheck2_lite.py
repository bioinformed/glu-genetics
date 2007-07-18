# -*- coding: utf-8 -*-
'''
File:          dupcheck2_lite.py

Authors:       Jun Lu          (lujun@mail.nih.gov)

Created:       2006-10-01

Abstract:      Detect expected and duplicate samples by genotype
               concordance (find_dups func was overwritten)

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys
import csv
import time

from glu.lib.genodata  import *
from scripts.dupcheck2 import *


def find_dups(samples,threshold=85,mincount=20,exp_dups=None,sample_phenos=None):
  threshold = float(threshold)/100
  samples = iter(samples)
  samples.next()
  samples = list(samples)
  sample_ids = set(map(itemgetter(0),samples))

  dup_pairs = []
  dup_sets  = union_find()

  dup_pairs.append( (None,None,None,None,None,1,1) )
  dup_pairs.append( (None,None,None,None,None,0,1) )
  dup_pairs.append( (None,None,None,None,None,1,0) )

  matchfunc = get_matchfunc(samples[0][1])

  if exp_dups is None:
    exp_dups = set()

  for (s1,genos1),(s2,genos2) in pair_generator(samples):
    s1,s2 = sorted([s1,s2])
    is_exp= (s1,s2) in exp_dups
    if not is_exp:
      continue
    matches,total = matchfunc(genos1, genos2)
    if len(genos1)< mincount or len(genos2) < mincount:
      raise ValueError, 'Error: the minimum counts (-m option) can not be greater that the number of loci'

    obs_dup = total and matches >= mincount and float(matches)/total >= threshold

    if not obs_dup and not is_exp:
      continue

    diff_pheno = 'NO|UNKNOWN'
    if sample_phenos is not None:
      s1_pheno,s2_pheno = sample_phenos.get(s1),sample_phenos.get(s2)
      if not s1_pheno or not s2_pheno:
        diff_pheno = 'UNKNOWN'
      elif s1_pheno != s2_pheno:
        diff_pheno = 'YES'
      else:
        diff_pheno = 'NO'

    if obs_dup:
      dup_sets.union(s1,s2)
      dup_pairs.append( (s1,s2,matches,total,diff_pheno,1, is_exp) )
    else:
      dup_pairs.append((s1,s2,matches,total,diff_pheno,0,1))

  return sample_ids,dup_pairs,dup_sets


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='A tab-delimited file with expected duplicates on each row, or with columns of sid,pid')
  parser.add_option('-d', '--dupout', dest='dupout', metavar='FILE',
                    help='Output of duplicate sets')
  parser.add_option('-f','--format',  dest='format', metavar='string', default='sdat',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat (default), trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', type='string', default='snp_acgt',
                    help='Input genotype representation.  Values=snp_acgt (default), snp_ab, snp_marker, or generic')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity shared between two individuals (default=85)')
  parser.add_option('-m', '--mincount', dest='mincount', metavar='N', type='int', default=20,
                    help='Minimum concordant genotypes to be considered informative for duplicate checking')
  parser.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows of data to N for testing purposes')
  #parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE',
  #                  help='Generate machine readable tabular output of results')

  parser.add_option('--phenofile', dest='phenofile', metavar='FILE',default=None,
                    help='A file containing sample phenotype information')
  parser.add_option('--phenos', dest='phenos', metavar='string',type='string',default=None,
                    help='Specify sample phenotypes (i.e. column headers in the phenotype file) to be included in output')
  parser.add_option('--pairformat', dest='pairformat', metavar='string',type='string',default='long',
                    help='Specify the output format (options: long or short)')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  print >> sys.stderr, 'Loading genotype data...'
  genorepr  = get_genorepr(options.genorepr)
  limit     = options.limit or None
  merger = VoteMerger()
  samples   = load_genostream(args[0], options.format, limit=limit, genorepr=genorepr, unique=False).as_sdat(merger)

  exp_dupsets = None
  if options.duplicates:
    print >> sys.stderr, 'Loading expected duplicates data...'
    exp_dupsets = load_expected_dupsets(options.duplicates)

  sample_phenos = None
  if options.phenos:
    if options.phenofile is not None:
      sample_phenos = load_sample_phenos(options.phenofile,options.phenos)
    else:
      print >> sys.stderr, 'Sample phenotype file is required'
      return

  print >> sys.stderr, 'Checking for duplicates...'
  sample_ids,dups,dup_sets = find_dups(samples,options.threshold,options.mincount,exp_dupsets,sample_phenos)

  if options.dupout:
    union_out = autofile(options.dupout, 'w')
    write_dupsets(options.dupout, dup_sets)

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  write_header(out,args[0],options.threshold,options.mincount)
  write_duplicate_sets(out,sample_ids,dup_sets,exp_dupsets,sample_phenos)
  write_dups_diff_phenos(out,dups,sample_phenos)
  write_duplicate_pairs(out,dups,sample_phenos,options.pairformat)


if __name__ == '__main__':
  main()
