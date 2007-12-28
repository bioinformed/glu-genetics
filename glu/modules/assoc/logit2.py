# -*- coding: utf-8 -*-
'''
File:          logit2.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Two-locus logistic genotype-phenotype association scan

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   itertools           import chain

from   scipy               import stats
from   scipy.linalg        import LinAlgError

from   glu.lib.fileutils   import autofile,hyphen,load_list
from   glu.lib.glm         import GLogit

from   glu.lib.association import print_results,build_models,format_pvalue,NULL,GENO,TREND


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of summary report')
  parser.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output of summary report')
  parser.add_option('-v', '--verbose', dest='verbose', metavar='LEVEL', type='int', default=1,
                    help='Verbosity level of diagnostic output.  O for none, 1 for some (default), 2 for exhaustive.')
  parser.add_option('-i', '--includesamples', dest='includesamples', metavar='FILE',
                    help='List of samples to include')
  parser.add_option('-d', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='List of samples to exclude')
  parser.add_option('-I', '--includeloci', dest='includeloci', metavar='FILE',
                    help='List of loci to include')
  parser.add_option('-D', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='List of loci to exclude')
  parser.add_option('--subset', dest='subset', metavar='FILE',
                    help='Only consider subset of 2-locus models around the list of SNPs in FILE')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--allowdups', dest='allowdups', action='store_true', default=False,
                    help='Allow duplicate individuals in the data (e.g., to accommodate weighting '
                         'or incidence density sampling)')
  parser.add_option('-w', '--window', dest='window', metavar='N', type='int', default=25,
                    help='Window size around each SNP')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')
  parser.add_option('--refalleles', dest='refalleles', metavar='FILE',
                    help='Mapping of locus name to the corresponding reference allele')
  return parser


def window(loci,n):
  '''
  >>> list(window(range(5),3))
  [([], 0, [1, 2, 3]), ([0], 1, [2, 3, 4]), ([], 2, [3, 4]), ([0, 1, 2], 3, [4]), ([1, 2, 3], 4, [])]
  '''
  import collections

  items = collections.deque()
  loci = iter(loci)

  # Fill queue to at most n+1 items
  for locus in loci:
    items.append(locus)
    if len(items) == n+1:
      break

  # Return first partition
  l = list(items)
  yield [],l[0],l[1:]

  # Return middle windows
  i = 0
  for locus in loci:
    if i==n:
      items.popleft()
    else:
      i+=1

    items.append(locus)
    l = list(items)
    yield l[:i],l[i],l[i+1:]

  # Return final partitions
  l = list(items)
  for i in xrange(i,len(l)-1):
    i += 1
    yield l[i-n:i],l[i],l[i+1:]


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout) or sys.stdout,'w')

  loci,models = build_models(args[0], args[1], options)

  subset = None
  if options.subset:
    subset = set(load_list(options.subset))

  if options.nullmodel:
    null_model = models.build_model(NULL(), {})

    # Obtain estimates of covariate effects under the null
    null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
    null.fit()

    print_results(out,null_model,null)

  headers = ['Ref_Locus',  'Alt_Locus',
             'Ref_Score',  'DF', 'p-value',
             'Alt_Score',  'DF', 'p-value',
             'Cond_Score', 'DF', 'p-value',
             'Joint_Score','DF', 'p-value']

  out.write('\t'.join(headers))
  out.write('\n')

  # For each locus
  for left,locus,right in window(loci,options.window):
    lname1 = locus[0]

    if subset is not None and lname1 not in subset:
      continue

    for mod1 in GENO,TREND:
      model1_term = mod1(lname1)
      model1 = models.build_model(model1_term, dict([locus]))
      if model1:
        break
    else:
      continue

    g  = GLogit(model1.y,model1.X,vars=model1.vars)
    n  = len(model1.vars)
    k1 = len(model1_term)

    indices = [ j*n+i for i in range(1,1+k1)
                      for j in range(len(g.categories)-1) ]

    try:
      st1,df1 = g.score_test(indices=indices).test()
    except LinAlgError:
      continue

    for other in chain(left,right):
      lmap = dict([locus,other])
      lname2 = other[0]

      for mod2 in GENO,TREND:
        model2_term = mod2(lname2)+NULL(lname1)
        model2 = models.build_model(model2_term, lmap)
        if model2:
          break
      else:
        continue

      g  = GLogit(model2.y,model2.X,vars=model2.vars)
      n  = len(model2.vars)
      k2 = len(model2_term)
      indices = [ j*n+i for i in range(1,1+k2)
                        for j in range(len(g.categories)-1) ]

      try:
        st2,df2 = g.score_test(indices=indices).test()
      except LinAlgError:
        continue

      interaction_term = TREND(lname1)*TREND(lname2)
      term = model1_term+model2_term+interaction_term
      model = models.build_model(term, lmap)

      if not model:
        interaction_term = NULL(lname1)+NULL(lname2)
        term  = model1_term+model2_term+interaction_term
        model = models.build_model(model1_term+model2_term, lmap)

      if not model:
        continue

      g = GLogit(model.y,model.X,vars=model.vars)
      n = len(model.vars)
      k = len(term)
      indices = [ j*n+i for i in range(1+k1,1+k)
                        for j in range(len(g.categories)-1) ]

      try:
        st,df = g.score_test(indices=indices).test()
      except LinAlgError:
        continue

      indices = [ j*n+i for i in range(1,1+k)
                        for j in range(len(g.categories)-1) ]

      try:
        stj,dfj = g.score_test(indices=indices).test()
      except LinAlgError:
        stj,dfj = 0,k

      sf = stats.distributions.chi2.sf
      out.write('\t'.join([lname1,lname2,
                       '%8.5f' % st1, '%d' % df1, format_pvalue(sf(st1,df1)),
                       '%8.5f' % st2, '%d' % df2, format_pvalue(sf(st2,df2)),
                       '%8.5f' % st,  '%d' % df,  format_pvalue(sf(st, df )),
                       '%8.5f' % stj, '%d' % dfj, format_pvalue(sf(stj,dfj))]))
      out.write('\n')


if __name__ == '__main__':
  main()
