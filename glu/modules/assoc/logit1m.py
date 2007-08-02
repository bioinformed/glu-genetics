# -*- coding: utf-8 -*-
'''
File:          logit1m.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Perform a logistic genotype-phenotype association scan with interactions

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   numpy               import exp
from   scipy               import stats

from   glu.lib.glm         import GLogit
from   glu.lib.association import build_models,print_results,TREND,NULL


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  parser.add_option('-i', '--includesamples', dest='includesamples', metavar='FILE',
                    help='List of samples to include')
  parser.add_option('-d', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='List of samples to exclude')
  parser.add_option('-I', '--includeloci', dest='includeloci', metavar='FILE',
                    help='List of loci to include')
  parser.add_option('-D', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='List of loci to exclude')
  parser.add_option('-x', '--fixedloci', dest='fixedloci', metavar='FILE',
                    help='List of loci to include in every model')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')
  parser.add_option('--refalleles', dest='refalleles', metavar='FILE',
                    help='Mapping of locus name to the corresponding reference allele')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  loci,models = build_models(args[0], args[1], options)

  if options.nullmodel:
    null_model = models.build_model(NULL(),{})

    # Obtain estimates of covariate effects under the null
    null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
    null.fit()

    print_results(sys.stdout,null_model,null)

  headers = ['Locus', 'Score1', 'DF', 'p-value', 'OR',
                      'Score2', 'DF', 'p-value', 'OR',
                      'Score3', 'DF', 'p-value', 'OR1', 'OR2', 'OR3',
                      'Score4', 'DF', 'p-value']

  print '\t'.join(headers)

  # FIXME: Load fixed loci
  fixed_loci  = dict(fixed_loci)
  fixed_trend = sum( (TREND(l) for l in fixed_loci), NULL())
  fixed_null  = sum( ( NULL(l) for l in fixed_loci), NULL())

  # For each locus
  for locus in loci:
    lname,genos = locus

    if lname in fixed_loci:
      continue

    lmap = fixed_loci.copy()
    lmap[lname] = genos

    scan_model = TREND(lname)
    model1 = models.build_model(fixed_null +scan_model, lmap)
    model2 = models.build_model(fixed_trend+scan_model, lmap)
    model3 = models.build_model(fixed_trend+scan_model+fixed_trend*scan_model, lmap)

    if None in (model1,model2,model3):
      continue

    sys.stdout.write(lname)

    # Compute 1 df score test w/o 8q24 SNPs
    g = GLogit(model1.y,model1.X)
    g.fit()
    st,df = g.score_test(indices=[1]).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f' % exp(g.beta[1,0]))

    # Compute 1 df score test w/ 8q24 SNPs
    g = GLogit(model2.y,model2.X)
    g.fit()
    st,df = g.score_test(indices=[3]).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f' % exp(g.beta[3,0]))

    # Compute 2 df score test on just the interaction terms
    g = GLogit(model3.y,model3.X)
    g.fit()
    st,df = g.score_test(indices=(4,5)).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f\t%.3f\t%.3f' % (exp(g.beta[3,0]),exp(g.beta[4,0]),exp(g.beta[5,0])))

    # Compute 3 df score test on the main effect and interaction terms
    st,df = g.score_test(indices=(3,4,5)).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\n')


if __name__ == '__main__':
  main()
