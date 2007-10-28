# -*- coding: utf-8 -*-
'''
File:          logit1f.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Quick and dirty code to implement a trend LRT

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   scipy               import stats

from   glu.lib.glm         import GLogit
from   glu.lib.association import print_results,build_models,TREND,NULL,format_pvalue


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
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
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--allowdups', dest='allowdups', action='store_true', default=False,
                    help='Allow duplicate individuals in the data (e.g., to accommodate weighting '
                         'or incidence density sampling)')
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

  loci        = dict(loci)
  trend_terms = sum( (TREND(l) for l in loci), NULL())
  null_terms  = sum( ( NULL(l) for l in loci), NULL())

  null_model = models.build_model(NULL(),{})
  null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
  null.fit()
  print_results(sys.stdout,null_model,null)

  null_model = models.build_model(null_terms,loci)
  null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
  null.fit()
  L0 = null.L
  print_results(sys.stdout,null_model,null)

  alt_model = models.build_model(trend_terms,loci)
  alt = GLogit(alt_model.y,alt_model.X,vars=alt_model.vars)
  alt.fit()
  L1 = alt.L
  print_results(sys.stdout,alt_model,alt)

  lrt = -2*(L0-L1)
  df  = len(trend_terms)
  p   = format_pvalue(stats.distributions.chi2.sf(lrt,df))
  print '-2(lnL0-lnL1) ~ X2 = %.2f, df=%d, p=%s' % (lrt,df,p)


if __name__ == '__main__':
  main()
