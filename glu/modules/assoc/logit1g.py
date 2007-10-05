# -*- coding: utf-8 -*-
'''
File:          logit1g.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Test a slew of logistic genotype-phenotype association models

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   numpy               import abs,bmat,exp
from   scipy               import stats
from   scipy.linalg        import LinAlgError

from   glu.lib.fileutils   import autofile,hyphen,load_list
from   glu.lib.genolib     import load_genostream
from   glu.lib.glm         import GLogit

from   glu.lib.association import print_results,build_models,format_pvalue, \
                                  NULL,GENO,ADOM,TREND,DOM,REC


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output summary results to FILE')
  parser.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output detailed results to FILE')
  parser.add_option('-p', '--detailsmaxp', dest='detailsmaxp', metavar='P', type='float', default=1.0,
                    help='Output detailed resutls for only pvalues below P threshold')
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
                    help='Minimum minor allele frequency filter (default=0.01)')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter (default=10)')
  parser.add_option('--refalleles', dest='refalleles', metavar='FILE',
                    help='Mapping of locus name to the corresponding reference allele')

  return parser


def eval_model(lname,k,model,details,detailsmaxp):
  assert 0<=k<=2

  res = ['']*([1,9,12][k])

  if model is None:
    return res

  try:
    g = GLogit(model.y,model.X,vars=model.vars)
    L,b,W = g.fit()
    ps = pw = format_pvalue(1)
    if k:
      st,df = g.score_test(indices=range(1,k+1)).test()
      wt,df = g.wald_test(indices=range(1,k+1)).test()
      sf = stats.distributions.chi2.sf
      ps = format_pvalue(sf(st,df))
      pw = format_pvalue(sf(wt,df))
  except LinAlgError:
    return res

  if k == 2:
    res = ['%.6f'  % L,
           '%8.5f' % st, '%d' % df, ps,
           '%8.5f' % wt, '%d' % df, pw,
           '%7.4f' % exp(b[1,0]), '%7.4f' % exp(b[2,0]),
            '%.6f' % W[1,1], '%.6f' % W[2,2], '%.6f' % W[1,2]]
  elif k == 1:
    res = ['%.6f'  % L,
           '%8.5f' % st, '%d' % df, ps,
           '%8.5f' % wt, '%d' % df, pw,
           '%7.4f' % exp(b[1,0]), '%.6f' % W[1,1]]
  elif k == 0:
    res = ['%.6f' % L]

  if k and details and min(ps,pw) <= detailsmaxp:
    details.write('\nRESULTS: %s\n\n' % lname)
    print_results(details,model,g)
    details.write('Score Test: X2=%6.3f, df=%d, p=%12.10f\n' % (st,df,ps))
    details.write('Wald  Test: X2=%6.3f, df=%d, p=%12.10f\n' % (wt,df,pw))
    details.write('-'*79)
    details.write('\n')

  return res


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout) or sys.stdout,'w')
  details = None
  if options.details:
    details = autofile(hyphen(options.details,sys.stdout),'w')
    if details is out:
      raise ValueError('Cannot send summary and detailed output to stdout')

  loci,models = build_models(args[0], args[1], options)

  if options.nullmodel:
    null_model = models.build_model(NULL(), {})

    # Obtain estimates of covariate effects under the null
    null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
    L_null,b_null,W_null = null.fit()

    print_results(out,null_model,null)
    if options.details:
      print_results(details,null_model,null)

  headers = ['Locus', 'lnL_null',
             'lnL_genotype',
             'Genotype Score',  'DF', 'p-value',
             'Genotype Wald',   'DF', 'p-value',
             'OR_het', 'OR_hom', 'var_het', 'var_hom', 'cov_homhet',
             'lnL_adddom',
             'AddDom Score',    'DF', 'p-value',
             'AddDom Wald',     'DF', 'p-value',
             'OR_add', 'OR_dom', 'var_add', 'var_dom', 'cov_adddom',
             'lnL_trend',
             'Trend Score',     'DF', 'p-value',
             'Trend Wald',      'DF', 'p-value',
             'OR_trend',        'var_trend',
             'lnL_recessive',
             'Recessive Score', 'DF', 'p-value',
             'Recessive Wald',  'DF', 'p-value',
             'OR_recessive',    'var_recessive',
             'lnL_dominant',
             'Dominant Score',  'DF', 'p-value',
             'Dominant Wald',   'DF', 'p-value',
             'OR_dominant',   ' var_dominant',]

  out.write('\t'.join(headers))
  out.write('\n')

  # For each locus
  for locus in loci:
    lname = locus[0]
    l = dict([locus])

    model_null  = models.build_model( NULL(lname), l)
    model_geno  = models.build_model( GENO(lname), l)
    model_adom  = models.build_model( ADOM(lname), l)
    model_trend = models.build_model(TREND(lname), l)
    model_rec   = models.build_model(  REC(lname), l)
    model_dom   = models.build_model(  DOM(lname), l)

    all_models  = [ model_null, model_geno, model_trend, model_dom, model_rec ]

    if all(m is None for m in all_models):
      if details:
        details.write('Skipping locus due to invalid model: %s\n\n' % lname)
      continue

    results  = [lname]
    results += eval_model(lname, 0, model_null,  details, options.detailsmaxp)
    results += eval_model(lname, 2, model_geno,  details, options.detailsmaxp)
    results += eval_model(lname, 2, model_adom,  details, options.detailsmaxp)
    results += eval_model(lname, 1, model_trend, details, options.detailsmaxp)
    results += eval_model(lname, 1, model_rec,   details, options.detailsmaxp)
    results += eval_model(lname, 1, model_dom,   details, options.detailsmaxp)

    out.write('\t'.join(results))
    out.write('\n')


if __name__ == '__main__':
  main()
