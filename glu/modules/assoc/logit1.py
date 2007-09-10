# -*- coding: utf-8 -*-
'''
File:          logit1.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Fit generalized logistic phenotype-genotype association models

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   numpy               import empty,exp,nan,abs,isfinite,vstack,vsplit
from   scipy               import stats

from   glu.lib.fileutils   import autofile,hyphen
from   glu.lib.glm         import GLogit

from   glu.lib.association import build_models,print_results,get_term,format_pvalue,NULL


#####
##### FIXME: This only works for GENO models.  Other models need to have custom OR columns
def get_ors(g,cats,k,indices):
  ors = exp(g.beta).take(indices).A.flatten()
  orsl = empty( (2*(cats-1),), dtype=float )
  orsl[:] = nan
  orsl[::3-k] = ors
  return orsl


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
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--genomodel', dest='genomodel', default='geno,trend', metavar='M1,M2,..',
                    help='Comma separated list of genetic models.  The first that can be fit will be used.  '
                         'Values: genotype/geno, adddom/adom, trend/multiplicative/mult, additive/add, '
                         'dominant/dom, recessive/rec, missing/miss, not_missing/not_miss, null.  Default=geno,trend')
  parser.add_option('--refalleles', dest='refalleles', metavar='FILE',
                    help='Mapping of locus name to the corresponding reference allele')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout) or sys.stdout,'w')
  if options.details:
    details = autofile(hyphen(options.details,sys.stdout),'w')
    if details is out:
      raise ValueError('Cannot send summary and detailed output to stdout')

  loci,models = build_models(args[0], args[1], options)

  null_model = models.build_model(NULL(),{})

  # Obtain estimates of covariate effects under the null
  null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
  cats = len(null.categories)
  null.fit()

  if options.nullmodel:
    print_results(out,null_model,null)
    if options.details:
      details.write('NULL MODEL: %s\n\n')
      print_results(details,null_model,null)

  headers = ['Locus']

  or_headers = []
  if len(null.categories) > 2:
    for i in range(1,len(null.categories)):
      or_headers += ['HetOR%d' % i, 'HomOR%d' % i]
  else:
    or_headers += ['HetOR', 'HomOR']

  headers += ['Score X2', 'p-value', 'Wald X2', 'p-value', 'LR X2', 'p-value', 'df']
  headers += or_headers

  out.write('\t'.join(headers))
  out.write('\n')

  initial_beta = [None]
  bs = vsplit(null.beta,cats-1)
  for i in range(1,3):
    bi = []
    for b in bs:
      bi.extend( [b[:1],[[0]]*i,b[1:]] )
    initial_beta.append(vstack(bi))

  terms = [ get_term(m) for m in options.genomodel.split(',') ]

  if not terms:
    raise ValueError('Must specify a genetic model to test')

  # For each locus
  for locus in loci:
    lname = locus[0]
    lmap  = dict([locus])

    # Scan through each model and find the first that can be built
    for term in terms:
      model_term = term(lname)
      model = models.build_model(model_term,lmap)
      if model:
        break
    else:
      # Otherwise, skip the locus
      continue

    if 0:
      f = csv.writer(file('%s.csv' % lname,'w'))
      f.writerow(model.vars)
      f.writerows(model.X.tolist())

    n = model.X.shape[1]
    k = len(model_term)

    if not k or not model.X.shape[0]:
      continue
    elif k>2:
      raise ValueError,'Unexpected number of parameters in model (n=%d,k=%d)' % (n,k)

    out.write(lname)

    ### SCORE TEST ###
    g = GLogit(model.y,model.X)

    # Construct genotype parameter indices
    indices = [ j*n+i for i in range(1,k+1)
                      for j in range(len(g.categories)-1) ]

    # FIXME: Can use initial_betas to possibly speed things up
    #g.fit(initial_beta=initial_beta[k])
    #scoretest = g.score_test(indices=indices,initial_beta=null.beta)

    g.fit()

    st,df_s = g.score_test(indices=indices).test()
    wt,df_w = g.wald_test(indices=indices).test()
    lt,df_l = g.lr_test(indices=indices).test()

    sp = format_pvalue(stats.distributions.chi2.sf(st,df_s))
    wp = format_pvalue(stats.distributions.chi2.sf(wt,df_w))
    lp = format_pvalue(stats.distributions.chi2.sf(lt,df_l))

    out.write('\t%.5f\t%s' % (st,sp))
    out.write('\t%.5f\t%s' % (wt,wp))
    out.write('\t%.5f\t%s' % (lt,lp))

    assert df_s == df_w == df_l
    out.write('\t%d' % df_s)

    ors = exp(g.beta).take(indices).A.flatten()

    ors = get_ors(g,cats,k,indices)
    orsstr = '\t'.join('%.4f' % orr if isfinite(orr) else '' for orr in ors)
    out.write('\t%s' % orsstr)

    out.write('\n')

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results(details,model,g)
      details.write('Score test           : X2=%9.5f, df=%d, p=%s\n' % (st,df_s,sp))
      details.write('Wald test            : X2=%9.5f, df=%d, p=%s\n' % (wt,df_w,wp))
      details.write('Likelihood ratio test: X2=%9.5f, df=%d, p=%s\n' % (lt,df_l,lp))
      details.write('\n')
      details.write('-'*79)
      details.write('\n')


if __name__ == '__main__':
  main()
