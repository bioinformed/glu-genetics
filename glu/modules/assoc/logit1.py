# -*- coding: utf-8 -*-
'''
File:          logit1.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Fit generalized logistic phenotype-genotype association models

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   numpy               import isfinite
from   scipy               import stats

from   glu.lib.fileutils   import autofile,hyphen,table_writer
from   glu.lib.glm         import GLogit,LinAlgError

from   glu.lib.association import build_models,print_results,get_term,format_pvalue,NULL


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output summary results to FILE')
  parser.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output detailed results to FILE')
  parser.add_option('-v', '--verbose', dest='verbose', metavar='LEVEL', type='int', default=1,
                    help='Verbosity level of diagnostic output.  O for none, 1 for some (default), 2 for exhaustive.')
  parser.add_option('-p', '--detailsmaxp', dest='detailsmaxp', metavar='P', type='float', default=1.0,
                    help='Output detailed results for only pvalues below P threshold')
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
  parser.add_option('--allowdups', dest='allowdups', action='store_true', default=False,
                    help='Allow duplicate individuals in the data (e.g., to accommodate weighting '
                         'or incidence density sampling)')
  parser.add_option('--genomodel', dest='genomodel', default='geno,trend', metavar='M1,M2,..',
                    help='Comma separated list of genetic models.  The first that can be fit will be used.  '
                         'Values: genotype/geno, adddom/adom, trend/multiplicative/mult, additive/add, '
                         'dominant/dom, recessive/rec, missing/miss, not_missing/not_miss, null.  Default=geno,trend')
  parser.add_option('--tests', dest='tests', default='score', metavar='T1,T2,..',
                    help='Comma separated list of tests to apply to each model.  Supported tests '
                         'include score, Wald, and likelihood ratio tests.  Values: score, wald, lrt.')
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

  if not null_model:
    raise ValueError('Cannot construct null model')

  # Obtain estimates of covariate effects under the null
  null = GLogit(null_model.y,null_model.X,vars=null_model.vars)

  null.fit()

  if options.nullmodel:
    print_results(out,null_model,null)
    if options.details:
      details.write('NULL MODEL:\n\n')
      print_results(details,null_model,null)

  terms = [ get_term(m) for m in options.genomodel.split(',') ]

  if not terms:
    raise ValueError('Must specify a genetic model to test')

  tests = set(t.strip().lower() for t in options.tests.split(','))
  tests.discard('')
  extra = tests - set(['score','wald','lrt'])
  if extra:
    raise ValueError('Unknown test(s) specified: %s' % ','.join(sorted(extra)))

  headers = ['Locus', 'Alleles', 'MAF', 'Geno Counts', 'Subjects']

  or_headers = []
  if len(null.categories) > 2:
    for i in range(1,len(null.categories)):
      or_headers += ['HetOR%d' % i, 'HetOR%d 95%% CI_l' % i, 'HetOR%d 95%% CI_u' % i,
                     'HomOR%d' % i, 'HomOR%d 95%% CI_l' % i, 'HomOR%d 95%% CI_u' % i ]
  else:
    or_headers += ['HetOR', 'HetOR 95% CI_l', 'HetOR 95% CI_u',
                   'HomOR', 'HomOR 95% CI_l', 'HomOR 95% CI_u']

  if 'score' in tests:
    headers += ['score X2', 'score p-value']
  if 'wald' in tests:
    headers += ['Wald X2',  'Wald p-value' ]
  if 'lrt' in tests:
    headers += ['LR X2',    'LR p-value'   ]
  if tests:
    headers += ['df']

  headers += or_headers

  out.write('\t'.join(headers))
  out.write('\n')

  # For each locus
  for locus in loci:
    lname = locus[0]
    lmap  = dict([locus])

    # Scan through each model and find the first that can be built
    for term in terms:
      model_term = term(lname)
      model = models.build_model(model_term,lmap)

      if not model:
        continue

      g = GLogit(model.y,model.X,vars=model.vars)

      try:
        g.fit()
      except LinAlgError:
        continue
      else:
        break
    else:
      # Otherwise, skip the locus
      continue

    # Desgin matrix debugging output
    if 0:
      f = table_writer('%s.csv' % lname,dialect='csv')
      f.writerow(model.vars)
      f.writerows(model.X.tolist())

    # R model verfication debugging code
    if 0:
      import rpy
      from   rpy   import r
      from   numpy import array,allclose

      vars = [ v.replace(':','.').replace('+','p').replace('-','m').replace('_','.') for v in model.vars[1:] ]
      frame = dict( (v,model.X[:,i+1].A.reshape(-1)) for i,v in enumerate(vars) )
      frame['y'] = model.y.A.reshape(-1)
      formula = 'y ~ ' + ' + '.join(v.replace(':','.') for v in vars)

      rpy.set_default_mode(rpy.NO_CONVERSION)
      mod = r.glm(r(formula),data=r.data_frame(**frame),family=r.binomial('logit'))
      rpy.set_default_mode(rpy.BASIC_CONVERSION)
      pmod = mod.as_py()

      coef = r.coefficients(mod)
      coef=array([coef['(Intercept)']] + [ coef[v] for v in vars ],dtype=float)
      coef2 = g.beta.A.reshape(-1)

      #assert allclose(coef,g.beta.A.reshape(-1),atol=1e-6)
      #assert allclose(r.vcov(mod),g.W,atol=1e-6)

    assert len(null.categories) >= len(g.categories)

    n = model.X.shape[1]
    k = len(model_term)
    c = len(g.categories)-1

    # Construct genotype parameter indices
    indices = [ j*n+i for j in range(c)
                      for i in model_term.indices() ]

    m = model.model_loci[lname]
    out.write('\t'.join([lname, ','.join(m.alleles),
                                '%.3f' % m.maf,
                                '|'.join(map(str,m.counts))]))
    out.write('\t')

    counts = [ str((g.y==cat).sum()) for cat in g.categories]
    out.write('|'.join(counts))

    sp = wp = lp = 1

    if 'score' in tests:
      st,df = g.score_test(indices=indices).test()
      sp    = stats.distributions.chi2.sf(st,df)
      sps   = format_pvalue(sp)
      out.write('\t%.5f\t%s' % (st,sps))

    if 'wald' in tests:
      wt,df = g.wald_test(indices=indices).test()
      wp    = stats.distributions.chi2.sf(wt,df)
      wps   = format_pvalue(wp)
      out.write('\t%.5f\t%s' % (wt,wps))

    if 'lrt' in tests:
      lt,df = g.lr_test(indices=indices).test()
      lp    = stats.distributions.chi2.sf(lt,df)
      lps   = format_pvalue(lp)
      out.write('\t%.5f\t%s' % (lt,lps))

    if tests:
      out.write('\t%d' % df)

    ors = []
    for cat in range(c):
      beta = g.beta[cat*n:(cat+1)*n,0]
      W    = g.W[cat*n:(cat+1)*n,:][:,cat*n:(cat+1)*n]
      for orr,(ci_l,ci_u) in zip(model_term.odds_ratios(beta),
                                 model_term.odds_ratio_ci(beta,W)):
        ors.extend( [orr,ci_l,ci_u] )

    orsstr = '\t'.join('%.4f' % orr if isfinite(orr) else '' for orr in ors)
    out.write('\t%s' % orsstr)

    # FIXME: This does not align the categories -- just the blank headers
    if len(null.categories) < len(g.categories):
      short = len(g.categories)-len(null.categories)
      out.write( '\t'*(short*6) )

    out.write('\n')

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results(details,model,g)

      if 'score' in tests:
        details.write('Score test           : X2=%9.5f, df=%d, p=%s\n' % (st,df,sps))
      if 'wald' in tests:
        details.write('Wald test            : X2=%9.5f, df=%d, p=%s\n' % (wt,df,wps))
      if 'lrt' in tests:
        details.write('Likelihood ratio test: X2=%9.5f, df=%d, p=%s\n' % (lt,df,lps))

      details.write('\n')
      details.write('-'*79)
      details.write('\n')


if __name__ == '__main__':
  main()
