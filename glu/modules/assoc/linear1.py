# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Fit single-SNP linear genotype-phenotype association models'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

import scipy.stats

from   numpy               import isfinite

from   glu.lib.fileutils   import autofile,hyphen,table_writer,table_options
from   glu.lib.glm         import Linear,LinAlgError

from   glu.lib.genolib     import geno_options
from   glu.lib.association import build_models,print_results_linear,format_pvalue


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('phenotypes',            help='tabular or delimited file of phenotypes and covariates')
  parser.add_argument('genotypes',  nargs='?', help='genotype file in any format supported by GLU')

  input = parser.add_argument_group('Input options')

  geno_options(input,input=True,filter=True)

  input.add_argument('--fixedloci', metavar='FILE',
                    help='Genotypes for fixed loci that can be included in the model')
  input.add_argument('--minmaf', metavar='N', default=0.01, type=float,
                    help='Minimum minor allele frequency filter')
  input.add_argument('--mingenos', metavar='N', default=10, type=int,
                    help='Minimum number of observed genotype filter.  default=10')

  table_options(input)

  analysis = parser.add_argument_group('Analysis options')

  analysis.add_argument('--model', metavar='FORMULA',
                      help='General formula for model to fit')
  analysis.add_argument('--test', metavar='FORMULA',
                      help='Formula terms to test.  Default is to test all genotype effects if a model is specified, '
                           'otherwise a 2df genotype test (GENO(locus)).')
  analysis.add_argument('--stats', default='score', metavar='T1,T2,..',
                      help='Comma separated list of test statistics to apply to each model.  Supported tests '
                           'include score, Wald, and likelihood ratio statistics.  Values: score, wald, lrt.')
  analysis.add_argument('--scan', metavar='NAME', default='locus',
                      help='Name of locus over which to scan, used in --model, --test and --display (default=locus)')
  analysis.add_argument('--pid', metavar='NAME', default='1',
                      help='Column name or number of subject in the phenotype file (default=1)')
  analysis.add_argument('--pheno', metavar='NAME', default='2',
                      help='Phenotype column name or number in the phenotype file (default=2)')
  analysis.add_argument('--refalleles', metavar='FILE',
                      help='Mapping of locus name to the corresponding reference allele')
  analysis.add_argument('--allowdups', action='store_true', default=False,
                      help='Allow duplicate individuals in the data (e.g., to accommodate weighting '
                           'or incidence density sampling)')

  output = parser.add_argument_group('Output options')

  output.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output summary results to FILE')
  output.add_argument('-O', '--details', metavar='FILE',
                    help='Output detailed results to FILE')
  output.add_argument('--display', metavar='FORMULA',
                      help='Formula terms to display in the summary output table.  Defaults to all test terms.')
  output.add_argument('--detailsmaxp', metavar='P', type=float, default=1.0,
                    help='Output detailed results for only pvalues below P threshold')
  output.add_argument('-v', '--verbose', metavar='LEVEL', type=int, default=1,
                    help='Verbosity level of diagnostic output.  O for none, 1 for some (default), 2 for exhaustive.')
  output.add_argument('--se', action='store_true', default=False,
                      help='Show standard errors around each estimate')
  output.add_argument('--ci', default=0, type=float, metavar='N',
                    help='Show confidence interval around each estimate of width N.  Set to zero to inhibit '
                         'output.  Default=0')

  return parser


def check_R(model,g):
  import rpy
  from   rpy   import r
  from   numpy import array,allclose

  vars = [ v.replace(':','.').replace('+','p').replace('-','m').replace('_','.') for v in model.vars[1:] ]
  frame = dict( (v,model.X[:,i+1].reshape(-1)) for i,v in enumerate(vars) )
  frame['y'] = model.y.reshape(-1)
  formula = 'y ~ ' + ' + '.join(v.replace(':','.') for v in vars)

  rpy.set_default_mode(rpy.NO_CONVERSION)
  mod = r.glm(r(formula),data=r.data_frame(**frame),family=r.binomial('logit'))
  rpy.set_default_mode(rpy.BASIC_CONVERSION)
  pmod = mod.as_py()

  coef  = r.coefficients(mod)
  coef  = array([coef['(Intercept)']] + [ coef[v] for v in vars ],dtype=float)
  coef2 = g.beta.reshape(-1)

  #assert allclose(coef,g.beta.reshape(-1),atol=1e-6)
  #assert allclose(r.vcov(mod),g.W,atol=1e-6)


def summary_header(options):
  ci = int(100*options.ci) if options.ci else None

  header = ['Locus', 'Alleles', 'MAF', 'Geno Counts', 'Subjects']

  if 'score' in options.stats:
    header += ['score X2', 'score p-value']
  if 'wald' in options.stats:
    header += ['Wald X2',  'Wald p-value' ]
  if 'lrt' in options.stats:
    header += ['LR X2',    'LR p-value'   ]
  if options.stats:
    header += ['df']

  # FIXME: SEs are optional
  for name in options.display.effect_names():
    if name.startswith('locus:'):
      name = name[6:]
    header += [ '%s' % name ]
    if options.se:
      header += [ '%s SE' % name ]
    if ci:
      header += [ '%s %d%% CI_l' % (name,ci),
                  '%s %d%% CI_u' % (name,ci) ]

  return header


def main():
  parser  = option_parser()
  options = parser.parse_args()
  phenos  = options.phenotypes
  genos   = options.genotypes

  if not (0<=options.ci<=1):
    parser.error('Confidence interval must be between 0 and 1')

  out = table_writer(options.output,hyphen=sys.stdout)
  if options.details:
    details = autofile(hyphen(options.details,sys.stdout),'w')
    if details is out:
      parser.error('Summary and detailed output cannot both be sent to stdout')

  loci,fixedloci,gterms,models = build_models(phenos, genos, options, deptype=float)

  if options.details:
    details.write('NULL MODEL:\n\n')

    null_model = models.build_model(options.null,fixedloci)

    if not null_model or not null_model.valid(minmaf=options.minmaf, mingenos=options.mingenos):
      raise ValueError('Cannot construct null model')

    # Obtain estimates of covariate effects under the null
    null = Linear(null_model.y,null_model.X,vars=null_model.vars)
    null.fit()

    print_results_linear(details,null_model,null)

  if not genos:
    return

  header = summary_header(options)
  out.writerow(header)

  # For each locus
  for lname,genos in loci:
    result = [lname]

    # Skip fixed terms
    if lname in fixedloci:
      out.writerow(result)
      continue

    lmap = fixedloci.copy()
    lmap[lname] = genos

    for t in gterms:
      t.name = lname

    model = models.build_model(options.model,lmap)

    if not model:
      out.writerow(result)
      continue

    m = model.model_loci[lname]
    n = model.X.shape[0]

    result += [','.join(m.alleles),
               '%.3f' % m.maf,
               '|'.join(map(str,m.counts)),
               str(n) ]

    if not model.valid(minmaf=options.minmaf, mingenos=options.mingenos):
      out.writerow(result)
      continue

    g = Linear(model.y,model.X,vars=model.vars)

    try:
      g.fit()
    except LinAlgError:
      out.writerow(result)
      continue

    # Design matrix debugging output
    if 0:
      f = table_writer('%s.csv' % lname,dialect='csv')
      f.writerow(model.vars)
      f.writerows(model.X.tolist())

    # R model verification debugging code
    if 0:
      check_R(model,g)

    # Construct genotype parameter indices
    test_indices = options.test.indices()

    sp = wp = lp = None
    df = None

    if 'score' in options.stats:
      try:
        st,df = g.score_test(indices=test_indices).test()
      except LinAlgError:
        result.extend( ['',''] )
      else:
        sp    = scipy.stats.distributions.chi2.sf(st,df)
        sps   = format_pvalue(sp)
        result.extend( ['%.5f' % st, sps ] )

    if 'wald' in options.stats:
      try:
        wt,df = g.wald_test(indices=test_indices).test()
      except LinAlgError:
        result.extend( ['',''] )
      else:
        wp    = scipy.stats.distributions.chi2.sf(wt,df)
        wps   = format_pvalue(wp)
        result.extend( ['%.5f' % wt, wps ] )

    if 'lrt' in options.stats:
      try:
        lt,df = g.lr_test(indices=test_indices).test()
      except LinAlgError:
        result.extend( ['',''] )
      else:
        lp    = scipy.stats.distributions.chi2.sf(lt,df)
        lps   = format_pvalue(lp)
        result.extend( ['%.5f' % lt, lps ] )

    if options.stats:
      result.append(df)

    ors  = options.display.odds_ratios(g.beta)
    ses  = options.display.standard_errors(g.W)
    cis  = options.display.odds_ratio_ci(g.beta,g.W,alpha=options.ci)

    res = []
    for i in range(len(ors)):
      res.append(ors[i])
      if options.se:
        res.append(ses[i])
      if options.ci:
        res += list(cis[i])

    result.extend('%.4f' % e if isfinite(e) else '' for e in res)

    out.writerow(result)

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results_linear(details,model,g)

      if options.stats:
        details.write('Testing: %s\n\n' % options.test.formula())

      if 'score' in options.stats and sp is not None:
        details.write('Score test           : X2=%9.5f, df=%d, p=%s\n' % (st,df,sps))
      if 'wald' in options.stats and wp is not None:
        details.write('Wald test            : X2=%9.5f, df=%d, p=%s\n' % (wt,df,wps))
      if 'lrt' in options.stats and lp is not None:
        details.write('Likelihood ratio test: X2=%9.5f, df=%d, p=%s\n' % (lt,df,lps))

      details.write('\n')
      details.write('-'*79)
      details.write('\n')


if __name__ == '__main__':
  main()
