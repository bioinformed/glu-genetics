# -*- coding: utf-8 -*-
'''
File:          linear1.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Fit linear genotype-phenotype association models

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

import scipy.stats

from   numpy               import isfinite

from   glu.lib.fileutils   import autofile,hyphen,table_writer
from   glu.lib.glm         import Linear,LinAlgError

from   glu.lib.association import build_models,print_results_linear,format_pvalue, \
                                  GENO,COMBINATION,PHENOTERM,GENOTERM,FormulaParser


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  input = optparse.OptionGroup(parser, 'Input options')

  input.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Input genotype format')
  input.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  input.add_option('-i', '--includesamples', dest='includesamples', metavar='FILE',
                    help='List of samples to include')
  input.add_option('-d', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='List of samples to exclude')
  input.add_option('-I', '--includeloci', dest='includeloci', metavar='FILE',
                    help='List of loci to include')
  input.add_option('-D', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='List of loci to exclude')
  input.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  input.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter.  default=10')

  analysis = optparse.OptionGroup(parser, 'Analysis options')

  analysis.add_option('--model', dest='model', metavar='F',
                      help='General formula for model to fit')
  analysis.add_option('--test', dest='test', metavar='F',
                      help='Formula terms to test.  Default is to test all genotype effects if a model is specified, '
                           'otherwise a 2df genotype test (GENO(locus)).')
  analysis.add_option('--stats', dest='stats', default='score', metavar='T1,T2,..',
                      help='Comma separated list of test statistics to apply to each model.  Supported tests '
                           'include score, Wald, and likelihood ratio statistics.  Values: score, wald, lrt.')
  analysis.add_option('--refalleles', dest='refalleles', metavar='FILE',
                      help='Mapping of locus name to the corresponding reference allele')
  analysis.add_option('--allowdups', dest='allowdups', action='store_true', default=False,
                      help='Allow duplicate individuals in the data (e.g., to accommodate weighting '
                           'or incidence density sampling)')

  output = optparse.OptionGroup(parser, 'Output options')

  output.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output summary results to FILE')
  output.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output detailed results to FILE')
  output.add_option('--display', dest='display', metavar='F',
                      help='Formula terms to display in the summary output table.  Defaults to all test terms.')
  output.add_option('-p', '--detailsmaxp', dest='detailsmaxp', metavar='P', type='float', default=1.0,
                    help='Output detailed results for only pvalues below P threshold')
  output.add_option('-v', '--verbose', dest='verbose', metavar='LEVEL', type='int', default=1,
                    help='Verbosity level of diagnostic output.  O for none, 1 for some (default), 2 for exhaustive.')
  output.add_option('--ci', dest='ci', default=0, type='float', metavar='N',
                    help='Show confidence interval around each estimate of width N.  Set to zero to inhibit '
                         'output.  Default=0')

  parser.add_option_group(input)
  parser.add_option_group(analysis)
  parser.add_option_group(output)

  return parser


def check_R(model,g):
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

  coef  = r.coefficients(mod)
  coef  = array([coef['(Intercept)']] + [ coef[v] for v in vars ],dtype=float)
  coef2 = g.beta.A.reshape(-1)

  #assert allclose(coef,g.beta.A.reshape(-1),atol=1e-6)
  #assert allclose(r.vcov(mod),g.W,atol=1e-6)


def parse_formulae(options,models):
  if options.model is not None:
    pheno,options.model = FormulaParser().parse(options.model)

  if options.test is not None:
    _,options.test = FormulaParser().parse(options.test)

  if options.display is not None:
    _,options.display = FormulaParser().parse(options.display)

  if options.model and not options.test:
    options.test = COMBINATION(t for t in options.model.terms() if not t.loci())

  if not options.test:
    options.test = GENO('locus')

  if not options.model:
    phenos = COMBINATION( PHENOTERM(pheno) for pheno in models.pheno_header[2:] )
    options.model = options.test + phenos

  if not options.display:
    options.display = options.test

  try:
    options.test = options.model.find(options.test)
  except KeyError:
    raise ValueError('Formula does not contain all terms to be tested')

  try:
    options.display = options.model.find(options.display)
  except KeyError:
    raise ValueError('Formula does not contain all terms to display')

  options.null = COMBINATION(t for t in options.model.terms()
                                if t not in options.test and not t.loci())

  if not any(t for t in options.test.terms() if t.loci()):
    raise ValueError('Test does not contain any locus terms')

  if any(t for t in options.null.terms() if t.loci()):
    raise ValueError('Null model may not contain loci to be scanned')

  options.stats = set(t.strip().lower() for t in options.stats.split(','))
  options.stats.discard('')
  extra = options.stats - set(['score','wald','lrt'])
  if extra:
    raise ValueError('Unknown test(s) specified: %s' % ','.join(sorted(extra)))


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
    header += [ '%s' % name, '%s SE' % name ]
    if ci:
      header += [ '%s %d%% CI_l' % (name,ci),
                  '%s %d%% CI_u' % (name,ci) ]

  return header


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  if options.ci < 0 or options.ci > 1:
    raise ValueError('Confidence interval must be between 0 and 1')

  out = table_writer(options.output,hyphen=sys.stdout)
  if options.details:
    details = autofile(hyphen(options.details,sys.stdout),'w')
    if details is out:
      raise ValueError('Cannot send summary and detailed output to stdout')

  loci,models = build_models(args[0], args[1], options, deptype=float)

  parse_formulae(options,models)

  if options.details:
    details.write('NULL MODEL:\n\n')

    null_model = models.build_model(options.null,{})

    if not null_model:
      raise ValueError('Cannot construct null model')

    # Obtain estimates of covariate effects under the null
    null = Linear(null_model.y,null_model.X,vars=null_model.vars)
    null.fit()

    print_results_linear(details,null_model,null)

  header = summary_header(options)
  out.writerow(header)

  gterms = []
  for t in options.model.expand_terms():
    if isinstance(t, GENOTERM) and t.name=='locus':
      gterms.append(t)

  if not gterms:
    raise ValueError('No genotype terms to scan')

  # For each locus
  for locus in loci:
    lname = locus[0]
    lmap  = dict([locus])

    for t in gterms:
      t.name = lname

    model = models.build_model(options.model,lmap)
    if not model:
      continue

    g = Linear(model.y,model.X,vars=model.vars)

    try:
      g.fit()
    except LinAlgError:
      # FIXME: Output bad fit line
      continue

    # Desgin matrix debugging output
    if 0:
      f = table_writer('%s.csv' % lname,dialect='csv')
      f.writerow(model.vars)
      f.writerows(model.X.tolist())

    # R model verfication debugging code
    if 0:
      check_R(model,g)

    m = model.model_loci[lname]
    n = model.X.shape[0]

    result = [lname,
              ','.join(m.alleles),
              '%.3f' % m.maf,
              '|'.join(map(str,m.counts)),
              str(n) ]

    # Construct genotype parameter indices
    test_indices = options.test.indices()

    sp = wp = lp = 1

    if 'score' in options.stats:
      st,df = g.score_test(indices=test_indices).test()
      sp    = scipy.stats.distributions.chi2.sf(st,df)
      sps   = format_pvalue(sp)
      result.extend( ['%.5f' % st, sps ] )

    if 'wald' in options.stats:
      wt,df = g.wald_test(indices=test_indices).test()
      wp    = scipy.stats.distributions.chi2.sf(wt,df)
      wps   = format_pvalue(wp)
      result.extend( ['%.5f' % wt, wps ] )

    if 'lrt' in options.stats:
      lt,df = g.lr_test(indices=test_indices).test()
      lp    = scipy.stats.distributions.chi2.sf(lt,df)
      lps   = format_pvalue(lp)
      result.extend( ['%.5f' % lt, lps ] )

    if options.stats:
      result.append('%d' % df)

    # FIXME: SEs are optional
    display = options.display
    estimates = []
    if options.ci:
      for b,se,(ci_l,ci_u) in zip(display.estimates(g.beta),
                                  display.se(g.W),
                                  display.estimate_ci(g.beta,g.W,options.ci)):
        estimates.extend( [b,se,ci_l,ci_u] )
    else:
      for b,se in zip(display.estimates(g.beta),
                                  display.se(g.W)):
        estimates.extend( [b,se] )

    result.extend('%.4f' % e if isfinite(e) else '' for e in estimates)

    out.writerow(result)

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results_linear(details,model,g)

      if options.stats:
        details.write('Testing: %s\n\n' % options.test.formula())

      if 'score' in options.stats:
        details.write('Score test           : X2=%9.5f, df=%d, p=%s\n' % (st,df,sps))
      if 'wald' in options.stats:
        details.write('Wald test            : X2=%9.5f, df=%d, p=%s\n' % (wt,df,wps))
      if 'lrt' in options.stats:
        details.write('Likelihood ratio test: X2=%9.5f, df=%d, p=%s\n' % (lt,df,lps))

      details.write('\n')
      details.write('-'*79)
      details.write('\n')


if __name__ == '__main__':
  main()
