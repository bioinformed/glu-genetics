# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Fit single-SNP generalized logistic genotype-phenotype association models'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   itertools           import izip

import scipy.stats

from   numpy               import zeros,isfinite

from   glu.lib.fileutils   import autofile,hyphen,table_writer
from   glu.lib.glm         import GLogit,LinAlgError

from   glu.lib.genolib     import geno_options
from   glu.lib.association import build_models,print_results,format_pvalue,TREND,NO_INTERCEPT


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes [genotypes]'
  parser = optparse.OptionParser(usage=usage)

  input = optparse.OptionGroup(parser, 'Input options')

  geno_options(input,input=True,filter=True)

  input.add_option('--fixedloci', dest='fixedloci', metavar='FILE',
                    help='Genotypes for fixed loci that can be included in the model')
  input.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  input.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter.  default=10')
  input.add_option('-c', '--categorical', dest='categorical', metavar='VAR', action='append',
                    help='Create indicator variables based on values of VAR')
  input.add_option('--includevar', dest='includevar', metavar='VAR=VAL', action='append',
                    help='Include only records with variable VAR equal to VAL')
  input.add_option('--excludevar', dest='excludevar', metavar='VAR=VAL', action='append',
                    help='Exclude all records with variable VAR equal to VAL')

  analysis = optparse.OptionGroup(parser, 'Analysis options')

  analysis.add_option('--model', dest='model', metavar='FORMULA',
                      help='General formula for model to fit')
  analysis.add_option('--test', dest='test', metavar='FORMULA',
                      help='Formula terms to test.  Default is to test all genotype effects if a model is specified, '
                           'otherwise a 2df genotype test (GENO(locus)).')
  analysis.add_option('--stats', dest='stats', default='score', metavar='T1,T2,..',
                      help='Comma separated list of test statistics to apply to each model.  Supported tests '
                           'include score, Wald, and likelihood ratio statistics.  Values: score, wald, lrt.')
  analysis.add_option('--scan', dest='scan', metavar='NAME', default='locus',
                      help='Name of locus over which to scan, used in --model, --test and --display (default=locus)')
  analysis.add_option('--pid', dest='pid', metavar='NAME', default=0,
                      help='Subject column name or number in the phenotype file (default=1)')
  analysis.add_option('--pheno', dest='pheno', metavar='NAME', default=1,
                      help='Phenotype column name or number in the phenotype file (default=2)')
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
  output.add_option('--display', dest='display', metavar='FORMULA',
                      help='Formula terms to display in the summary output table.  Defaults to all test terms.')
  output.add_option('--detailsmaxp', dest='detailsmaxp', metavar='P', type='float', default=1.0,
                    help='Output detailed results for only pvalues below P threshold')
  output.add_option('-v', '--verbose', dest='verbose', metavar='LEVEL', type='int', default=1,
                    help='Verbosity level of diagnostic output.  O for none, 1 for some (default), 2 for exhaustive.')
  output.add_option('--se', dest='se', action='store_true', default=False,
                    help='Show standard errors around each estimate')
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


def summary_header(options,null):
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
  if len(null.categories) > 2:
    for i in range(1,len(null.categories)):
      for name in options.display.effect_names():
        if name.startswith('locus:'):
          name = name[6:]
        header.append( '%s%d OR' % (name,i) )
        if options.se:
          header += ['%s%d SE' % (name,i) ]
        if ci:
          header += [ '%s%d OR %d%% CI_l' % (name,i,ci),
                      '%s%d OR %d%% CI_u' % (name,i,ci) ]
  else:
    for name in options.display.effect_names():
      if name.startswith('locus:'):
        name = name[6:]
      header.append('%s OR' % name)
      if options.se:
        header += ['%s SE' % name]
      if ci:
        header += [ '%s OR %d%% CI_l' % (name,ci),
                    '%s OR %d%% CI_u' % (name,ci) ]

  return header


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) < 1:
    parser.print_help(sys.stderr)
    return

  phenos = args[0]
  genos  = args[1] if len(args)==2 else None

  if options.ci < 0 or options.ci > 1:
    raise ValueError('Confidence interval must be between 0 and 1')

  out = table_writer(options.output,hyphen=sys.stdout)
  if options.details:
    details = autofile(hyphen(options.details,sys.stdout),'w')
    if details is out:
      raise ValueError('Cannot send summary and detailed output to stdout')

  loci,fixedloci,gterms,models = build_models(phenos, genos, options)

  null_model = models.build_model(options.null,fixedloci)

  if not null_model:
    raise ValueError('Cannot construct null model')

  # Obtain estimates of covariate effects under the null
  null = GLogit(null_model.y,null_model.X,vars=null_model.vars)
  null.fit()

  if options.details:
    details.write('NULL MODEL:\n\n')
    print_results(details,null_model,null)

  if not genos:
    return

  header = summary_header(options,null)
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

    base_model = NO_INTERCEPT()+TREND(lname)
    for fixed in fixedloci:
      base_model += TREND(fixed)
    base = models.build_model(base_model,lmap,mingenos=0)

    if not base:
      out.writerow(result)
      continue

    counts = zeros( (len(null.categories),3), dtype=int )
    for pheno,geno in izip(base.y,base.X):
      counts[pheno[0,0],geno[0,0]] += 1

    mafs = (counts*[0.,0.5,1.]).sum(axis=1)/counts.sum(axis=1)

    m = base.model_loci[lname]
    result += ['|'.join(m.alleles),
               '|'.join('%.3f' % m for m in mafs),
               '/'.join('|'.join('%d' % c for c in row) for row in counts),
               '|'.join('%d' % c for c in counts.sum(axis=1))]

    model = models.build_model(options.model,lmap)

    if not model:
      out.writerow(result)
      continue

    g = GLogit(model.y,model.X,vars=model.vars)

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

    # R model verfication debugging code
    if 0:
      check_R(model,g)

    assert len(null.categories) >= len(g.categories)

    n = model.X.shape[1]

    # Construct genotype parameter indices
    test_indices = [ j*n+i for j in range(len(g.categories)-1)
                           for i in options.test.indices() ]

    sp = wp = lp = None
    df = ''

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

    res = []
    for categ in null.categories[1:]:
      try:
        cat  = g.categories.tolist().index(categ)-1
      except ValueError:
        res += ['']*(1+bool(options.se)+2*bool(options.ci))
      else:
        beta = g.beta[cat*n:(cat+1)*n,0]
        W    = g.W[cat*n:(cat+1)*n,:][:,cat*n:(cat+1)*n]

        ors  = options.display.odds_ratios(beta)
        ses  = options.display.standard_errors(W)
        cis  = options.display.odds_ratio_ci(beta,W,alpha=options.ci)

        for i in range(len(ors)):
          res.append(ors[i])
          if options.se:
            res.append(ses[i])
          if options.ci:
            res += list(cis[i])

    result.extend('%.4f' % v if isfinite(v) else '' for v in res)

    out.writerow(result)

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results(details,model,g)

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
