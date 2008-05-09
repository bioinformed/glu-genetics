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

from   numpy               import isfinite
from   scipy               import stats
from   scipy.linalg        import LinAlgError

from   glu.lib.fileutils   import autofile,hyphen
from   glu.lib.glm         import Linear
from   glu.lib.association import build_models,get_term,NULL,print_results_linear,format_pvalue


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

  loci,models = build_models(args[0], args[1], options, deptype=float)

  if options.nullmodel:
    null_model = models.build_model(NULL(),{})
    g = Linear(null_model.y,null_model.X,vars=null_model.vars)
    g.fit()
    print_results_linear(out,null_model,g)
    if options.details:
      print_results_linear(details,null_model,g)

  tests = set(t.strip().lower() for t in options.tests.split(','))
  tests.discard('')
  extra = tests - set(['score','wald','lrt'])
  if extra:
    raise ValueError('Unknown test(s) specified: %s' % ','.join(sorted(extra)))

  headers = ['Locus', 'Alleles', 'MAF', 'Geno Counts', 'Subjects']

  if 'score' in tests:
    headers += ['score X2', 'score p-value']
  if 'wald' in tests:
    headers += ['Wald X2',  'Wald p-value' ]
  if 'lrt' in tests:
    headers += ['LR X2',    'LR p-value'   ]

  if tests:
    headers += ['df']

  headers += ['Het Estimate', 'Het SE',
              'Hom Estimate', 'Hom SE']

  out.write('\t'.join(headers))
  out.write('\n')

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

    k = len(model_term)
    n = model.X.shape[0]

    g = Linear(model.y,model.X,vars=model.vars)
    try:
      g.fit()
    except LinAlgError:
      # Replace with something more robust ASAP
      continue

    m = model.model_loci[lname]
    out.write('\t'.join([lname, ','.join(m.alleles),
                                '%.3f' % m.maf,
                                '|'.join(map(str,m.counts)),
                                str(n)]))

    sp = wp = lp = 1
    if 'score' in tests:
      st,df = g.score_test(indices=model_term.indices()).test()
      sp    = stats.distributions.chi2.sf(st,df)
      sps   = format_pvalue(sp)
      out.write('\t%.5f\t%s' % (st,sps))

    if 'wald' in tests:
      wt,df = g.wald_test(indices=model_term.indices()).test()
      wp    = stats.distributions.chi2.sf(wt,df)
      wps   = format_pvalue(wp)
      out.write('\t%.5f\t%s' % (wt,wps))

    if 'lrt' in tests:
      lt,df = g.lr_test(indices=model_term.indices()).test()
      lp    = stats.distributions.chi2.sf(lt,df)
      lps   = format_pvalue(lp)
      out.write('\t%.5f\t%s' % (lt,lps))

    if tests:
      out.write('\t%d' % df)

    estimates = []
    for b,se in zip(model_term.estimates(g.beta),
                    model_term.se(g.W)):
      estimates.extend( [b,se] )

    estimatesstr = '\t'.join('%.4f' % e if isfinite(e) else '' for e in estimates)
    out.write('\t%s' % estimatesstr)

    out.write('\n')

    if options.details and min(sp,wp,lp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results_linear(details,model,g)

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
