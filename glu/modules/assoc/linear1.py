# -*- coding: utf-8 -*-
'''
File:          linear1.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Fit linear genotype-phenotype association models

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   scipy               import stats
from   scipy.linalg        import LinAlgError

from   glu.lib.utils       import pick,tally
from   glu.lib.fileutils   import autofile,hyphen
from   glu.lib.glm         import Linear
from   glu.lib.association import build_models,get_term,NULL,print_results_linear,format_pvalue


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

  loci,models = build_models(args[0], args[1], options, deptype=float)

  if options.nullmodel:
    null_model = models.build_model(NULL(),{})
    g = Linear(null_model.y,null_model.X,vars=null_model.vars)
    g.fit()
    print_results_linear(out,null_model,g)
    if options.details:
      print_results_linear(details,null_model,g)

  out.write('\t'.join( ['Locus','Wald X2', 'df', 'p-value', 'Score X2', 'df', 'p-value'] ))
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

    g = Linear(model.y,model.X,vars=model.vars)
    try:
      g.fit()
      w,df = g.wald_test(indices=range(1,k+1)).test()
      s,df = g.score_test(indices=range(1,k+1)).test()
    except LinAlgError:
      # Replace with something more robust ASAP
      continue

    sf = stats.distributions.chi2.sf
    pw = format_pvalue(sf(w,df))
    ps = format_pvalue(sf(s,df))
    out.write('\t'.join( [lname,
                          '%6.3f' % w, '%d' % df ,pw,
                          '%6.3f' % s, '%d' % df ,ps] ))
    out.write('\n')

    if options.details and min(pw,ps) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results_linear(details,model,g)
      details.write('Wald Test:  X2=%6.3f, df=%d, p=%s\n' % (w,df,pw))
      details.write('Score Test: X2=%6.3f, df=%d, p=%s\n' % (s,df,ps))
      details.write('-'*79)
      details.write('\n')


if __name__ == '__main__':
  main()
