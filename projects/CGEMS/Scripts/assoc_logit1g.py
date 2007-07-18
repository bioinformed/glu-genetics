import csv
import sys

from   numpy                import abs,bmat,exp
from   scipy                import stats
from   scipy.linalg         import LinAlgError

from   biozilla.fileutils   import autofile,hyphen,load_list
from   biozilla.genodata    import load_genostream
from   biozilla.logit       import GLogit

from   biozilla.association import contingency_table,contingency_analysis,print_results,load_phenos, \
                                   LocusModelBuilder,NULL,GENO,ADOM,TREND,DOM,REC


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output summary results to FILE')
  parser.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output detailed results to FILE')
  parser.add_option('-p', '--detailsmaxp', dest='detailsmaxp', metavar='P', type='float', default=1.0,
                    help='Output detailed resutls for only pvalues below P threshold')
  parser.add_option('-d', '--dropsamples', dest='dropsamples', metavar='FILE',
                    help='List of samples to drop')
  parser.add_option('-D', '--droploci', dest='droploci', metavar='FILE',
                    help='List of loci to drop')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter (default=0.01)')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter (default=10)')

  return parser


def eval_model(lname,k,model,details,detailsmaxp):
  assert 0<=k<=2

  res = ['']*([1,9,12][k])

  if model is None:
    return res

  try:
    g = GLogit(model.y,model.X,add_mean=False,vars=model.vars)
    L,b,W = g.fit()
    ps = pw = 1
    if k:
      st,df = g.score_test(indices=range(1,k+1)).test()
      wt,df = g.wald_test(indices=range(1,k+1)).test()
      sf = stats.distributions.chi2.sf
      ps = sf(st,df)
      pw = sf(wt,df)
  except LinAlgError:
    return res

  if k == 2:
    res = ['%.6f'  % L,
           '%8.5f' % st, '%d' % df, '%9.8f' % ps,
           '%8.5f' % wt, '%d' % df, '%9.8f' % pw,
           '%7.4f' % exp(b[1,0]), '%7.4f' % exp(b[2,0]),
            '%.6f' % W[1,1], '%.6f' % W[2,2], '%.6f' % W[1,2]]
  elif k == 1:
    res = ['%.6f'  % L,
           '%8.5f' % st, '%d' % df, '%9.8f' % ps,
           '%8.5f' % wt, '%d' % df, '%9.8f' % pw,
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

  header,phenos = load_phenos(args[0])
  loci = load_genostream(args[1],options.format).as_ldat()

  if options.droploci:
    loci = loci.transformed(exclude_loci=options.droploci)

  if options.dropsamples:
    drop   = set(load_list(options.dropsamples))
    phenos = (p for p in phenos if p[0] not in drop)

  if loci.samples:
    samples = set(loci.samples)
    phenos  = (p for p in phenos if p[0] in samples)

  phenos = list(phenos)

  assert loci.samples
  models = LocusModelBuilder(loci.samples,header,phenos,loci.genorepr,
                             minmaf=options.minmaf,mingenos=options.mingenos)

  if options.nullmodel:
    null_model = models.build_model(NULL(), {})

    # Standardize all columns by 1/max(abs(column))
    X_null = null_model.X
    null_model.X = X_null.A/abs(X_null).max(axis=0)

    # Obtain estimates of covariate effects under the null
    null = GLogit(null_model.y,null_model.X,add_mean=False,vars=null_model.vars)
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

  loci    = iter(loci)
  samples = loci.next()

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

    all_models = [ model_null, model_geno, model_trend, model_dom, model_rec ]

    if all_models.count(None) == len(all_models):
      continue

    results = [lname]
    results.extend( eval_model(lname, 0, model_null,  details, options.detailsmaxp) )
    results.extend( eval_model(lname, 2, model_geno,  details, options.detailsmaxp) )
    results.extend( eval_model(lname, 2, model_adom,  details, options.detailsmaxp) )
    results.extend( eval_model(lname, 1, model_trend, details, options.detailsmaxp) )
    results.extend( eval_model(lname, 1, model_rec,   details, options.detailsmaxp) )
    results.extend( eval_model(lname, 1, model_dom,   details, options.detailsmaxp) )

    out.write('\t'.join(results))
    out.write('\n')


if __name__ == '__main__':
  if 0:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
  else:
    main()
