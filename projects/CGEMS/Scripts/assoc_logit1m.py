import csv
import sys

from   itertools            import islice

from   numpy                import abs,bmat
from   scipy                import stats

from   biozilla.fileutils   import load_list
from   biozilla.genodata    import load_genostream
from   biozilla.logit       import GLogit,exp

from   biozilla.association import contingency_table,contingency_analysis,print_results,load_phenos, \
                                   LocusModelBuilder,TREND,NULL


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=ldat, sdat, hapmap, genotriple')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  parser.add_option('-d', '--dropsamples', dest='dropsamples', metavar='FILE',
                    help='List of samples to drop')
  parser.add_option('-D', '--droploci', dest='droploci', metavar='FILE',
                    help='List of loci to drop')
  parser.add_option('-x', '--fixedloci', dest='fixedloci', metavar='FILE',
                    help='List of loci to include in every model')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  header,phenos = load_phenos(args[0])
  loci = load_genostream(args[1],options.format).as_ldat()

  fixed_loci = []
  if options.fixedloci:
    fixed_loci = islice(load_genostream(options.fixedloci,options.format).as_ldat(),1,None)

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
    null_model = models.build_model(NULL(),{})

    # Standardize all columns by 1/max(column)
    X_null = null_model.X
    null_model.X = X_null.A/abs(X_null).max(axis=0)

    # Obtain estimates of covariate effects under the null
    null = GLogit(null_model.y,null_model.X,add_mean=False,vars=null_model.vars)
    null.fit()

    print_results(sys.stdout,null_model,null)

  headers = ['Locus', 'Score1', 'DF', 'p-value', 'OR',
                      'Score2', 'DF', 'p-value', 'OR',
                      'Score3', 'DF', 'p-value', 'OR1', 'OR2', 'OR3',
                      'Score4', 'DF', 'p-value']

  print '\t'.join(headers)

  loci    = iter(loci)
  samples = loci.next()

  fixed_loci  = dict(fixed_loci)
  fixed_trend = sum( (TREND(l) for l in fixed_loci), NULL())
  fixed_null  = sum( ( NULL(l) for l in fixed_loci), NULL())

  # For each locus
  for locus in loci:
    lname,genos = locus

    if lname in fixed_loci:
      continue

    lmap = fixed_loci.copy()
    lmap[lname] = genos

    scan_model = TREND(lname)
    model1 = models.build_model(fixed_null +scan_model, lmap)
    model2 = models.build_model(fixed_trend+scan_model, lmap)
    model3 = models.build_model(fixed_trend+scan_model+fixed_trend*scan_model, lmap)

    if None in (model1,model2,model3):
      continue

    sys.stdout.write(lname)

    # Compute 1 df score test w/o 8q24 SNPs
    g = GLogit(model1.y,model1.X,add_mean=False)
    g.fit()
    st,df = g.score_test(indices=[1]).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f' % exp(g.beta[1,0]))

    # Compute 1 df score test w/ 8q24 SNPs
    g = GLogit(model2.y,model2.X,add_mean=False)
    g.fit()
    st,df = g.score_test(indices=[3]).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f' % exp(g.beta[3,0]))

    # Compute 2 df score test on just the interaction terms
    g = GLogit(model3.y,model3.X,add_mean=False)
    g.fit()
    st,df = g.score_test(indices=(4,5)).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\t%.3f\t%.3f\t%.3f' % (exp(g.beta[3,0]),exp(g.beta[4,0]),exp(g.beta[5,0])))

    # Compute 3 df score test on the main effect and interaction terms
    st,df = g.score_test(indices=(3,4,5)).test()
    sys.stdout.write('\t%8.5f\t%d\t%9.7f' % (st,df,stats.distributions.chi2.sf(st,df)))
    sys.stdout.write('\n')


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
