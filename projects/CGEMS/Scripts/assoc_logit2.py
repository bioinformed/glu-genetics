import csv
import sys

from   itertools            import islice,izip,chain

from   numpy                import abs,bmat,exp
from   scipy                import stats
from   scipy.linalg         import LinAlgError

from   biozilla.fileutils   import load_list
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
                    help='Output of summary report')
  parser.add_option('-O', '--details', dest='details', metavar='FILE',
                    help='Output of summary report')
  parser.add_option('-d', '--dropsamples', dest='dropsamples', metavar='FILE',
                    help='List of samples to drop')
  parser.add_option('-D', '--droploci', dest='droploci', metavar='FILE',
                    help='List of loci to drop')
  parser.add_option('--subset', dest='subset', metavar='FILE',
                    help='Only consider subset of 2-locus models around the list of SNPs in FILE')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('-w', '--window', dest='window', metavar='N', type='int', default=25,
                    help='Window size around each SNP')
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')
  return parser


def window(loci,n):
  '''
  >>> list(window(range(10),3))

  '''
  import collections

  items = collections.deque()
  loci = iter(loci)

  for locus in loci:
    items.append(locus)
    if len(items) == n+1:
      break

  l = list(items)
  yield [],l[0],l[1:]

  i = 0
  for locus in loci:
    if i==n:
      items.popleft()
    else:
      i+=1

    items.append(locus)
    l = list(items)
    yield l[:i],l[i],l[i+1:]

  l = list(items)
  for i in xrange(i,len(l)-1):
    i += 1
    yield l[i-n:i],l[i],l[i+1:]


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

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

  subset = None
  if options.subset:
    subset = set(load_list(options.subset))

  if options.nullmodel:
    null_model = models.build_model(NULL(), {})

    # Standardize all columns by 1/max(abs(column))
    X_null = null_model.X
    null_model.X = X_null.A/abs(X_null).max(axis=0)

    # Obtain estimates of covariate effects under the null
    null = GLogit(y,X_null,add_mean=False,vars=null_model.vars)
    null.fit()

    print_results(sys.stdout,null_model,null)

  headers = ['Ref_Locus',  'Alt_Locus',
             'Ref_Score',  'DF', 'p-value',
             'Alt_Score',  'DF', 'p-value',
             'Cond_Score', 'DF', 'p-value']

  print '\t'.join(headers)

  loci    = iter(loci)
  samples = loci.next()

  # For each locus
  for left,locus,right in window(loci,options.window):
    lname1 = locus[0]

    if subset is not None and lname1 not in subset:
      continue

    model1_term = GENO(lname1)
    model1 = models.build_model(model1_term, dict([locus]))

    if not model1:
      model1_term = TREND(lname1)
      model1 = models.build_model(model1_term, dict([locus]))

    if not model1:
      continue

    k1 = len(model1_term)

    try:
      g = GLogit(model1.y,model1.X,add_mean=False,vars=model1.vars)
      st1,df1 = g.score_test(indices=range(1,1+k1)).test()
    except LinAlgError:
      continue

    for other in chain(left,right):
      lmap = dict([locus,other])
      lname2 = other[0]

      model2_term = GENO(lname2)+NULL(lname1)
      model2 = models.build_model(model2_term, lmap)

      if not model2:
        model2_term = TREND(lname2)+NULL(lname1)
        model2 = models.build_model(model2_term, lmap)

      if not model2:
        continue

      k2 = len(model2_term)

      try:
        g = GLogit(model2.y,model2.X,add_mean=False,vars=model2.vars)
        st2,df2 = g.score_test(indices=range(1,1+k2)).test()
      except LinAlgError:
        continue

      interaction_term = TREND(lname1)*TREND(lname2)
      term = model1_term+model2_term+interaction_term
      model = models.build_model(term, lmap)

      if not model:
        interaction_term = NULL(lname1)+NULL(lname2)
        term = model1_term+model2_term+interaction_term
        model = models.build_model(model1_term+model2_term, lmap)

      if not model:
        continue

      k  = len(term)

      try:
        g = GLogit(model.y,model.X,add_mean=False,vars=model.vars)
        st,df = g.score_test(indices=range(1+k1,1+k)).test()
      except LinAlgError:
        continue

      sf = stats.distributions.chi2.sf
      print '\t'.join([lname1,lname2,
                       '%8.5f' % st1, '%d' % df1, '%9.8f' % sf(st1,df1),
                       '%8.5f' % st2, '%d' % df2, '%9.8f' % sf(st2,df2),
                       '%8.5f' % st,  '%d' % df,  '%9.8f' % sf(st, df)])


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
