import csv
import sys

from   numpy                import array,empty,zeros_like,exp,nan,abs,isfinite,vstack,vsplit
from   scipy                import stats

from   biozilla.fileutils   import autofile,hyphen,load_list
from   biozilla.genodata    import load_genostream
from   biozilla.logit       import GLogit

from   biozilla.association import contingency_table,contingency_analysis,table_trend,load_phenos,print_results, \
                                   LocusModelBuilder,get_term,NULL

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
  parser.add_option('--minmaf', dest='minmaf', metavar='N', default=0.01, type='float',
                    help='Minimum minor allele frequency filter')
  parser.add_option('--mingenos', dest='mingenos', metavar='N', default=10, type='int',
                    help='Minimum number of observed genotype filter')
  parser.add_option('--skipunadjusted', dest='skipunadjusted', action='store_true', default=True,
                    help='Skip unadjusted tests')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')
  parser.add_option('--genomodel', dest='genomodel', default='geno,trend', metavar='M1,M2,..',
                    help='Comma separated list of genetic models.  The first that can be fit will be used.  '
                         'Values: genotype/geno, adddom/adom, trend/multiplicative/mult, additive/add, '
                         'dominant/dom, recessive/rec, missing/miss, not_missing/not_miss, null.  Default=geno,trend')
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

  header,phenos = load_phenos(args[0])
  loci = load_genostream(args[1],options.format).as_ldat()

  if options.droploci:
    loci = loci.transformed(exclude_loci=options.droploci)

  if options.dropsamples:
    drop = set(load_list(options.dropsamples))
    phenos = (p for p in phenos if p[0] not in drop)

  if loci.samples:
    samples = set(loci.samples)
    phenos = (p for p in phenos if p[0] in samples)

  phenos = list(phenos)

  models = LocusModelBuilder(loci.samples,header,phenos,loci.genorepr,
                             minmaf=options.minmaf,mingenos=options.mingenos)

  null_model = models.build_model(NULL(),{})

  # Standardize all columns by 1/max(column)
  X_null = null_model.X
  null_model.X = X_null.A/abs(X_null).max(axis=0)

  # Obtain estimates of covariate effects under the null
  null = GLogit(null_model.y,null_model.X,add_mean=False,vars=null_model.vars)
  cats = len(null.categories)
  null.fit()

  if options.nullmodel:
    print_results(out,null_model,null)
    if options.details:
      print_results(details,null_model,null)

  headers = ['Locus']

  or_headers = []
  if len(null.categories) > 2:
    for i in range(1,len(null.categories)):
      or_headers += ['HetOR%d' % i, 'HomOR%d' % i]
  else:
    or_headers += ['HetOR', 'HomOR']

  if not options.skipunadjusted:
    headers += ['Unadjusted Score', 'DF', 'p-value']
    headers += or_headers

  headers += ['Adjusted Score', 'DF',  'p-value']
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

  loci    = iter(loci)
  samples = loci.next()

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

    y    = model.y
    X    = model.X
    vars = model.vars

    if 0:
      f = csv.writer(file('%s.csv' % lname,'w'))
      f.writerow(vars)
      f.writerows(X.tolist())

    n = X.shape[1]
    k = len(model_term)

    if not k or not X.shape[0]:
      continue
    elif k>2:
      raise ValueError,'Unexpected number of parameters in model (n=%d,k=%d)' % (n,k)

    # Construct genotype classes and genotype parameter indices
    X_trend = zeros_like(y)
    indices = []
    for i in range(1,k+1):
      X_trend += i*X[:,i]
      for j in range(cats-1):
        indices.append(j*n+i)

    indices.sort()

    out.write(lname)

    cp = 1
    if not options.skipunadjusted:
      ### Contingency table test ###
      counts = contingency_table(y,X_trend)

      if options.genomodel == 'unconstrained':
        ct,df_c,ors = contingency_analysis(counts)
      elif options.genomodel == 'trend':
        ct,df_c,ors = table_trend(counts),1,array([])
      else:
        raise ValueError('Unknown genotype model')

      cp = stats.distributions.chi2.sf(ct,df_c)
      out.write('\t%8.5f\t%d\t%9.7f' % (ct,df_c,cp))

      orsl = empty( (2*(cats-1),), dtype=float )
      orsl[:] = nan
      if ors.size:
        orsl[::3-ors.shape[1]] = ors.flatten()

      orsstr = '\t'.join('%7.4f' % orr if isfinite(orr) else '' for orr in orsl)
      out.write('\t%s' % orsstr)

    if 1:
      ### SCORE TEST ###
      g = GLogit(y,X,add_mean=False)
      st,df_s = g.score_test(indices=indices,initial_beta=null.beta).test()
      sp = stats.distributions.chi2.sf(st,df_s)
      out.write('\t%8.5f\t%d\t%9.7f' % (st,df_s,sp))

      g.fit(initial_beta=initial_beta[k])
      ors = exp(g.beta).take(indices).A.flatten()

      orsl = empty( (2*(cats-1),), dtype=float )
      orsl[:] = nan
      orsl[::3-k] = ors
      orsstr = '\t'.join('%7.4f' % orr if isfinite(orr) else '' for orr in orsl)
      out.write('\t%s' % orsstr)

    out.write('\n')

    if options.details and min(cp,sp) <= options.detailsmaxp:
      details.write('\nRESULTS: %s\n\n' % lname)
      print_results(details,model,g)
      if not options.skipunadjusted:
        details.write('Unadjusted Score Test: X2=%6.3f, df=%d, p=%12.10f\n' % (ct,df_c,cp))
      details.write('Adjusted   Score Test: X2=%6.3f, df=%d, p=%12.10f\n' % (st,df_s,sp))
      details.write('-'*79)
      details.write('\n')


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
