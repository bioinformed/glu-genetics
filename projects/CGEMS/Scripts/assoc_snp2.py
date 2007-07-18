import csv
import sys
import random

from   numpy          import array,matrix,vstack,compress
from   numpy.linalg   import LinAlgError
from   itertools      import izip,imap,repeat,islice
from   operator       import itemgetter

from   biozilla.utils import pick,tally,autofile
from   biozilla.association import *

from   assoc_snp1     import load_phenos,build_genomap,print_results


def build_genomap2(locus1,locus2,genos1,genos2,geno_indices):
  genos1 = set(genos1)
  genos2 = set(genos2)

  genos = zip(pick(locus1,geno_indices.itervalues()),
              pick(locus2,geno_indices.itervalues()))

  genocounts = tally(g for g in genos if '' not in g and '  ' not in g)
  i1 = {}
  i2 = {}
  for g1,g2 in genocounts:
    i1.setdefault(g1,len(i1))
    i2.setdefault(g2,len(i2))

  i1 = sorted(i1.items(), key=itemgetter(1))
  i2 = sorted(i2.items(), key=itemgetter(1))

  if 0:
    print
    print '--------------------------------'
    print genocounts
    print i1
    print i2
    print map(itemgetter(0),i1)
    print map(itemgetter(0),i2)
    c = [ [g1]+[genocounts.get((g1,g2),0) for g2,i in i2 ] for g1,i in i1 ]
    print
    print c
    print genos1
    print genos2

  if len(i1) < 2 or len(i2) < 2:
    return {}, []

  if 1:
    gs1 = set(g[0] for g in genocounts)
    gs2 = set(g[1] for g in genocounts)
    base1 = (gs1 - genos1).pop()
    base2 = (gs2 - genos2).pop()

    dropg1 = set( g for g in gs1 if not genocounts.get( (g,base2) ) )
    dropg2 = set( g for g in gs2 if not genocounts.get( (base1,g) ) )
    genos1 -= dropg1
    genos2 -= dropg2

  genocounts = [ (g,n) for g,n in genocounts.iteritems() if g[0] in genos1 and g[1] in genos2 ]
  genocounts = sorted(genocounts, key=itemgetter(1))

  #genocounts = [ (g,n) for g,n in genocounts if n > 4 ]

  #print locus1[0],locus2[0],genocounts
  counts = [(g,i) for i,(g,n) in enumerate(genocounts)]
  return dict(counts),[ g for g,i in counts ]


def build_two_locus_models(loci,genodata,header,phenos):
  genodata = iter(genodata)
  geno_header = genodata.next()
  genodata = dict( (locus[0],locus) for locus in genodata)

  pids = [row[0] for row in phenos]
  covs = [ map(float,row[2:]) for row in phenos ]

  pidset = set(pids)
  geno_indices = dict( (pid,i) for i,pid in enumerate(geno_header) if i and pid in pidset )

  pids = [ row[0]             for row in phenos if row[0] in geno_indices ]
  stats= [ [int(row[1])]      for row in phenos if row[0] in geno_indices ]
  covs = [ map(float,row[2:]) for row in phenos if row[0] in geno_indices ]

  for lname1,lname2 in loci:
    locus1 = genodata[lname1]
    locus2 = genodata[lname2]

    genomap1,genonames1 = build_genomap(locus1,geno_indices)
    genomap2,genonames2 = build_genomap(locus2,geno_indices)

    if len(genomap1)<2 or len(genomap2)<2:
      print
      print '------------------'
      print genomap1
      print genomap2

    if len(genomap1)<2 or len(genomap2)<2:
      continue

    genomap12,genonames12 = build_genomap2(locus1,locus2,genonames1[:-1],genonames2[:-1],geno_indices)

    k1  = len(genomap1)-1
    k2  = len(genomap2)-1
    k12 = len(genomap12)

    #print genonames1[:-1],genonames2[:-1],genonames12

    X = []
    y = []
    for pid,stat,cov in izip(pids,stats,covs):
      g1 = locus1[geno_indices[pid]]
      g2 = locus2[geno_indices[pid]]
      g  = g1,g2

      if '' in g or '  ' in g:
        continue

      if g1 not in genomap1 or g2 not in genomap2:
        continue

      genos = [0]*(k1+k2+k12)
      i = genomap1[g1]
      if i!=k1:
        genos[i] = 1

      j = genomap2[g2]
      if j!=k2:
        genos[k1+j] = 1

      # Add interactions
      l = genomap12.get( (g1,g2), None )
      if l is not None:
        genos[k1+k2+l] = 1

      #print pid,g1,g2,' '.join(map(str,genos))

      y.append(stat)
      X.append( [1] + genos + cov )

    y = matrix(y, dtype=int)
    X = matrix(X, dtype=float)

    vars = (['_mean'] + [ '%s[%s]' % (lname1,g) for g,n in genonames1[:-1] ]
                      + [ '%s[%s]' % (lname2,g) for g,n in genonames2[:-1] ]
                      + header[2:])
    yield lname1,lname2,len(genonames1)-1,len(genonames2)-1,len(genomap12),vars,y,X


def geno_vectors(X,k,offset):
  n,m = X.shape
  offset2 = m+offset
  if k == 2:
    g1,g2   = X.T[offset].T==1,X.T[offset+1].T==1
    g       = -(g1|g2),g1,g2
    indices = (offset,offset+1,offset2,offset2+1)
  elif k == 1:
    g1,g2   = X.T[offset].T==1,0
    g       = -g1,g1,0
    indices = (offset,offset2)
  else:
    raise ValueError,'Unexpected number of parameters in model (k=%d)' % k

  return g,indices


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  parser.add_option('--skipunadjusted', dest='skipunadjusted', action='store_true', default=False, help='Skip unadjusted tests')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')

  return parser




def assoc_snp2(loci,genodata,phenoheader,phenos,out):
  #loci = [('rs1447295','rs7837688')]
  models = build_two_locus_models(loci,genodata,phenoheader,phenos)

  if 1:
    y      = matrix( [ [int(row[1])]         for row in phenos], dtype=int   )
    X_null = matrix( [[1]+map(float,row[2:]) for row in phenos], dtype=float )

    # Obtain estimates of covariate effects under the null
    L,b_null,W = logit3(y,X_null,add_mean=False)

    #vars = ['_mean']+header[2:]
    #print_results(vars,X_null,L,b_null,W)

    # Obtain estimates of mean with no covariates under the null
    L,b_null_nocov,W = logit3(y,X_null.take([0], axis=1),add_mean=False)

    #vars = ['_mean']
    #print_results(vars,X_null,L,b_null,W)

    # Build starting estimates for other models
    m,n = X_null.shape
    bg  = [None]
    bg += [ vstack( (b_null[0],[[0]]*k,b_null[1:n],
                     b_null[n],[[0]]*k,b_null[n+1:]) ) for k in range(1,5) ]

    bng = [None]
    bng+= [ vstack( (b_null_nocov[0],[[0]]*k,
                     b_null_nocov[1],[[0]]*k) ) for k in range(1,3) ]

  out.write('\t'.join( ['Ref Locus','Comp Locus',
                        'Contingency Ref', 'df',
                        'Contingency Comp', 'df',
                        'Stratified Contingency', 'df',
                        #'Two locus LR',
                        #'Two locus Score', 'df',
                       ] ))
  out.write('\n')

  # For each locus
  for lname1,lname2,k1,k2,k12,vars,y,X in models:
    #if lname2 != 'rs2026339':
    #  continue

    if 0:
      f = csv.writer(file('%s.csv' % lname,'w'))
      f.writerow(vars)
      f.writerows(X.tolist())

    if not k1 or not k2:
      continue

    n = X.shape[1]
    k = n-X_null.shape[1]

    g1,indices1 = geno_vectors(X,k1,1)
    g2,indices2 = geno_vectors(X,k2,1+k1)

    ys = y==0,y==1,y==2

    out.write('%s\t%s' % (lname1,lname2))
    df = 2*k

    gg1 = array((g1[1] + 2*g1[2]).flat)
    gg2 = array((g2[1] + 2*g2[2]).flat)
    yy  = array((ys[1] + 2*ys[2]).flat)

    if 1:
      ### Contingency table test ###
      counts = contingency_table3(yy,gg1)
      ct,df,ors = contingency_analysis(counts)
      #out.write('\t%8.5f\t%9.7f' % (ct,1-stats.chi2.cdf(ct,df)))
      out.write('\t%9.7f\t%d' % (1-stats.chi2.cdf(ct,df),df))

      counts = contingency_table3(yy,gg2)
      ct,df,ors = contingency_analysis(counts)
      #out.write('\t%8.5f\t%9.7f' % (ct,1-stats.chi2.cdf(ct,df)))
      out.write('\t%9.7f\t%d' % (1-stats.chi2.cdf(ct,df),df))

    if 1:
      ct = 0
      df = 0
      for i in range(k1+1):
        counts = contingency_table3(compress(gg1==i,yy), compress(gg1==i,gg2))
        t,d,ors = contingency_analysis(counts)
        #print t,d,1-stats.chi2.cdf(t,d)
        if t is not nan:
          ct += t
          df += d

      #print ct
      #out.write('\t%8.5f\t%9.7f' % (ct,1-stats.chi2.cdf(ct,df)))
      out.write('\t%9.7f\t%d' % (1-stats.chi2.cdf(ct,df),df))

    if 0:
      ### LR Model with covariates ###
      indices = range(k1+1)+range(k1+k2+k12+1,n)
      lr,df = lr_test3(y,X,indices)
      #df = 2*(k2+k12)
      #out.write('\t%8.5f\t%9.7f' % (lr,1-stats.chi2.cdf(lr,df)))
      out.write('\t%9.7f' % 1-stats.chi2.cdf(lr,df))

    if 0:
      ### Model without covariates ###
      X_nocov = X.take( range(0,k+1), axis=1 )
      lr_nocov,df = lr_test3(y, X_nocov, [0], b_initial=bng[k], b_initial_null=b_null_nocov)
      out.write('\t%8.5f\t%9.7f' % (lr_nocov,1-stats.chi2.cdf(lr_nocov,df)))

    if 0:
      ### SCORE TEST ###

      # Obtain estimates of covariate effects under the null
      indices2 = (range(  k1+1,  k1+k2+k12+1)
               +  range(n+k1+1,n+k1+k2+k12+1))
      indices = [i for i in range(X.shape[1]) if i not in indices2]

      X_null = X.take(indices, axis=1)
      L,b_null,W = logit3(y,X_null,add_mean=False)

      b_null = list(b_null.flat)
      b_null[k1+1:k1+1]     = [0]*(k2+k12)
      b_null[n+k1+1:n+k1+1] = [0]*(k2+k12)
      b_null = matrix(b_null).T
      #print b_null.T

      s1,s2,w11,w12,w22 = build_scores3(ys,X,b_null)
      st,df = do_score_test3(s1,s2,w11,w12,w22,X,g2,indices2)
      #out.write('\t%8.5f\t%9.7f' % (st,1-stats.chi2.cdf(st,df)))
      out.write('\t%9.7f' % 1-stats.chi2.cdf(st,df))

    # FIXME: Permutation is not right for loci in LD
    if 0:
      n = 1000
      m = 0
      for y,X in islice(perms,n):
        g1,indices1 = geno_vectors(X,k1,1)
        g2,indices2 = geno_vectors(X,k2,1+k1)
        stp,df = do_score_test3(s1,s2,w11,w12,w22,X,g,indices)
        if stp >= st:
          m+= 1
      out.write('\t%9.7f' % (float(m)/n))

    if 0:
      out.write('\t%d' % df)

    out.write('\n')


def main():
  parser = option_parser()
  options,args = parser.parse_args()
  options.skipunadjusted = True

  if len(args) != 2:
    parser.print_help()
    return

  header,phenos = load_phenos(args[0])
  phenos = list(phenos)
  genodata = list(csv.reader(autofile(args[1]),dialect='excel-tab'))
  reflocus = args[1].split('/')[-1].split('.')[0]
  loci = ( (reflocus,locus[0]) for locus in islice(genodata,1,None) if locus[0]!=reflocus )

  assoc_snp2(loci,genodata,header,phenos,sys.stdout)


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
