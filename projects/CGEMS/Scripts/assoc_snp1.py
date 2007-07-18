import csv
import sys
import random

from   numpy          import array,matrix,exp,vsplit,zeros,where,hstack,vstack,log,exp,maximum,nan,abs
from   itertools      import izip,imap,repeat,islice

from   biozilla.utils import pick,tally,autofile
from   biozilla.logit import logit3

from   biozilla.association    import *


def normalize(row,n):
  m = len(row)

  if n == m:
    return row
  elif n < m:
    return row[:n]
  else:
    return row + ['']*(n-m)


def load_phenos(filename):
  phenos = csv.reader(autofile(filename),dialect='excel-tab')
  header = phenos.next()
  return header,imap(normalize, phenos, repeat(len(header)))


def compute_maf(counts):
  hom1 = hom2 = het = 0

  for g,n in counts.iteritems():
    if g[0] != g[1]:
      het = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  n = (hom1+het+hom2)*2

  if not n:
    return 0.0

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)

  maf = float(2*hom1 + het)/n
  return maf


def build_genomap(locus,geno_indices,minmaf=0.01,mingenos=15):
  genocounts = tally(pick(locus,geno_indices.itervalues()))

  genocounts.pop('',None)
  genocounts.pop('  ',None)

  maf = compute_maf(genocounts)

  #print locus[0],maf

  if maf < minmaf:
    {},None

  genocounts = sorted( ( (n,g) for g,n in genocounts.iteritems() if n>=mingenos ), reverse=True)

  if len(genocounts) < 2:
    {},None

  homs = [ g for (n,g) in genocounts if g[0]==g[1]]
  hets = [ g for (n,g) in genocounts if g[0]!=g[1]]

  if len(homs) + len(hets) < 2:
    {},None

  testgenos = homs[:1] + hets + homs[1:]
  genomap = dict( (g,i) for i,g in enumerate(testgenos) )

  return genomap,testgenos


def build_locus_models(loci,header,phenos,minmaf=0.01,mingenos=15):
  loci = iter(loci)
  geno_header = loci.next()

  pids = [row[0] for row in phenos]
  covs = [ map(float,row[2:]) for row in phenos ]

  pidset = set(pids)
  #random.shuffle(geno_header)
  geno_indices = dict( (pid,i) for i,pid in enumerate(geno_header) if i and pid in pidset )

  pids = [ row[0]             for row in phenos if row[0] in geno_indices ]
  stats= [ [int(row[1])]      for row in phenos if row[0] in geno_indices ]
  covs = [ map(float,row[2:]) for row in phenos if row[0] in geno_indices ]

  for locus in loci:
    lname = locus[0]

    genomap,testgenos = build_genomap(locus,geno_indices,minmaf,mingenos)

    if len(genomap) < 2:
      continue

    k = max(genomap.itervalues())

    X = []
    y = []

    for pid,stat,cov in izip(pids,stats,covs):
      g = locus[geno_indices[pid]]

      if not g.strip():
        continue

      if g not in genomap:
        continue

      genos = [0]*k

      i = genomap[g]

      if i:
        genos[i-1] = 1

      y.append(stat)
      X.append( [1] + genos + cov )

    y = matrix(y, dtype=int)
    X = matrix(X, dtype=float)

    vars = ['_mean'] + testgenos[:-1] + header[2:]
    yield lname,vars,y,X


def print_results(vars,X,L,b,W):
  b = b.T
  stde = W.diagonal().A**.5
  z = b.A/stde
  oddsr = exp(b)

  print 'Multinomial logistic regression'
  print 'Observations =',len(X)
  print 'Log likelihood =',L
  print
  print 'TYPE       Variable             Coef.   Std.Err.     OR       Z'
  print '---------- -------------- ----------- ---------- ---------- -----'
  for i,cat in enumerate(['EARLY','ADVANCED']):
    i = i*len(vars)
    for j,var in enumerate(vars):
      k = j+i
      print '%-10s %-15s %10.6f %10.6f %10.6f %5.2f' % (cat,var,b[0,k],stde[0,k],oddsr[0,k],z[0,k])
    print
  print


def option_parser():
  import optparse

  usage = 'usage: %prog [options] phenotypes genotypes'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  parser.add_option('--skipunadjusted', dest='skipunadjusted', action='store_true', default=False, help='Skip unadjusted tests')
  parser.add_option('--nullmodel', dest='nullmodel', action='store_true', default=False, help='Show null model')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()
  options.skipunadjusted = False

  if len(args) != 2:
    parser.print_help()
    return

  header,phenos = load_phenos(args[0])
  loci   = csv.reader(autofile(args[1]),dialect='excel-tab')

  phenos = list(phenos)
  models = build_locus_models(loci,header,phenos)

  y      = matrix( [ [int(row[1])]         for row in phenos], dtype=int   )
  X_null = matrix( [[1]+map(float,row[2:]) for row in phenos], dtype=float )

  # Standardize all columns by 1/max(column)[B
  X_null = X_null.A/abs(X_null).max(axis=0)

  # Obtain estimates of covariate effects under the null
  L,b_null,W = logit3(y,X_null,add_mean=False)

  if 1: #options.nullmodel:
    vars = ['_mean']+header[2:]
    print_results(vars,X_null,L,b_null,W)

  # Obtain estimates of mean with no covariates under the null
  L,b_null_nocov,W = logit3(y,X_null.take([0], axis=1),add_mean=False)

  #vars = ['_mean']
  #print_results(vars,X_null,L,b_null,W)

  # Build starting estimates for other models
  m,n = X_null.shape
  bg  = [None]
  bg += [ vstack( (b_null[0],[[0]]*k,b_null[1:n],
                   b_null[n],[[0]]*k,b_null[n+1:]) ) for k in range(1,3) ]

  bng = [None]
  bng+= [ vstack( (b_null_nocov[0],[[0]]*k,
                   b_null_nocov[1],[[0]]*k) ) for k in range(1,3) ]

  if 0:
    print '\t'.join( ['Locus','Contingency Table',     'p-value',
                              'LR no covariates',      'p-value',
                              'LR with covariates',    'p-value',
                              'Score with covariates', 'p-value',
                              'DF'] )
  else:
    print '\t'.join( ['Locus','Contingency Table',     'p-value',
                              'HetOR1',
                              'HomOR1',
                              'HetOR2',
                              'HomOR2',
                              'Score with covariates', 'p-value',
                              'HetOR1',
                              'HomOR1',
                              'HetOR2',
                              'HomOR2',
                              'DF'] )

  # For each locus
  for lname,vars,y,X in models:
    if 0:
      f = csv.writer(file('%s.csv' % lname,'w'))
      f.writerow(vars)
      f.writerows(X.tolist())

    n = X.shape[1]
    k = n-X_null.shape[1]

    if k == 2:
      g1,g2   = X.T[1].T==1,X.T[2].T==1
      g       = -(g1|g2),g1,g2
      gs      = where(X.T[1].T,1,0) + where(X.T[2].T,2,0)
      indices = (1,2,n+1,n+2)
    elif k == 1:
      g1,g2   = X.T[1].T==1,0
      g       = -g1,g1
      gs      = where(X.T[1].T,1,0)
      indices = (1,n+1)
    elif k == 0:
      continue
    else:
      raise ValueError,'Unexpected number of parameters in model (n=%d,k=%d)' % (n,k)

    ys = y==0,y==1,y==2

    sys.stdout.write(lname)
    df = 2*k

    if 1: # not options.skipunadjusted:
      ### Contingency table test ###
      counts = contingency_table(y,gs)

      if 0:
        print
        print counts
        print counts.sum(axis=0)
        print counts.sum(axis=1)
        print

      ct,df_c,ors = contingency_analysis(counts)
      sys.stdout.write('\t%8.5f\t%9.7f' % (ct,1-stats.chi2.cdf(ct,df_c)))

      ors = [ '%7.4f' % orr for orr in ors.flatten() ]
      if len(ors) < 4:
        ors += ['']*(4-len(ors))
      orsstr = '\t'.join(ors)

      sys.stdout.write('\t%s' % orsstr)

    if 1:
      ### SCORE TEST ###
      s1,s2,w11,w12,w22 = build_scores3(ys,X,bg[k])
      st,df_s = do_score_test3(s1,s2,w11,w12,w22,X,g,indices)
      sys.stdout.write('\t%8.5f\t%9.7f' % (st,1-stats.chi2.cdf(st,df_s)))

    if 1:
      try:
        L_model,b_model,W_model = logit3(y,X,initial_beta=bg[k],add_mean=False)
        ors = exp(b_model).take(indices).A.flatten()
      except LinAlgError:
        L_model,b_model,W_model = nan,matrix([]),matrix([[]])
        ors = array([])

      if len(ors) and ors.max() > 1000:
        ors = []
        st = nan
      else:
        ors = [ '%7.4f' % orr for orr in ors ] + ['']*(4-len(ors))

      ors += ['']*(4-len(ors))
      orsstr = '\t'.join(ors)

      sys.stdout.write('\t%s' % orsstr)

    if 0:
      n = 1000
      m = 0
      for y,X in islice(perms,n):
        if k == 2:
          g1,g2   = X.T[1].T==1,X.T[2].T==1
          g       = -(g1|g2),g1,g2
        elif k == 1:
          g1,g2   = X.T[1].T==1,0
          g       = -g1,g1
        stp = do_score_test3(s1,s2,w11,w12,w22,X,g,indices)
        if stp >= st:
          m+= 1
      sys.stdout.write('\t%9.7f' % (float(m)/n))

    if 0:
      ### LR Model with covariates ###
      lr = lr_test3(y,X, [0]+range(1+k,n), b_initial=bg[k], b_initial_null=b_null)
      sys.stdout.write('\t%8.5f\t%9.7f' % (lr,1-stats.chi2.cdf(lr,df)))

      ### Model without covariates ###
      X_nocov = X.take( range(0,k+1), axis=1 )
      lr_nocov = lr_test3(y, X_nocov, [0], b_initial=bng[k], b_initial_null=b_null_nocov)
      sys.stdout.write('\t%8.5f\t%9.7f' % (lr_nocov,1-stats.chi2.cdf(lr_nocov,df)))

    sys.stdout.write('\t%d\n' % df)


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
