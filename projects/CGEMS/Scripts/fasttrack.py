import csv
import sys
import random

from   numpy          import array,matrix,exp,vsplit,zeros,where,hstack,vstack,log,exp,nan
from   itertools      import izip,imap,repeat,islice

from   biozilla.utils import pick,tally,autofile
from   biozilla.logit import mlogit3

from   association    import *


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


def make_index(header):
  return dict( (f,i) for i,f in enumerate(header) )


def filter_pheno_columns(header, phenos, keep):
  indices = [ header.index(k) for k in keep ]

  def _filter_phenos():
    for pheno in phenos:
      yield [ pheno[i] for i in indices ]

  return keep,_filter_phenos()


def expand_phenos(header, phenos):
  new_header = header + ['YR','STAT','PREVALENT']
  idx = make_index(new_header)

  ccmap = {'1':'0','2':'1','3':'2','4':'2'}
  drops = set(r[0] for r in csv.reader(file('drop_samples'),dialect='excel-tab'))

  def _expand_phenos():
    for pheno in phenos:
      cc    = ccmap[pheno[idx['CC_OVERALL']]]
      yrs   = pheno[idx['YR1']],pheno[idx['YR2']],pheno[idx['YR3']]
      stats = pheno[idx['STAT1']],pheno[idx['STAT2']],pheno[idx['STAT3']]

      pid = pheno[idx['PID']]

      if pid in drops:
        continue

      for yr,stat in zip(yrs,stats):
        prev = '1' if yr=='0' else '0'
        if stat == '1':
          yield pheno + [yr,'0',prev]
        elif stat == '2':
          assert cc in ('1','2')
          yield pheno + [yr,cc,prev]
          break

  return new_header,_expand_phenos()


def expand_phenos_simple(header, phenos):
  new_header = header + ['YR','STAT','PREVALENT']
  idx = make_index(new_header)

  ccmap = {'1':'0','2':'1','3':'2','4':'2'}

  drops = set(r[0] for r in csv.reader(file('drop_samples'),dialect='excel-tab'))

  def _expand_phenos():
    for pheno in phenos:
      cc    = ccmap[pheno[idx['CC_OVERALL']]]
      yrs   = pheno[idx['YR1']],pheno[idx['YR2']],pheno[idx['YR3']]
      stats = pheno[idx['STAT1']],pheno[idx['STAT2']],pheno[idx['STAT3']]

      prev = 0
      for yr,stat in zip(yrs,stats):
        if yr == 0:
          prev = '1'
          break

      pid = pheno[idx['PID']]
      if pid not in drops:
        yield pheno + ['0',cc,prev]

  return new_header,_expand_phenos()


def filter_pheno_records(header,phenos):
  def _filter_pheno_records():
    idx = make_index(header)

    for pheno in phenos:
      consent = pheno[idx['HAS_CONSENT']]
      samples = pheno[idx['GENOTYPED_SAMPLE_IDS']]
      center  = pheno[idx['CENTER']]

      if consent == '1' and samples and center != '3':
        yield pheno

  return header,_filter_pheno_records()


def expand_nominal(header,phenos,field,valuemap):
  values = sorted(set(valuemap.values()))
  header = header + values
  pheno_idx = make_index(header)
  value_idx = make_index(values)

  def _expand_nominal():
    for pheno in phenos:
      vals = ['0']*len(values)
      value_field = valuemap[pheno[pheno_idx[field]]]
      vals[ value_idx[value_field] ] = '1'
      yield pheno + vals

  return header,_expand_nominal()


def compute_maf(counts):
  hom1 = hom2 = het = 0

  for g,n in counts.iteritems():
    if g[0] != g[1]:
      het = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)

  maf = float(2*hom1 + het)/(hom1+het+hom2)/2
  return maf


def build_locus_models(loci,header,phenos,mingenos=15):
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

  perms = {}
  for row in phenos:
    pid = row[0]

    if pid not in geno_indices:
      continue

    try:
      j = row[2:].index('1')
    except ValueError:
      j = 0
    perms.setdefault(j,[]).append(pid)


  for locus in loci:
    lname = locus[0]

    genocounts = tally(pick(locus,geno_indices.itervalues()))

    genocounts.pop('',None)
    genocounts.pop('  ',None)

    maf = compute_maf(genocounts)

    #print lname,maf

    if maf < 0.01:
      continue

    poolgenos  = [ g for g,n in genocounts.iteritems() if n<mingenos ]
    genocounts = sorted( ( (n,g) for g,n in genocounts.iteritems() if n>=mingenos ))

    if len(genocounts) < 2:
      continue

    homs = [ g for (n,g) in genocounts if g[0]==g[1]]
    hets = [ g for (n,g) in genocounts if g[0]!=g[1]]

    if len(homs) + len(hets) < 2:
      continue

    testgenos = homs[:1] + hets + homs[1:]

    genomap = dict( (g,i) for i,g in enumerate(testgenos) )

    for g in poolgenos:
      genomap[g] = 1

    def _gen():
      k = max(genomap.itervalues())

      X = []
      y = []
      valid_pids = []

      for pid,stat,cov in izip(pids,stats,covs):
        g = locus[geno_indices[pid]]

        if not g.strip():
          continue

        if g not in genomap:
          continue

        valid_pids.append(pid)

        genos = [0]*k

        i = genomap[g]

        if i!=k:
          genos[i] = 1

        y.append(stat)
        X.append( [1] + genos + cov )

      y = matrix(y, dtype=int)
      X = matrix(X, dtype=float)

      yield y,X

      permsets = [ [p for p in ps if locus[geno_indices[p]]] for ps in perms.itervalues() ]

      X = array(X)
      z,o = array(0.0),array(1.0)

      while 1:
        pidmap = {}
        for permset in permsets:
          permutedset = permset[:]
          random.shuffle(permutedset)
          pidmap.update( izip(permset,permutedset) )

        for i,pid in enumerate(valid_pids):
          g = locus[geno_indices[pid]]
          if not g:
            continue
          pid = pidmap[pid]
          g = locus[geno_indices[pid]]

          j = genomap[g]
          X[i,1:k+1] = z
          if j!=k:
            X[i,j+1] = o

        yield y,matrix(X)

    vars = ['_mean']+ [ '%s[%s]' % (lname,g) for g,n in genocounts[:-1] ] + header[2:]
    gen = _gen()
    y,X = gen.next()
    yield lname,vars,y,X,gen


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


def main():
  header,phenos = load_phenos(sys.argv[1])
  header = map(str.upper,header)
  header,phenos = filter_pheno_records(header,phenos)

  if 1:
    header,phenos = expand_phenos(header,phenos)
  else:
    header,phenos = expand_phenos_simple(header,phenos)

  agemap = dict( (str(i),'AGELEVEL%s' % i) for i in range(4) )
  header,phenos = expand_nominal(header,phenos,'AGELEVEL',agemap)

  centermap = dict( (str(i),'CENTER%02d' % i) for i in range(12) )
  header,phenos = expand_nominal(header,phenos,'CENTER',centermap)

  yrmap = dict( (str(i),'YRIND%d' % i) for i in range(9) )
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  yrmap = dict( (str(i),'INCIDENT') for i in range(9) )
  yrmap['0'] = 'PREV'
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  keep = ['PID','STAT',
          'CENTER01','CENTER02','CENTER04','CENTER05','CENTER06',
          'CENTER08','CENTER09','CENTER10', 'PREV',
          'AGELEVEL1','AGELEVEL2','AGELEVEL3']

  #keep = ['PID','STAT']
  header,phenos = filter_pheno_columns(header,phenos,keep)
  phenos = list(phenos)

  if 0:
    pids = set(csv.reader(autofile(sys.argv[2]),dialect='excel-tab').next()[1:])
    phenos = [ pheno for pheno in phenos if pheno[0] in pids ]
    f = csv.writer(file('phenos.csv','w'))
    f.writerow(header)
    f.writerows(phenos)
    return

  loci   = csv.reader(autofile(sys.argv[2]),dialect='excel-tab')
  models = build_locus_models(loci,header,phenos)

  y      = matrix( [ [int(row[1])]         for row in phenos], dtype=int   )
  X_null = matrix( [[1]+map(float,row[2:]) for row in phenos], dtype=float )

  # Obtain estimates of covariate effects under the null
  L,b_null,W = mlogit3(y,X_null,add_mean=False)

  #vars = ['_mean']+header[2:]
  #print_results(vars,X_null,L,b_null,W)

  # Obtain estimates of mean with no covariates under the null
  L,b_null_nocov,W = mlogit3(y,X_null.take([0], axis=1),add_mean=False)

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
  for lname,vars,y,X,perms in models:
    if 0:
      f = csv.writer(file('%s.csv' % lname,'w'))
      f.writerow(vars)
      f.writerows(X.tolist())

    n = X.shape[1]
    k = n-X_null.shape[1]

    if k == 2:
      g1,g2   = X.T[1].T==1,X.T[2].T==1
      g       = -(g1|g2),g1,g2
      indices = (1,2,n+1,n+2)
    elif k == 1:
      g1,g2   = X.T[1].T==1,0
      g       = -g1,g1
      indices = (1,n+1)
    elif k == 0:
      continue
    else:
      raise ValueError,'Unexpected number of parameters in model (n=%d,k=%d)' % (n,k)

    ys = y==0,y==1,y==2

    sys.stdout.write(lname)
    df = 2*k

    if 1:
      ### Contingency table test ###
      ors,ct,df_c = contingency_analysis(ys,g)
      sys.stdout.write('\t%8.5f\t%9.7f' % (ct,chdtrc(df_c,ct)))

      ors = [ '%7.4f' % orr for orr in ors.flatten() ]
      if len(ors) < 4:
        ors += ['']*(4-len(ors))
      orsstr = '\t'.join(ors)

      sys.stdout.write('\t%s' % orsstr)

      ### SCORE TEST ###
      s1,s2,w11,w12,w22 = build_scores(ys,X,bg[k],indices)
      st = do_score_test(s1,s2,w11,w12,w22,X,g,indices)


      L_model,b_model,W_model = mlogit3(y,X,initial_beta=bg[k],add_mean=False)

      ors = exp(b_model).take(indices).A.flatten()

      if ors.max() > 1000:
        ors = []
        st = nan
      else:
        ors = [ '%7.4f' % orr for orr in ors ] + ['']*(4-len(ors))

      ors += ['']*(4-len(ors))
      orsstr = '\t'.join(ors)

      sys.stdout.write('\t%8.5f\t%9.7f' % (st,chdtrc(df,st)))
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
        stp = do_score_test(s1,s2,w11,w12,w22,X,g,indices)
        if stp >= st:
          m+= 1
      sys.stdout.write('\t%9.7f' % (float(m)/n))

    if 0:
      ### LR Model with covariates ###
      lr = lr_test(y,X, [0]+range(1+k,n), b_initial=bg[k], b_initial_null=b_null)
      sys.stdout.write('\t%8.5f\t%9.7f' % (lr,chdtrc(df,lr)))

      ### Model without covariates ###
      X_nocov = X.take( range(0,k+1), axis=1 )
      lr_nocov = lr_test(y, X_nocov, [0], b_initial=bng[k], b_initial_null=b_null_nocov)
      sys.stdout.write('\t%8.5f\t%9.7f' % (lr_nocov,chdtrc(df,lr_nocov)))


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
