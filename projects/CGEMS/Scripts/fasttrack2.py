import csv
import sys
import random
import cephes

from   math           import sqrt
from   numpy          import array,matrix,vsplit,zeros,where,hstack,vstack,log,exp,\
                             compress,sometrue,corrcoef
from   numpy.linalg   import LinAlgError
from   itertools      import izip,imap,repeat,islice
from   operator       import itemgetter

from   biozilla.utils import pick,tally,autofile
from   biozilla.xtab  import xtab_list as xtab
from   logit          import mlogit3,stack2x2


def chdtrc(df,xx):
  if xx < 0 or xx>100 or not df:
    return 1.
  return cephes.chdtrc(df,xx)


def chdtri(df,p):
  if p >= 1:
    p = 0.999
  return cephes.chdtri(df,p)


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

      pid = pheno[idx['PID']]
      if pid not in drops:
        yield pheno + ['0',cc,'0']

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


def build_genomap(locus,geno_indices):
  genocounts = tally(pick(locus,geno_indices.itervalues()))
  genocounts.pop('',None)

  def genokey(i):
    g,n = i
    return g[0]==g[1],n
  genocounts = sorted(genocounts.items(), key=genokey)

  #print locus[0],[genocounts[-1]]+genocounts[:-1]
  counts = [(g,i) for i,(g,n) in enumerate(genocounts) ]
  return dict(counts),[ g for g,i in counts ]


def build_genomap2(locus1,locus2,genos1,genos2,geno_indices):
  genos1 = set(genos1)
  genos2 = set(genos2)

  genos = zip(pick(locus1,geno_indices.itervalues()),
              pick(locus2,geno_indices.itervalues()))

  genocounts = tally(g for g in genos if '' not in g)
  i1 = {}
  i2 = {}
  for g1,g2 in genocounts:
    i1.setdefault(g1,len(i1))
    i2.setdefault(g2,len(i2))

  i1 = sorted(i1.items(), key=itemgetter(1))
  i2 = sorted(i2.items(), key=itemgetter(1))

  if 0:
    print
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

  pids = [row[0] for row in phenos]
  covs = [ map(float,row[2:]) for row in phenos ]

  pidset = set(pids)
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

  genos = dict( (locus[0],locus) for locus in genodata)

  for lname1,lname2 in loci:
    locus1 = genos[lname1]
    locus2 = genos[lname2]

    genomap1,genonames1   = build_genomap(locus1,geno_indices)
    genomap2,genonames2   = build_genomap(locus2,geno_indices)
    genomap12,genonames12 = build_genomap2(locus1,locus2,genonames1[:-1],genonames2[:-1],geno_indices)

    if not genomap1 or not genomap2:
      continue

    k1  = len(genomap1)-1
    k2  = len(genomap2)-1
    k12 = len(genomap12)

    #print genonames1[:-1],genonames2[:-1],genonames12
    gxx = []
    def _gen():

      X = []
      y = []
      valid_pids = []

      for pid,stat,cov in izip(pids,stats,covs):
        g1 = locus1[geno_indices[pid]]
        g2 = locus2[geno_indices[pid]]

        if not g1 or not g2:
          continue

        valid_pids.append(pid)

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

      yield y,X

      #### FIXME: Needs to be updated for 2 locus models ####
      raise NotImplementedError,'Needs to be updated for 2 locus models'

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

    vars = (['_mean'] + [ '%s[%s]' % (lname1,g) for g,n in genonames1[:-1] ]
                      + [ '%s[%s]' % (lname2,g) for g,n in genonames2[:-1] ]
                      + header[2:])
    gen = _gen()
    y,X = gen.next()
    yield lname1,lname2,len(genonames1)-1,len(genonames2)-1,len(genomap12),vars,y,X,gen


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


def contingency_table(ys,g):
  n,m = max(ys)+1,max(g)+1
  c = zeros( (n,m) )

  for i in range(n):
    y = ys==i
    for j in range(m):
      c[i,j] = (y & (g==j) ).sum()

  #print
  #print c

  c = c.compress(sometrue(c),axis=1).compress(sometrue(c,axis=1),axis=0)
  return c


def contingency_analysis(ys,g,mincell=5):
  c = contingency_table(ys,g)

  m,n = c.shape

  if m<2 or n<2:
    return None,0

  df = (n-1)*(m-1)

  if c.min() <= mincell:
    import rpy
    import Numeric

    cc = Numeric.array(c.tolist())
    p = rpy.r.fisher_test(cc, conf_int=False, workspace=20000000)['p.value']
    t = chdtri(df,p)
  else:
    m1 = c.sum(axis=0)
    m2 = c.sum(axis=1)
    e  = matrix(m2).T*matrix(m1)/float(sum(m1))
    t  = float(((c-e).A**2/e).sum())

  n,m = c.shape
  return t,max(0,(n-1)*(m-1))


def score_test(ys,g,X,b,indices):
  s1,s2,w11,w12,w22 = build_scores(ys,X,b)
  return do_score_test(s1,s2,w11,w12,w22,X,g,indices)


def build_scores(ys,X,b):
  n = X.shape[1]

  # Form linear predictors and expected values for the current data
  # using the estimates under the null
  eta1,eta2 = X*b[:n],X*b[n:]
  mu0 = 1/(1+exp(eta1)+exp(eta2)).A
  mu1 = exp(eta1).A*mu0
  mu2 = exp(eta2).A*mu0

  # Form the inverse information matrix
  w11 = mu1 * (1-mu1)
  w12 = mu1 *   -mu2
  w22 = mu2 * (1-mu2)

  # Compute score vector
  s1  = ys[1] - mu1
  s2  = ys[2] - mu2

  return s1,s2,w11,w12,w22


def do_score_test(s1,s2,w11,w12,w22,X,g,indices):
  df = len(indices)
  k = df/2

  s = matrix([ where(X.T[i].T, s1, 0).sum() for i in indices[:k] ]
          +  [ where(X.T[i].T, s2, 0).sum() for i in indices[:k] ])

  # Form the inverse information matrix
  W11 = X.T*(w11*X.A)
  W12 = X.T*(w12*X.A)
  W22 = X.T*(w22*X.A)

  try:
    V = stack2x2(W11, W12,
                 W12, W22).I
  except LinAlgError:
    return 0,df

  # Extract the rows and columns that deal with the genotype parameters
  VV = V.take(indices).take(indices,axis=1)

  # Compute test statistic
  t = float((s*VV*s.T)[0,0])

  if t>100:
    t = 0

  return t,df


def lr_test(y,X,indices, b_initial=None,b_initial_null=None):
  X_null = X.take(indices, axis=1)
  df = 2*(X.shape[1]-X_null.shape[1])

  try:
    ### Model with covariates ###
    # Fit the full model
    L1,b_a,W = mlogit3(y,X,initial_beta=b_initial,add_mean=False)
    #print_results(vars,X,L1,b_a,W)

    # Fit the reduced model
    L0,b_n,W = mlogit3(y,X_null,initial_beta=b_initial_null,add_mean=False)

    # Compute the likelihood ratio statistic
    lr = float(2*(L1-L0))
  except LinAlgError:
    lr = 0

  return lr,df


def permutations_needed(p_hat, w, g):
  '''
  Estimate the number of permutations required to estimate an empirical
  p-value of a test statistic using a Monte Carlo permutation procedure with N
  replicate permutations.  We choose N such that the estimated empirical
  p-value, p_hat, is within a proportion w (the width parameter) of its true
  p-value, p, with predetermined confidence probability g. That is, we want
  the standard deviation sigma_p_hat of p_hat to be proportional to p_hat.
  This permutation process can be viewed as a set of N independent Bernoulli
  trials each with success probability p.

  Using a Normal approximation for the distribution of p_hat, we obtain

       1 - p_hat
  N = -----------  Phi^-2[ (g+1)/2 ]
      w^2 * p_hat

  where Phi^-2 is the squared inverse of the standard normal cumulative
  distribution function.

  For example, to estimate an empirical p-value that is within 20% of its
  true value with 95% confidence (w=.2, g=.95), then N is approximately
  N ~= 100*(1-p_hat)/p_hat.  For suitably small values of p_hat, this
  further simplifies to N ~= 100/p_hat.
  '''
  import cephes


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


def main():
  header,phenos = load_phenos(sys.argv[1])
  header = map(str.upper,header)
  header,phenos = filter_pheno_records(header,phenos)

  if 1:
    header,phenos = expand_phenos(header,phenos)
    keep = ['PID','STAT',
            'CENTER01','CENTER02','CENTER04','CENTER05','CENTER06',
            'CENTER08','CENTER09','CENTER10','PREV',
            'AGELEVEL1','AGELEVEL2','AGELEVEL3']
  else:
    header,phenos = expand_phenos_simple(header,phenos)
    keep = ['PID','STAT',
            'CENTER01','CENTER02','CENTER04','CENTER05','CENTER06',
            'CENTER08','CENTER09','CENTER10',
            'AGELEVEL1','AGELEVEL2','AGELEVEL3']
    #keep = ['PID','STAT']

  agemap = dict( (str(i),'AGELEVEL%s' % i) for i in range(4) )
  header,phenos = expand_nominal(header,phenos,'AGELEVEL',agemap)

  centermap = dict( (str(i),'CENTER%02d' % i) for i in range(12) )
  header,phenos = expand_nominal(header,phenos,'CENTER',centermap)

  yrmap = dict( (str(i),'YRIND%d' % i) for i in range(9) )
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  yrmap = dict( (str(i),'INCIDENT') for i in range(9) )
  yrmap['0'] = 'PREV'
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  if 0:
    pids = set(csv.reader(autofile(sys.argv[2]),dialect='excel-tab').next()[1:])
    phenos = [ pheno for pheno in phenos if pheno[0] in pids ]
    f = csv.writer(file('phenos_fasttrack.csv','w'))
    f.writerow(header)
    f.writerows(phenos)
    return

  header,phenos = filter_pheno_columns(header,phenos,keep)
  phenos = list(phenos)

  genodata = list(csv.reader(autofile(sys.argv[2]),dialect='excel-tab'))
  reflocus = sys.argv[2].split('/')[-1].split('.')[0]
  loci = ( (reflocus,locus[0]) for locus in islice(genodata,1,None) if locus[0]!=reflocus )

  #loci = [('rs1447295','rs7837688')]
  models = build_two_locus_models(loci,genodata,header,phenos)

  if 1:
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
                     b_null[n],[[0]]*k,b_null[n+1:]) ) for k in range(1,5) ]

    bng = [None]
    bng+= [ vstack( (b_null_nocov[0],[[0]]*k,
                     b_null_nocov[1],[[0]]*k) ) for k in range(1,3) ]

  print '\t'.join( ['Ref Locus','Comp Locus',
                    'Contingency Ref', 'df',
                    'Contingency Comp', 'df',
                    'Stratified Contingency Hybrid', 'df',
                    'Two locus LR',
                    'Two locus Score', 'df'] )

  # For each locus
  for lname1,lname2,k1,k2,k12,vars,y,X,perms in models:
    #if lname2 != 'rs6991990':
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

    sys.stdout.write('%s\t%s' % (lname1,lname2))
    df = 2*k

    if 1:
      ### Contingency table test ###
      gg1 = array((g1[1] + 2*g1[2]).flat)
      gg2 = array((g2[1] + 2*g2[2]).flat)
      yy  = array((ys[1] + 2*ys[2]).flat)

      ct,df = contingency_analysis(yy,gg1)
      #sys.stdout.write('\t%8.5f\t%9.7f' % (ct,chdtrc(df,ct)))
      sys.stdout.write('\t%9.7f\t%d' % (chdtrc(df,ct),df))

      ct,df = contingency_analysis(yy,gg2)
      #sys.stdout.write('\t%8.5f\t%9.7f' % (ct,chdtrc(df,ct)))
      sys.stdout.write('\t%9.7f\t%d' % (chdtrc(df,ct),df))

      ct = 0
      df = 0
      for i in range(k1+1):
        t,d = contingency_analysis( compress(gg1==i,yy), compress(gg1==i,gg2) )
        #print t,d,chdtrc(d,t)
        if t is not None:
          ct += t
          df += d

      #print ct
      #sys.stdout.write('\t%8.5f\t%9.7f' % (ct,chdtrc(df,ct)))
      sys.stdout.write('\t%9.7f\t%d' % (chdtrc(df,ct),df))

    # FISHER EXACT TEST
    if 0:
      ct = 0
      df = 0
      for i in range(k1+1):
        pp = contingency_analysis_exact( compress(gg1==i,yy), compress(gg1==i,gg2) )
        if pp is not None:
          ct += log(pp)
          df += 2
      ct = -2*float(ct)
      sys.stdout.write('\t%9.7f' % chdtrc(df,ct))

    if 1:
      ### LR Model with covariates ###
      indices = range(k1+1)+range(k1+k2+k12+1,n)
      lr,df = lr_test(y,X,indices)
      #df = 2*(k2+k12)
      #sys.stdout.write('\t%8.5f\t%9.7f' % (lr,chdtrc(df,lr)))
      sys.stdout.write('\t%9.7f' % chdtrc(df,lr))

    if 0:
      ### Model without covariates ###
      X_nocov = X.take( range(0,k+1), axis=1 )
      lr_nocov,df = lr_test(y, X_nocov, [0], b_initial=bng[k], b_initial_null=b_null_nocov)
      sys.stdout.write('\t%8.5f\t%9.7f' % (lr_nocov,chdtrc(df,lr_nocov)))

    if 1:
      ### SCORE TEST ###

      # Obtain estimates of covariate effects under the null
      indices2 = (range(  k1+1,  k1+k2+k12+1)
               +  range(n+k1+1,n+k1+k2+k12+1))
      indices = [i for i in range(X.shape[1]) if i not in indices2]

      X_null = X.take(indices, axis=1)
      L,b_null,W = mlogit3(y,X_null,add_mean=False)

      b_null = list(b_null.flat)
      b_null[k1+1:k1+1]     = [0]*(k2+k12)
      b_null[n+k1+1:n+k1+1] = [0]*(k2+k12)
      b_null = matrix(b_null).T
      #print b_null.T

      s1,s2,w11,w12,w22 = build_scores(ys,X,b_null)
      st,df = do_score_test(s1,s2,w11,w12,w22,X,g2,indices2)
      #sys.stdout.write('\t%8.5f\t%9.7f' % (st,chdtrc(df,st)))
      sys.stdout.write('\t%9.7f' % chdtrc(df,st))

    # FIXME: Permutation is not right for loci in LD
    if 0:
      n = 1000
      m = 0
      for y,X in islice(perms,n):
        g1,indices1 = geno_vectors(X,k1,1)
        g2,indices2 = geno_vectors(X,k2,1+k1)
        stp,df = do_score_test(s1,s2,w11,w12,w22,X,g,indices)
        if stp >= st:
          m+= 1
      sys.stdout.write('\t%9.7f' % (float(m)/n))

    sys.stdout.write('\t%d' % df)
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
