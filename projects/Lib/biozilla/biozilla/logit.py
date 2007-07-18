from math      import pi
from itertools import izip

# Import Python numerical libraries (from http://numeric.scipy.org/)
from numpy import array, matrix, asmatrix, zeros, zeros_like, ones, where, exp, log, vsplit, unique, bmat

from utils import tally

CONV = 1e-8

def grand_mean(X):
  return bmat([ ones((len(X),1)), X ])


def linreg(y, X, add_mean=False):
  # Add type specific means, if requested
  if add_mean:
    X = grand_mean(X)

  # n=observations and m=covariates
  n,m = X.shape

  W  = (X.T*X).I
  b  = W*X.T*y
  ss = ((y - X*b).A**2).sum()/(n-m)
  return b,W,ss


def logit(y, X, initial_beta=None, add_mean=True, max_iterations=50):
  # Add type specific means, if requested
  if add_mean:
    X = grand_mean(X)

  # n=observations and m=covariates
  n,m = X.shape

  # Initial estimate for mu1_0 and mu2_0
  m1 = float(y.sum())/n

  # Initial estimates of beta (b1,b2)
  if initial_beta is not None:
    b = initial_beta
  else:
    b = asmatrix(zeros((m,1), dtype=float))
    b[0,0] = log(m1/(1-m1))

  # Initialize likelihood to -infinity
  L = -1e308*1e308

  # Iterate until convergence
  for i in range(max_iterations):
    eta1 = X*b

    # Compute expected values and linear predictors
    mu0  = 1/(1+exp(eta1)).A
    mu1  = exp(eta1).A*mu0
    eta0 = log(mu0)

    # Compute likelihood
    newL = eta0.sum() + y.choose( (0,eta1) ).sum()

    d,L = newL-L,newL

    # Stop when likelihood converges
    if d < CONV:
      break

    # Form weights and quadrants of the expected information matrix
    w = mu1 * (1-mu1)

    # Build information matrix and compute the inverse
    W = (X.T*(w*X.A)).I

    # Update regression estimates
    b = W*X.T*(w*eta1.A + y - mu1)

  else:
    raise RuntimeError,'logit estimator failed to converge'

  # Return log-likelihood, regression estimates, and covariance matrix
  return L,b,W


class GLogit(object):

  def __init__(self, y, X, ref=None, vars=None, add_mean=True, max_iterations=50):
    y = asmatrix(y,dtype=int)
    X = asmatrix(X,dtype=float)

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    # Add type specific means, if requested
    if add_mean:
      X = grand_mean(X)

    categories = unique(y.A.ravel())

    if ref is None:
      ref = categories[0]
    else:
      if ref not in categories:
        raise ValueError, 'Reference class for dependent variable not observed in data'
      categories = [ref] + [ r for r in categories if r!=ref ]

    # Recode y with ordinal categories from 0..k-1, where 0 is reference
    # FIXME: This should be a scalar loop
    y_ord = zeros_like(y)
    for i,v in enumerate(categories):
      y_ord += where(y==v, i, 0)

    self.y     = y
    self.y_ord = y_ord
    self.X     = X
    self.vars  = vars
    self.categories = categories

    self.L    = None
    self.beta = None
    self.W    = None


  def fit(self, initial_beta=None, max_iterations=50):
    # Form the dependent variable vector for classes y==1 and y==2
    y = self.y_ord
    X = self.X
    k = len(self.categories)-1

    # n=observations and m=covariates
    n,m = X.shape

    # Initial estimate for mu1_0 and mu2_0
    mu = array( [(y==i).sum() for i in self.categories[1:] ], dtype=float)/n

    # Initial estimates of beta (b1,b2)
    if initial_beta is None:
      b = []
      for i in range(k):
        bi = asmatrix(zeros(m), dtype=float).T
        bi[0,0] = log(mu[i]/(1-sum(mu)))
        b.append(bi)
    else:
      b = vsplit(initial_beta,k)

    # Initialize log likelihood to -infinity
    L = -1e308*1e308

    # Iterate until convergence
    for iteration in range(max_iterations):
      eta = [ X*bi for bi in b ]

      # Compute expected values and linear predictors
      mu0  = 1/(1+sum(exp(etai) for etai in eta)).A
      eta0 = log(mu0)
      mu   = [ exp(eta[i]).A*mu0 for i in range(k) ]

      # Compute likelihood
      newL = eta0.sum() + self.y_ord.choose( [0]+eta ).sum()

      d,L = newL-L,newL

      # Stop when likelihood converges
      if d < CONV:
        break

      # Form weights and blocks of the expected information matrix
      w = self.weights(mu)

      # Build information matrix and compute the inverse
      W = self.information_matrix(w).I

      # Update regression estimates
      beta = self.update_estimates(w,W,mu,eta)

      b = vsplit(beta,k)

    else:
      raise RuntimeError,'glogit estimator failed to converge'

    # Return log-likelihood, regression estimates, and covariance matrix
    self.L    = L
    self.beta = beta
    self.W    = W

    return L,beta,W


  def weights(self,mu):
    k = len(mu)
    w = [ [None]*k for i in range(k) ]
    for i in range(0,k):
      for j in range(i,k):
        if i == j:
          w[i][i] = mu[i]*(1-mu[i])
        else:
          w[i][j] = w[j][i] = -mu[i]*mu[j]
    return w


  def information_matrix(self,w):
    X = self.X
    k = len(w)
    W = [ [None]*k for i in range(k) ]
    for i in range(0,k):
      for j in range(i,k):
        W[i][j] = W[j][i] = X.T*(w[i][j]*X.A)
    return bmat(W)


  def update_estimates(self,w,W,mu,eta):
    X = self.X
    k = len(mu)
    v = []
    for i in range(k):
      t = (self.y_ord==i+1) - mu[i]
      for j in range(k):
        t += w[i][j]*eta[j].A
      v += [[X.T*t]]
    return W * bmat(v)

  def score_test(self,beta=None,initial_beta=None,indices=None):
    return GLogitScoreTest(self,beta=beta,initial_beta=initial_beta,indices=indices)

  def wald_test(self,indices=None):
    return GLogitWaldTest(self,indices=indices)


class GLogitScoreTest(object):
  def __init__(self,model,beta=None,initial_beta=None,indices=None):
    y = model.y
    X = model.X
    k = len(model.categories)-1
    n = X.shape[1]

    if indices is None:
      indices = []
      for i in range(k):
        s = i*n
        indices += range(s+1,s+n)

    if beta is None:
      # Form null design matrix by finding the rows of the design matrix
      # that are to be scored and exclusing those from the null model
      design_indices = tally(i%n for i in indices)
      assert design_indices
      assert all(m == k for m in design_indices.itervalues())

      null_indices = [ i for i in range(n) if i not in design_indices ]
      X_null = X[:,null_indices]

      # Fit null model
      L,beta_null,W = GLogit(y,X_null,add_mean=False).fit(initial_beta=initial_beta)

      # Augment null beta with zeros for all parameters to be tested
      null_indices = (i for i in range(n*k) if i not in indices)
      beta = asmatrix(zeros((X.shape[1]*k,1), dtype=float))
      for b,i in izip(beta_null.A.ravel(),null_indices):
        beta[i,0] = b

    self.model   = model
    self.beta    = beta
    scores,w,W   = self.build_scores()
    self.indices = indices
    self.scores  = scores
    self.w = w
    self.W = W


  def build_scores(self):
    y = self.model.y_ord
    X = self.model.X
    k = len(self.model.categories)-1

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    n = X.shape[1]

    # Form linear predictors and expected values for the current data
    # using the estimates under the null
    b = vsplit(self.beta,k)

    eta = [ X*bi for bi in b ]

    # Compute expected values and linear predictors
    mu0  = 1/(1+sum(exp(etai) for etai in eta)).A
    eta0 = log(mu0)
    mu   = [ exp(eta[i]).A*mu0 for i in range(k) ]

    # Form weights and blocks of the expected information matrix
    w = self.model.weights(mu)

    # Build information matrix and compute the inverse
    W = self.model.information_matrix(w).I

    # Compute scores
    # FIXME: Can scores and mu be further vectorized?
    scores = bmat([ (y==i+1) - mu[i] for i in range(k) ])

    return scores,w,W


  def test(self):
    scores  = self.scores
    indices = self.indices
    W       = self.W
    m       = len(scores)
    df      = len(indices)
    X       = self.model.X
    n       = X.shape[1]

    # Compute score vectors with the full design matrix
    indices = array(indices, dtype=int)
    s       = asmatrix(( X[:,indices%n].A*scores[:,indices//n].A ).sum(axis=0))

    # Extract the covariance matrix of the score parameters
    V = W.take(indices,axis=0).take(indices,axis=1)

    # Compute score test statistic
    x2 = float((s*V*s.T)[0,0])

    return x2,df


class GLogitWaldTest(object):
  def __init__(self,model,indices=None):
    if indices is None:
      n = X.shape[1]
      k = len(model.categories)-1
      indices = []
      for i in range(k):
        s = i*n
        indices += range(s+1,s+n)

    self.model   = model
    self.indices = indices

  def test(self):
    indices = self.indices
    beta    = self.model.beta
    W       = self.model.W
    df      = len(indices)

    b  = beta[indices,:]
    w  = W[indices,:][:,indices]
    x2 = float((b.T*w.I*b)[0,0])

    return x2,df


class Linear(object):

  def __init__(self, y, X, vars=None, add_mean=True):
    y = asmatrix(y,dtype=float)
    X = asmatrix(X,dtype=float)

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    # Add type specific means, if requested
    if add_mean:
      X = grand_mean(X)

    self.y     = y
    self.X     = X
    self.vars  = vars

    self.L    = None
    self.beta = None
    self.W    = None
    self.ss   = None

  def fit(self):
    y = self.y
    X = self.X

    # n=observations and m=covariates
    n,m = X.shape

    W      = (X.T*X).I
    beta   = W*X.T*y
    e      = y - X*beta
    L      = -(n/2.0)*(1 + log(2*pi) - log(e.T*e/n))
    ss     = (e.A**2).sum()/(n-m)

    # Return log-likelihood, regression estimates, and covariance matrix
    self.L    = L
    self.beta = beta
    self.W    = W
    self.ss   = ss

    return L,beta,W,ss

  def score_test(self,beta=None,indices=None):
    return LinearScoreTest(self,beta=beta,indices=indices)

  def wald_test(self,indices=None):
    return LinearWaldTest(self,indices=indices)


class LinearScoreTest(object):
  def __init__(self,model,beta=None,indices=None):
    y = model.y
    X = model.X
    n = X.shape[1]

    if indices is None:
      indices = range(1,n)

    if beta is None:
      design_indices = set(indices)
      null_indices = [ i for i in range(n) if i not in design_indices ]
      X_null = X[:,null_indices]

      # Fit null model
      null = Linear(y,X_null,add_mean=False)
      null.fit()

      # Augment null beta with zeros for all parameters to be tested
      beta = asmatrix(zeros((n,1), dtype=float))
      for b,i in izip(null.beta.A.ravel(),null_indices):
        beta[i,0] = b

      self.ss = null.ss
    else:
      # FIXME: residual variance should always be estimated from the null model
      self.ss = model.ss

    self.model   = model
    self.beta    = beta
    self.indices = indices
    self.scores  = model.y - X*beta


  def test(self):
    scores  = self.scores
    indices = self.indices
    W       = self.model.W
    df      = len(indices)
    X       = self.model.X
    ss      = self.ss

    # Compute score vectors with the full design matrix
    indices = array(indices, dtype=int)
    s = asmatrix(( X[:,indices].A*scores.A ).sum(axis=0))

    # Extract the covariance matrix of the score parameters
    V = W[indices,:][:,indices]/ss

    # Compute score test statistic
    x2 = float((s*V*s.T)[0,0])

    return x2,df


class LinearWaldTest(object):
  def __init__(self,model,indices=None):
    if indices is None:
      n = X.shape[1]
      indices = range(1,n)

    self.model   = model
    self.indices = indices

  def test(self):
    indices = self.indices
    beta    = self.model.beta
    W       = self.model.W
    ss      = self.model.ss
    df      = len(indices)

    b  = beta[indices,:]
    w  = W[indices,:][:,indices]*ss
    x2 = float((b.T*w.I*b)[0,0])

    return x2,df


def main():
  # Load and parse test data
  dfile = file('test/cancer.dat')
  D = []
  indices = [4,5,3,6]
  for line in dfile:
    fields = line.split()
    row = [ fields[i] for i in indices ]
    if '.' in row:
      continue
    row = map(float,row)
    D.append(row)

  D=asmatrix(D,dtype=float)
  #D = asmatrix(bmat( [[D]]*50 )) # Inflate data set to check scalability

  # Extract dependent variable (first column) and design matrix (following columns)
  y,X = D[:,0],D[:,1:]

  # Compute the polytomous logit model for 3 type classes
  if 1:
    L,b,W = GLogit(y,X,ref=0).fit()
  elif 1:
    L,b,W = logit3(y,X)
  else:
    y = y==2
    L,b,W = logit(y,X)

  # Compute standard errors
  stde = W.diagonal().A**.5

  # Compute Z-statistic
  z = b.T.A/stde

  # Show results
  print 'obs =',len(D)
  print 'logL =',L
  print '-2logL =',-2*L
  print 'beta =',b.T
  print 'ss =',stde
  print 'Z =',z
  print 'X2 =',z**2
  print 'score test=',GLogit(y,X).score_test().test()


if __name__ == '__main__':
  main()
