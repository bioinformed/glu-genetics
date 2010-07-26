import numpy as np
import scipy as sp
import scipy.sparse, scipy.linalg

import glu.lib.glm._glmnet

from   glu.lib.fileutils import table_reader


class GLMnetResultsBase(object):
  def _error(self):
    n = self.error_code

    # No error
    if n == 0:
      return

    # Fatal errors
    if 0 < n < 7777:
      raise RuntimeError('glmnet: memory allocation error')
    elif n == 7777:
      raise RuntimeError('glmnet: all predictors have zero variance')
    elif 8000 <= n < 9000:
      raise RuntimeError('glmnet: null probability for class %d < 1.0e-5' % (n-8000))
    elif 9000 <= n < 10000:
      raise RuntimeError('glmnet: null probability for class %d > 1.0-1.0e-5' % (n-9000))
    elif n == 10000:
      raise RuntimeError('glmnet: all penalty factors <= 0')
    elif n>0:
      raise RuntimeError('glmnet: unknown fatal error')

    # Non-fatal errors
    # XXX: Should be warnings
    elif -10000<=n<0:
      raise RuntimeError('glmnet: convergence for %d-th lambda value not reached after %d iterations; solutions for larger lambdas returned'
                                   % (-n-1,self.maxit))
    elif n < -10000:
      raise RuntimeError('glmnet: number of nonzero coefficients along the path exceeds %d at %d-th lambda value; solutions for larger lambdas returned'
                                   % (-n-10001,self.pmax))
    elif n < 0:
      raise RuntimeError('glmnet: unknown non-fatal error')


  def indices(self, i=None):
    if i is None:
      i = -1

    return self.variables[:self.compressed_size[i]]-1



class GLMnetGaussianResults(GLMnetResultsBase):
  def __init__(self, x, y, weights, results, nobs, nvars, maxit, pmax):
    lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr = results

    self.nobs             = nobs
    self.nvars            = nvars
    self.solutions        = lmu
    self.intercept        = a0
    self.compressed_coef  = ca
    self.variables        = ia
    self.compressed_size  = nin
    self.r2               = rsq
    self.lambdas          = alm
    self.iters            = nlp
    self.error_code       = jerr
    self.maxit            = maxit
    self.pmax             = pmax
    self.cov              = None

    self._error()

    m                     = self.compressed_size[-1]+1
    self.resids           = y-self.predict(x)
    self.ss               = (self.resids**2).sum()/(nobs-m)
    self.r2_unconstrained = 1-self.resids.var()/y.var()

    self.compute_covariance(x, weights)

  def coef(self, i=None):
    if i is None:
      i = -1

    n    = self.compressed_size[i]
    b    = np.zeros(self.nvars+1)
    b[0] = self.intercept[i]
    b[self.variables[:n]] = self.compressed_coef[:n,i]

    return b

  def predict(self, x, i=None):
    if i is None:
      i = -1

    n  = self.compressed_size[i]
    f  = self.intercept[i]
    f += np.dot(x[:,self.variables[:n]-1],self.compressed_coef[:n,i])

    return f

  def compute_covariance(self, x, weights=None):
    indices = self.indices()
    xx = np.hstack( (np.ones( (self.nobs,1) ), x[:,indices]))

    if weights is not None:
      xx *= np.sqrt(weights[:,np.newaxis])

    # XXX: Linreg may be faster
    u,s,vt = scipy.linalg.svd(xx,full_matrices=False)

    # Truncated inverse squared singular values
    mask = (s<0)
    s2 = s**-2
    s2[mask] = 0
    rank = (s2!=0).sum()

    # Form covariance matrix based on SVD using truncated singular values:
    #   X     = U*S*V'                 U is unitary, so U'*U=U*U'=I;
    #  X'X    = V*S^2*V'               V is unitary, so V'*V=V*V'=I;
    #   I     = V*S^2*V'*(V*S^-2*V')   S is diagonal positive semi-definite
    # (X'X).I = V*S**-2*V'
    self.cov = np.dot(vt.T,(s2*vt.T).T)


class GLMnetMultinomialResults(GLMnetResultsBase):
  def __init__(self, x, y, weights, results, nobs, nvars, categories, maxit, pmax):
    lmu,a0,ca,ia,nin,dev,alm,nlp,jerr = results

    self.nobs               = nobs
    self.nvars              = nvars
    self.categories         = categories
    self.solutions          = lmu
    self.intercept          = a0
    self.compressed_coef    = ca
    self.variables          = ia
    self.compressed_size    = nin
    self.deviance_explained = dev
    self.lambdas            = alm
    self.iters              = nlp
    self.error_code         = jerr
    self.maxit              = maxit
    self.pmax               = pmax
    self.cov                = None

    self._error()

    self.compute_covariance(x, weights)


  def coef(self, cat, i=None, ref=0):
    if i is None:
      i = -1

    n    = self.compressed_size[i]
    p    = self.variables[:n]
    b    = np.zeros(self.nvars+1)
    b[0] = self.intercept[cat,i]
    b[p] = self.compressed_coef[:n,cat,i]

    if ref is not None:
      b[0] -= self.intercept[ref,i]
      b[p] -= self.compressed_coef[:n,ref,i]

    return b

  def predict_logit(self, cat, x, i=None, ref=0):
    if i is None:
      i = -1

    n = self.compressed_size[i]
    p = self.variables[:n]-1
    f = self.intercept[cat,i]

    coef = self.compressed_coef[:n,cat,i]

    if ref is not None:
      f     -= self.intercept[ref,i]
      # Copy needed or else we mutate the originals
      coef  = coef-self.compressed_coef[:n,ref,i]

    f += np.dot(x[:,p],coef)

    return f

  def predict(self, cat, x, i=None, ref=0):
    p = self.predict_logit(cat, x, i=i, ref=ref)
    return 1/(1+np.exp(-p))

  def compute_covariance(self, x, weights=None):
    # XXX: Handle only logistic covariance -- multinomial is a bit hairier
    if len(self.categories) != 2:
      return

    indices = self.indices()
    xx = np.hstack( (np.ones( (self.nobs,1) ), x[:,indices]))

    indices = np.append([0], indices+1)
    beta = self.coef(1)[indices]
    eta1 = np.dot(xx, beta)

    # Form weights
    mu1  = 1/(1+np.exp(-eta1))
    w    = mu1 * (1-mu1)
    xx  *= np.sqrt( (w*weights)[:,np.newaxis] )

    # Compute SVD of weighted design matrix
    # XXX: Linreg may be faster
    u,s,vt = scipy.linalg.svd(xx,full_matrices=False)

    # Truncated inverse squared singular values
    mask = (s<0)
    s2 = s**-2
    s2[mask] = 0
    rank = (s2!=0).sum()

    # Form covariance matrix based on SVD using truncated singular values:
    #   X     = U*S*V'                 U is unitary, so U'*U=U*U'=I;
    #  X'X    = V*S^2*V'               V is unitary, so V'*V=V*V'=I;
    #   I     = V*S^2*V'*(V*S^-2*V')   S is diagonal positive semi-definite
    # (X'X).I = V*S**-2*V'
    self.cov = np.dot(vt.T,(s2*vt.T).T)


def glmnet(x,y,weights=None,family='gaussian',alpha=1.0,nlambda=None,lambda_min=None,
               lambdas=None,standardize=True,thresh=1e-4,dfmax=None,pmax=None,exclude=None,
               penalty_factor=None,maxit=100,hessian_exact=False,update='covariance'):

  if family not in ('gaussian','binomial','multinomial'):
    raise ValueError('Unknown family: %s' % family)

  if update not in ('covariance','naive'):
    raise ValueError('Unknown update: %s' % update)

  if x.ndim != 2:
    raise ValueError('expected x to be two-dimensional array')

  nobs,nvars = x.shape

  if y.shape[0] != nobs:
    raise ValueError('x and y have a different number of observations')

  if weights is None:
    weights = np.ones(nobs)

  if lambda_min is None:
    if nobs < nvars:
      lambda_min = 5e-2
    else:
      lambda_min = 1e-4

  if pmax is None:
    pmax = min((nvars+1)*1.2, nvars)

  if family == 'gaussian':
    if y.ndim != 1:
      raise ValueError('expected y to be one-dimensional array')

  else:
    ka = {'covariance':1, 'naive':2}[update]

    if y.ndim==2 and y.shape[1]==1:
      y = y[:,0]

    if y.ndim==1:
      categories = sorted(set(y.flat))
      y = np.asarray(np.vstack([ y==cat for cat in sorted(categories) ]).T, dtype=float)
    elif y.ndim != 2:
      raise ValueError('expected y to be two-dimensional array')
    else:
      categories = range(y.shape[1])

    if family=='binomial' and y.shape[1] > 2:
      raise ValueError('More than two classes detected; use multinomial family instead')

    nc = y.shape[1]

    valid = weights>0
    if not np.all(valid):
      y    = y[valid]
      x    = x[valid]
      nobs = y.shape[0]

    y = np.multiply(y,weights[:,np.newaxis])

  if exclude:
    jd = np.unique(exclude,dtype=int)
    if jd[0] < 0 or jd[-1] >= nvars:
      raise ValueError('Excluded variables out of range')
  else:
    jd = np.zeros(0,dtype=int)

  if penalty_factor is None:
    vp = np.ones(nvars)
  else:
    vp  = np.asarray(penalty_factor,dtype=float)

  isd = bool(standardize)

  if nlambda is None:
    nlambda = min(100,nvars)

  if lambdas is None:
    if lambda_min >= 1:
      raise ValueError('lambda_min should be less than 1')
    flmin = lambda_min
    # FIXME: Should not need any ulams
    ulam = np.zeros(nlambda)
  else:
    if np.any(lambdas<0):
      raise ValueError('lambdas must be non-negative')
    flmin = 1
    nlamdas = len(lambdas)
    ulam  = np.array(sorted(lambdas,reverse=True),dtype=float)

  nx = int(pmax)

  if dfmax is None:
    ne = nvars+1
  else:
    ne = dfmax

  if family=='gaussian':
    res = _glmnet.elnet(x,y,weights,jd,vp,nx,flmin,ulam,nlam=nlambda,thr=thresh)
    results = GLMnetGaussianResults(x, y, weights, res, nobs, nvars, maxit, nx)

    print 'nvars:',results.nvars
    print 'nobs:',results.nobs
    print 'Solutions:',results.solutions
    print 'Intercepts:',results.intercept
    print 'coefficients:',results.compressed_coef
    print 'variables:',results.variables
    print 'degrees of freedom:',results.compressed_size
    print 'r2_constrained:',results.r2
    print 'r2_unconstrained:',results.r2_unconstrained
    print 'lambdas',results.lambdas
    print 'total passes over the data:',results.iters
    print 'error flag:',results.error_code
    print 'coef',results.coef()
    print 'pred',results.predict(x)
    print 'ss  ',results.ss
    print 'cov ',results.cov
    print 'var ',results.cov.diagonal()
    print 'stde',results.cov.diagonal()**.5

  else:
    res = _glmnet.lognet(x,y,jd,vp,nx,flmin,ulam,kopt=hessian_exact,nlam=nlambda,thr=thresh)
    results = GLMnetMultinomialResults(x, y, weights, res, nobs, nvars, categories, maxit, nx)

    print 'nvars:',results.nvars
    print 'nobs:',results.nobs
    print 'categories:',results.categories
    print 'Solutions:',results.solutions
    print 'Intercepts:',results.intercept
    print 'coefficients:',results.compressed_coef
    print 'variables:',results.variables
    print 'degrees of freedom:',results.compressed_size
    print 'deviance explained:',results.deviance_explained
    print 'lambdas',results.lambdas
    print 'total passes over the data:',results.iters
    print 'error flag:',results.error_code

    for cat in range(len(results.categories)):
      print 'coef%d' % cat,results.coef(cat)
      print 'f%d' % cat, results.predict_logit(cat,x,ref=None)
      print 'l%d' % cat, results.predict(cat,x,ref=None)
      print 'y%d' % cat, y[:,cat]
      print 'dev%d' % cat, ((y[:,cat]-results.predict(cat,x,ref=None))**2).sum()

    print 'cov ',results.cov
    print 'var ',results.cov.diagonal()
    print 'stde',results.cov.diagonal()**.5

  return results


def main():
  x = np.array([[ 0.60, 1.20, 3.90],
                [ 5.00, 4.00, 2.50],
                [ 1.00,-4.00,-5.50],
                [-1.00,-2.00,-6.50],
                [-4.20,-8.40,-4.80]])
  y = np.array([3, 4, -1, -5, -1]).T

  glmnet(x,y)


def main():
  x = table_reader('ramaswamy_x.csv',columns='2-',skip=1)
  x = np.array( [ map(float,row) for row in x ], dtype=float)

  y = table_reader('ramaswamy_y.csv',columns='2-',skip=1)
  y = np.array( [ map(float,row) for row in y ], dtype=float)

  glmnet(x,y,family='multinomial',nlambda=200)


def main3():
  D = []
  indices = [0,1,2]
  for line in file('heart.dat'):
    fields = line.split()
    row = [ fields[i] for i in indices ]
    row = map(float,row)
    D.append(row)

  D=np.asarray(D,dtype=float)

  y,x = D[:,0].astype(int),D[:,1:]

  glmnet(x,y,family='binomial',thresh=1e-6)

  print 'expected beta=[-6.36346994 -1.02410605  0.11904492]'
  print 'expected stde=[ 3.21389766  1.17107845  0.05497905]'


if __name__=='__main__':
  main()
