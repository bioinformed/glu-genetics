# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Generalized linear models'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

from   math         import pi
from   itertools    import izip

import numpy as np

import scipy.stats

from   numpy        import dot
from   scipy.linalg import svd,cholesky,norm,inv,schur,rsf2csf,LinAlgError

CONV = 1e-8
COND = 1e-8
EPS  = 1e-5


def bmat2d(data):
  '''
  Construct a block matrix from a 2d list of equal sized components

  >>> x = np.array([[4.,1.,1.],
  ...               [2.,4.,1.],
  ...               [0.,1.,4.]])
  >>> y = np.array([[1.,0.,0.],
  ...               [0.,1.,0.],
  ...               [0.,0.,1.]])
  >>> bmat2d([[x,y],[y,x]])
  array([[ 4.,  1.,  1.,  1.,  0.,  0.],
         [ 2.,  4.,  1.,  0.,  1.,  0.],
         [ 0.,  1.,  4.,  0.,  0.,  1.],
         [ 1.,  0.,  0.,  4.,  1.,  1.],
         [ 0.,  1.,  0.,  2.,  4.,  1.],
         [ 0.,  0.,  1.,  0.,  1.,  4.]])
  '''
  #n = [ d.shape[0]    for d in data[0] ]
  #m = [ d[0].shape[1] for d in data    ]

  n,m = np.atleast_2d(data[0][0]).shape
  if any( d.shape != (n,m) for row in data for d in row ):
    raise ValueError('blocks must all be of same size')

  b = np.zeros( (n*len(data[0]),m*len(data)), dtype=data[0][0].dtype )

  for i in xrange(len(data[0])):
    for j in xrange(len(data)):
      b[i*n:(i+1)*n,:][:,j*m:(j+1)*m] = data[i][j]

  return np.asarray(b)


# Derived from SciPy's linalg.lstsq routine, though heavily modified and
# extended.
def linear_least_squares(a, b, weights=None, sqrtweights=False, cond=None):
  '''linear_least_squares(a, b, cond=None) -> x,ss,cov,resids,rank,s,vt

  Return least-squares solution of a * x = b.

  Inputs:
    a    -- An M x N matrix.
    b    -- An M x nrhs matrix or M vector.
    cond -- Used to determine effective rank of a.

  Outputs:
    x      -- The solution (N x nrhs matrix) to the minimization problem:
              2-norm(| b - a * x |) -> min
    resids -- The residual sum-of-squares for the solution matrix x
    rank   -- The effective rank of a.
    s      -- Singular values of a in decreasing order. The condition number
              of a is abs(s[0]/s[-1]).

  Full-rank Example
  -----------------

  >>> X = np.array([[ 0.60, 1.20, 3.90],
  ...               [ 5.00, 4.00, 2.50],
  ...               [ 1.00,-4.00,-5.50],
  ...               [-1.00,-2.00,-6.50],
  ...               [-4.20,-8.40,-4.80]])
  >>> y = np.array([3, 4, -1, -5, -1]).T

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,cond=5e-4)

  Results:

  >>> beta.T
  array([ 0.95333333, -0.84333333,  0.90666667])
  >>> ss
  0.16999999999999976
  >>> cov
  array([[ 0.06222222, -0.04222222,  0.01333333],
         [-0.04222222,  0.05444444, -0.02888889],
         [ 0.01333333, -0.02888889,  0.02666667]])
  >>> resids
  0.33999999999999952

  Rank, singular Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 15.,   6.,   3.])
  >>> vt
  array([[ 0.33333333,  0.66666667,  0.66666667],
         [ 0.66666667,  0.33333333, -0.66666667],
         [ 0.66666667, -0.66666667,  0.33333333]])

  Rank Deficient Example
  ----------------------

  >>> X = np.array([[0.05,  0.05, 0.25, -0.25],
  ...                [0.25,  0.25, 0.05, -0.05],
  ...                [0.35,  0.35, 1.75, -1.75],
  ...                [1.75,  1.75, 0.35, -0.35],
  ...                [0.30, -0.30, 0.30,  0.30],
  ...                [0.40, -0.40, 0.40,  0.40]])
  >>> y = np.array([1, 2, 3, 4, 5, 6]).T

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,cond=5e-4)

  Results:

  >>> beta.T
  array([ 4.96666667, -2.83333333,  4.56666667,  3.23333333])
  >>> ss
  0.82666666666666699
  >>> cov
  array([[ 0.34027778, -0.15972222,  0.21527778,  0.28472222],
         [-0.15972222,  0.34027778, -0.28472222, -0.21527778],
         [ 0.21527778, -0.28472222,  0.34027778,  0.15972222],
         [ 0.28472222, -0.21527778,  0.15972222,  0.34027778]])
  >>> resids
  2.4800000000000009

  Rank, singular Values, and right-hand singular vectors:

  >>> rank
  3
  >>> np.allclose(s, np.array([  3.00000000e+00,   2.00000000e+00,   1.00000000e+00, 6.93889390e-18]))
  True
  >>> np.allclose(vt.sum(axis=0), [0,2,0,0])
  True
  >>> np.allclose(vt.sum(axis=1), [1,1,-1,1])
  True

  Weighted Example (diagonal)
  ---------------------------

  >>> X = np.array([[-0.57, -1.28,  -0.39],
  ...                [-1.93,  1.08,  -0.31],
  ...                [ 2.30,  0.24,  -0.40],
  ...                [-0.02,  1.03,  -1.43]])
  >>> y = np.array([1.31,-4.01,5.56,3.22]).T
  >>> w = np.array([0.5,1,2,5])

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,weights=w)

  Results:

  >>> beta.T
  array([ 2., -1., -3.])
  >>> np.allclose(ss, 2.36769765e-31)
  True
  >>> cov
  array([[ 0.07494261,  0.05450618,  0.04577265],
         [ 0.05450618,  0.55082517,  0.39779859],
         [ 0.04577265,  0.39779859,  0.38118818]])
  >>> np.allclose(resids, 2.36769765e-31)
  True

  Rank, singular Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 4.11371205,  3.81196479,  1.0665822 ])
  >>> vt
  array([[-0.13309546, -0.61438512,  0.77769951],
         [-0.98717469,  0.15197524, -0.04888411],
         [ 0.0881574 ,  0.77423153,  0.62673265]])

  Weighted Example (general)
  ---------------------------

  >>> X = np.array([[-0.57, -1.28,  -0.39],
  ...                [-1.93,  1.08,  -0.31],
  ...                [ 2.30,  0.24,  -0.40],
  ...                [-0.02,  1.03,  -1.43]])
  >>> y = np.array([1.31,-4.01,5.56,3.22]).T
  >>> w = np.diag([0.5,1,2,5])

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,weights=w)

  Results:

  >>> beta.T
  array([ 2., -1., -3.])
  >>> np.allclose(ss, 2.36769765e-31)
  True
  >>> cov
  array([[ 0.07494261,  0.05450618,  0.04577265],
         [ 0.05450618,  0.55082517,  0.39779859],
         [ 0.04577265,  0.39779859,  0.38118818]])
  >>> np.allclose(resids, 2.36769765e-31)
  True

  Rank, singular Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 4.11371205,  3.81196479,  1.0665822 ])
  >>> vt
  array([[-0.13309546, -0.61438512,  0.77769951],
         [-0.98717469,  0.15197524, -0.04888411],
         [ 0.0881574 ,  0.77423153,  0.62673265]])
  '''
  from scipy.linalg import lapack,calc_lwork

  a1,b1 = map(np.asarray,(a,b))

  if a1.ndim != 2:
    raise ValueError('incompatible design matrix dimensions')

  m,n = a1.shape

  if b1.ndim == 2:
    nrhs = b1.shape[1]
  else:
    nrhs = 1

  if m != b1.shape[0]:
    raise ValueError('incompatible dimensions')

  # Apply weights via 'whitening' or 'preconditioning'.  This is not always
  # stable for wildly variable weights, but works satisfactorily when
  # weights are well-conditioned.
  if weights is not None:
    # Weight vector
    k = weights.ndim
    if k==1 or (k==2 and 1 in weights.shape):
      weights = np.asarray(weights,dtype=float).reshape(-1)

      if weights.size != m:
        raise ValueError('incompatible weight vector dimensions')

      if not sqrtweights:
        weights = np.sqrt(weights)

      a1 = (weights*a1.T).T
      b1 = (weights*b1.T).T

    # Generalized covariance
    elif weights.ndim == 2:
      if weights.shape[0] != weights.shape[1] != m:
        raise ValueError('incompatible weight matrix dimensions')

      if not sqrtweights:
        weights = cholesky(weights,lower=False)

      a1 = np.dot(weights,a1)
      b1 = np.dot(weights,b1)

    else:
      raise ValueError('incompatible weight dimensions')

  gelss, = lapack.get_lapack_funcs(('gelss',),(a1,b1))

  if n>m:
    # extend b matrix as it will be filled with a larger solution
    b2 = np.zeros((n,nrhs), dtype=gelss.dtype)
    if b1.ndim==2:
      b2[:m,:] = b1
    else:
      b2[:m,0] = b1
    b1 = b2

  if gelss.module_name.startswith('flapack'):
    try:
      work = gelss(a1, b1, lwork=-1)[4]
      lwork = work[0].real.astype(np.int)
      v,x,s,rank,work,info = gelss(a1, b1, cond=cond, lwork=lwork,
                                           overwrite_a=0, overwrite_b=0)

    # There has to be a better way to do this...
    except Exception, e:
      lwork = calc_lwork.gelss(gelss.prefix,m,n,nrhs)[1]
      v,x,s,rank,info = gelss(a1, b1, cond=cond, lwork=lwork,
                              overwrite_a=0, overwrite_b=0)
  else:
    raise NotImplementedError('gelss not available from %s' % gelss.module_name)

  if info>0:
    raise LinAlgError('SVD did not converge')
  elif info<0:
    raise ValueError('illegal value in %d-th argument of %s.gelss' % (-info.gelss.module_name))

  resids = np.array([],   dtype=x.dtype)
  cov    = np.array([[]], dtype=x.dtype)
  ss     = np.array([],   dtype=x.dtype)
  vt     = v[:min(m,n)]
  beta   = x[:n]

  # FIXME: These should be un-squared residuals or else why bother?
  if n<m and rank==n:
    resids = (x[n:]**2).sum(axis=0)
  else:
    resids = ((b1-np.dot(a1,beta))**2).sum(axis=0)

  # Estimated residual variance
  if m>rank:
    ss = resids/(m-rank)

  # Truncated inverse squared singular values (vector)
  s2 = s**-2
  s2[rank:] = 0

  # Form covariance matrix based on SVD using truncated singular values:
  #   X     = U*S*V'                   U is unitary, so U'*U=U*U'=I;
  #  X'X    = V*S**2*V'                V is unitary, so V'*V=V*V'=I;
  #   I     = V*S**2*V'*(V*S**-2*V')   S is diagonal positive semi-definite
  # (X'X).I = V*S**-2*V'
  cov = np.dot(vt.T,(s2*vt.T).T)

  return beta,ss,cov,resids,rank,s,vt


def sqrtm(x):
  '''
  Compute the square root y of x such that y*y = x, where x is a general
  positive semi-definite matrix and * is matrix multiplication.

  Algorithm by Nicholas J. Higham, 'A New sqrtm for Matlab' (1999),
  Manchester Centre for Computational Mathematics Numerical Analysis
  Reports.  http://www.maths.man.ac.uk/~nareports/narep336.ps.gz

  Near-singular example:

  >>> x = ([[ 1.,      0.5,     0.3333,  0.25  ],
  ...       [ 0.5,     0.3333,  0.25,    0.2   ],
  ...       [ 0.3333,  0.25,    0.2,     0.1667],
  ...       [ 0.25,    0.2,     0.1667,  0.1429]])
  >>> y = sqrtm(x)
  >>> y
  array([[ 0.91145927,  0.33907405,  0.19263836,  0.13100093],
         [ 0.33907405,  0.35263109,  0.24780888,  0.18047398],
         [ 0.19263836,  0.24780888,  0.24041034,  0.20900738],
         [ 0.13100093,  0.18047398,  0.20900738,  0.22244957]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = np.array([[4.,1.,1.],
  ...              [2.,4.,1.],
  ...              [0.,1.,4.]])
  >>> y = sqrtm(x)
  >>> y
  array([[ 1.97119712,  0.23914631,  0.23914631],
         [ 0.51131184,  1.95468751,  0.2226367 ],
         [-0.03301922,  0.25565592,  1.98770673]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True
  '''
  A = np.asanyarray(x)
  if A.ndim!=2:
    raise ValueError('Non-matrix input to matrix function')

  # Refactor into complex Schur form
  T,Z = schur(A)
  T,Z = rsf2csf(T,Z)
  n,n = T.shape

  # Compute upper triangular square root R of T a column at a tim
  R = np.diag(np.sqrt(np.diag(T)))
  for j in xrange(n):
    for i in xrange(j-1,-1,-1):
      k = slice(i+1,j)
      R[i,j] = (T[i,j] - np.dot(R[i,k],R[k,j]))/(R[i,i] + R[j,j])

  X = np.dot(np.dot(Z,R),np.conjugate(Z.T)).astype(A.dtype.char)

  # XXX: Should be a warning
  if 0:
    nzeig = (np.diag(T)==0).any()
    if nzeig:
      print 'Matrix is singular and may not have a square root.'

  return X


def sqrtm_svd(x):
  '''
  Compute the square root y of x such that y*y = x, where x is a general
  positive semi-definite matrix and * is standard matrix multiplication.
  Computation is performed by singular value decomposition.

  Will not work for defective inputs (and no error is generated)!

  Near-singular example:

  >>> x = ([[ 1.,      0.5,     0.3333,  0.25  ],
  ...       [ 0.5,     0.3333,  0.25,    0.2   ],
  ...       [ 0.3333,  0.25,    0.2,     0.1667],
  ...       [ 0.25,    0.2,     0.1667,  0.1429]])
  >>> y = sqrtm_svd(x)
  >>> y
  array([[ 0.91145927,  0.33907405,  0.19263836,  0.13100093],
         [ 0.33907405,  0.35263109,  0.24780888,  0.18047398],
         [ 0.19263836,  0.24780888,  0.24041034,  0.20900738],
         [ 0.13100093,  0.18047398,  0.20900738,  0.22244957]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = np.array([[4.,1.,1.],
  ...               [2.,4.,1.],
  ...               [0.,1.,4.]])
  >>> y = sqrtm_svd(x)
  >>> y
  array([[ 1.94181569,  0.11454081,  0.37488366],
         [ 0.63982574,  1.94044655,  0.2169288 ],
         [-0.16206229,  0.28998723,  1.97233742]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  False
  '''
  u,s,vt = svd(x)
  return np.dot(u,((s**0.5)*vt.T).T)


def sqrtm_eig(x):
  '''
  Compute the square root y of x such that y*y = x, where x is a general
  positive semi-definite matrix and * is standard matrix multiplication.
  Computation is performed by eigenvalue decomposition.

  Will not work for defective inputs (and no error is generated)!

  Near-singular example:

  >>> x = ([[ 1.,      0.5,     0.3333,  0.25  ],
  ...       [ 0.5,     0.3333,  0.25,    0.2   ],
  ...       [ 0.3333,  0.25,    0.2,     0.1667],
  ...       [ 0.25,    0.2,     0.1667,  0.1429]])
  >>> y = sqrtm_eig(x)
  >>> y
  array([[ 0.91145927,  0.33907405,  0.19263836,  0.13100093],
         [ 0.33907405,  0.35263109,  0.24780888,  0.18047398],
         [ 0.19263836,  0.24780888,  0.24041034,  0.20900738],
         [ 0.13100093,  0.18047398,  0.20900738,  0.22244957]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = np.array([[4.,1.,1.],
  ...               [2.,4.,1.],
  ...               [0.,1.,4.]])
  >>> y = sqrtm_eig(x)
  >>> y
  array([[ 0.76018647,  1.01358196,  0.50679098],
         [ 1.01358196,  3.08349342, -1.0563295 ],
         [ 0.50679098, -1.0563295 ,  2.06991146]])
  >>> error = norm(x-np.dot(y,y))
  >>> error > 1e-14
  True
  '''
  from scipy.linalg import eig

  d,e = eig(x)
  d = (d**0.5).astype(float)
  return np.dot(e,(d*e).T).astype(float)


def sqrtm_symmetric(x,cond=1e-7):
  '''
  Compute the square root y of x such that y*y = x, where x is a symmetric
  positive semi-definite matrix and * is standard matrix multiplication.
  Computation is performed by symmetric eigenvalue decomposition.

  Near-singular example:

  >>> x = ([[ 1.,      0.5,     0.3333,  0.25  ],
  ...       [ 0.5,     0.3333,  0.25,    0.2   ],
  ...       [ 0.3333,  0.25,    0.2,     0.1667],
  ...       [ 0.25,    0.2,     0.1667,  0.1429]])
  >>> y = sqrtm_symmetric(x)
  >>> y
  array([[ 0.91145927,  0.33907405,  0.19263836,  0.13100093],
         [ 0.33907405,  0.35263109,  0.24780888,  0.18047398],
         [ 0.19263836,  0.24780888,  0.24041034,  0.20900738],
         [ 0.13100093,  0.18047398,  0.20900738,  0.22244957]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True
  '''
  from scipy.linalg import eigh

  d,e = eigh(x)
  d[d<cond] = 0
  return np.dot(e,((d**0.5)*e).T).astype(float)


def sqrtm_symmetric2(x):
  '''
  Compute the square root y of x such that y*y = x, where x is a symmetric
  positive semi-definite matrix and * is standard matrix multiplication.
  Computation is performed by singular value decomposition of the lower
  Cholesky factorization of x (algorithm by Golub and Van Loan).

  Near-singular example:

  >>> x = ([[ 1.,      0.5,     0.3333,  0.25  ],
  ...       [ 0.5,     0.3333,  0.25,    0.2   ],
  ...       [ 0.3333,  0.25,    0.2,     0.1667],
  ...       [ 0.25,    0.2,     0.1667,  0.1429]])
  >>> y = sqrtm_symmetric2(x)
  >>> y
  array([[ 0.91145927,  0.33907405,  0.19263836,  0.13100093],
         [ 0.33907405,  0.35263109,  0.24780888,  0.18047398],
         [ 0.19263836,  0.24780888,  0.24041034,  0.20900738],
         [ 0.13100093,  0.18047398,  0.20900738,  0.22244957]])
  >>> error = norm(x-np.dot(y,y))
  >>> error < 1e-14
  True
  '''
  l=cholesky(x,lower=1)
  u,s,vt = svd(l)
  return np.dot(u,(s*u).T)


def bdiag_to_mat(x):
  '''
  Expand a block matrix of block diagonal matrices stored as vectors into a
  dense matrix.  Each block must be of equal size.

  >>> w11 = np.array([16.,36.])
  >>> w22 = np.array([64.,25.])
  >>> w12 = np.array([-4.,-16.])
  >>> bdiag_to_mat([[w11,w12],[w12,w22]])
  array([[ 16.,   0.,  -4.,   0.],
         [  0.,  36.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])
  '''
  n = len(x)
  m = [ [None]*n for i in xrange(n) ]
  for i in xrange(n):
    for j in xrange(n):
      m[i][j] = np.diag(np.asarray(x[i][j],dtype=float).reshape(-1))
  return bmat2d(m)


def bdiag_copy(x):
  '''
  Return a copy of the specified block matrix of block diagonal matrices
  stored as vectors.

  >>> w11 = np.array([16.,36.])
  >>> w22 = np.array([64.,25.])
  >>> w12 = np.array([-4.,-16.])
  >>> w=[[w11,w12],[w12,w22]]
  >>> x=bdiag_copy(w)
  >>> w[0][0] *= 2
  >>> bdiag_to_mat(x)  # Not changed
  array([[ 16.,   0.,  -4.,   0.],
         [  0.,  36.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])
  >>> bdiag_to_mat(w)  # Changed
  array([[ 32.,   0.,  -4.,   0.],
         [  0.,  72.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])
  '''
  k = len(x)
  return [ [ x[i][j].copy() for j in xrange(k) ] for i in xrange(k) ]


def block_cholesky(w,lower=True):
  '''
  Efficiently compute the lower Cholesky decomposition of a block matrix
  composed entirely of diagonal matrices.  Input is given as a 2-dimensional
  list of diagonal vectors.  As the results are also in the form of a block
  matrix composed of diagonal vectors, output is also a 2-dimensional list
  of diagonal vectors.

  This function employed a method that is essentially the classical scalar
  Cholesky-Crout algorithm adapted to exploit the block structure of the
  input and the diagonal nature of each block element.

  The traditional algorithm costs n^3*m^3/3 FLOPS, where n is the number
  of row/column blocks and m is the diagonal length of each block.  This
  algorithm requires only n^3*m/3 FLOPS, a reduction over the usual
  algorithm of a factor of m^2 for suitably structured inputs.

  Example of a 2x2 block matrix of 2x2 diagonal matrices:

  >>> w11 = np.array([16.,36.])
  >>> w22 = np.array([64.,25.])
  >>> w12 = np.array([-4.,-16.])
  >>> w   = [[w11,w12],[w12,w22]]
  >>> x   = bdiag_to_mat(w)
  >>> x
  array([[ 16.,   0.,  -4.,   0.],
         [  0.,  36.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])

  Perform a lower Cholesky decomposition on the expanded form:

  >>> l = np.array(cholesky(x.astype(float),lower=1))
  >>> l
  array([[ 4.        ,  0.        ,  0.        ,  0.        ],
         [ 0.        ,  6.        ,  0.        ,  0.        ],
         [-1.        ,  0.        ,  7.93725393,  0.        ],
         [ 0.        , -2.66666667,  0.        ,  4.22952585]])

  Verify x=l*l':

  >>> norm(np.dot(l,l.T)-x) < 1e-14
  True

  Now do the same thing using our algorithm, still using the packed format:

  >>> l2 = bdiag_to_mat(block_cholesky(w))
  >>> l2
  array([[ 4.        ,  0.        ,  0.        ,  0.        ],
         [ 0.        ,  6.        ,  0.        ,  0.        ],
         [-1.        ,  0.        ,  7.93725393,  0.        ],
         [ 0.        , -2.66666667,  0.        ,  4.22952585]])

  Verify the factorization matches the expanded version above and solves x=l*l':

  >>> norm(l2-l) < 1e-14
  True
  >>> norm(np.dot(l2,l2.T)-x) < 1e-14
  True

  Derivation:

  Let A be a symmetric block matrix with Aij as the i,j'th block of size
  nxn.  Each Aij element is a vector of length m, representing the diagonal
  of a mxm matrix.  For n=3,

      |               |          |               |   |               |
      | A11  A12  A13 |          | L11   0    0  |   | L11  L21  L31 |
      |               |          |               |   |               |
  A = | A21  A22  A23 | = L*L' = | L21  L22   0  | * |  0   L22  L32 |
      |               |          |               |   |               |
      | A31  A32  A33 |          | L31  L32  L33 |   |  0    0   L33 |
      |               |          |               |   |               |

      |                                                   |
      | L11*L11      L11*L12              L11*L13         |
      |                                                   |
    = | L21*L11  L21*L21+L22*L22      L21*L31+L22*L32     |
      |                                                   |
      | L31*L11  L31*L12+L32*L22  L31+L31+L32*L32+L33*L33 |
      |                                                   |

  where + and * are defined as element-wise addition and multiplication over
  vectors.

  Elements of L can be determined in a manner akin to the Cholesky-Crout
  algorithm:

  L11 = sqrt(A11)
  L21 = A21/L11
  L31 = A31/L11
  L22 = sqrt(A22-L21*L21)
  L32 = (A23-L31*L12)/L22
  L33 = sqrt(A33-L31+L31+L32*L32)

  This 3x3 case can be generalized to the nxn, where computation starts from
  the upper left corner of the matrix A and proceeds to calculate the matrix
  column by column by:

             min(i,j)-1
  Lij = (Aij - sum      Lik*Ljk)/Ljj     for i>j
               k=1

                   i-1
  Lii = sqrt(Aii - sum Lik^2)
                   k=1
  '''
  n = len(w)
  m = len(w[0][0])
  c = [ [np.zeros(m)]*n for i in xrange(n) ]

  for j in xrange(n):
    for i in xrange(j,n):
      ww = np.asarray(w[i][j],dtype=float).reshape(-1)
      s = ww - np.sum(c[i][k]*c[j][k] for k in xrange(min(i,j)))
      if i==j:
        c[i][j] = s**0.5
      else:
        c[i][j] = s/c[j][j]

  if not lower:
    for i in xrange(n):
      for j in xrange(i):
        c[i][j],c[j][i] = c[j][i],c[i][j]

  return c


def block_inverse_from_triangular(u,lower=True):
  '''
  Efficiently compute the inverse of a block matrix composed entirely of
  diagonal matrices.  Input is the upper or lower Cholesky factor of the
  desired matrix, given as a 2-dimensional list of diagonal vectors.  As the
  results are also in the form of a block matrix composed of diagonal
  vectors, output is also a 2-dimensional list of diagonal vectors.

  NOTE: u is destroyed by this function and replaced by its inverse

  The algorithm employed is essentially the classical inverse of a
  triangular system adapted to exploit the block structure of the input and
  the diagonal nature of each block element.

  >>> w11 = np.array([16.,36.])
  >>> w22 = np.array([64.,25.])
  >>> w12 = np.array([-4.,-16.])
  >>> w   = [[w11,w12],[w12,w22]]
  >>> x   = bdiag_to_mat(w)
  >>> x
  array([[ 16.,   0.,  -4.,   0.],
         [  0.,  36.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])

  Find the inverse using the expanded form:

  >>> y = inv(x)
  >>> y
  array([[ 0.06349206,  0.        ,  0.00396825,  0.        ],
         [ 0.        ,  0.03881988,  0.        ,  0.02484472],
         [ 0.00396825,  0.        ,  0.01587302,  0.        ],
         [ 0.        ,  0.02484472,  0.        ,  0.05590062]])
  >>> norm(np.dot(x,y)-np.eye(4)) < 1e-14
  True

  Now with our algorithm, starting with a packed lower block cholesky
  decomposition:

  >>> l = block_cholesky(w)
  >>> z = bdiag_to_mat(block_inverse_from_triangular(l))
  >>> z
  array([[ 0.06349206,  0.        ,  0.00396825,  0.        ],
         [ 0.        ,  0.03881988,  0.        ,  0.02484472],
         [ 0.00396825,  0.        ,  0.01587302,  0.        ],
         [ 0.        ,  0.02484472,  0.        ,  0.05590062]])

  Verify the solution matches above and solves x*y=I:

  >>> norm(z-y) < 1e-14
  True
  >>> norm(np.dot(x,z)-np.eye(4)) < 1e-14
  True
  '''
  m = len(u[0][0])
  n = len(u)

  if lower:
    for i in xrange(n):
      for j in xrange(i):
        u[i][j],u[j][i] = u[j][i],u[i][j]

  # Invert upper triangular system, updating u in place
  for j in xrange(n):
    ujj     = 1/u[j][j]
    u[j][j] = ujj
    t       = [ u[k][j] for k in xrange(j) ]
    for i in xrange(j):
      u[i][j] = -ujj*u[i][i]*t[i]+np.sum(u[i][k]*t[k] for k in xrange(i+1,j-1))

  # Form a*a' = w^-1
  c = [ [np.zeros(m)]*n for i in xrange(n) ]
  for i in xrange(n):
    for j in xrange(i,n):
      c[i][j] = c[j][i] = np.sum(u[i][k]*u[j][k] for k in xrange(j,n))

  return c


def block_inverse(w):
  '''
  Efficiently compute the inverse of a block matrix composed entirely of
  diagonal matrices.  Input is given as a 2-dimensional list of diagonal
  vectors.  As the results are also in the form of a block matrix composed
  of diagonal vectors, output is also a 2-dimensional list of diagonal
  vectors.

  This is essentially the classical inverse of a triangular system adapted
  to exploit the block structure of the input and the diagonal nature of
  each block element.

  >>> w11 = np.array([16.,36.])
  >>> w22 = np.array([64.,25.])
  >>> w12 = np.array([-4.,-16.])
  >>> w   = [[w11,w12],[w12,w22]]
  >>> x   = bdiag_to_mat(w)
  >>> x
  array([[ 16.,   0.,  -4.,   0.],
         [  0.,  36.,   0., -16.],
         [ -4.,   0.,  64.,   0.],
         [  0., -16.,   0.,  25.]])

  Find the inverse using the expanded form:

  >>> y = inv(x)
  >>> y
  array([[ 0.06349206,  0.        ,  0.00396825,  0.        ],
         [ 0.        ,  0.03881988,  0.        ,  0.02484472],
         [ 0.00396825,  0.        ,  0.01587302,  0.        ],
         [ 0.        ,  0.02484472,  0.        ,  0.05590062]])
  >>> norm(np.dot(x,y)-np.eye(4)) < 1e-14
  True

  Now with our algorithm, starting with the packed input:

  >>> z = bdiag_to_mat(block_inverse(w))
  >>> z
  array([[ 0.06349206,  0.        ,  0.00396825,  0.        ],
         [ 0.        ,  0.03881988,  0.        ,  0.02484472],
         [ 0.00396825,  0.        ,  0.01587302,  0.        ],
         [ 0.        ,  0.02484472,  0.        ,  0.05590062]])

  Verify the solution matches above and solves x*y=I:

  >>> norm(z-y) < 1e-14
  True
  >>> norm(np.dot(x,z)-np.eye(4)) < 1e-14
  True
  '''
  a = block_cholesky(w,lower=False)
  return block_inverse_from_triangular(a,lower=False)


def grand_mean(X):
  return np.bmat([ np.ones((len(X),1)), X ]).A


def linreg(y, X, add_mean=False):
  '''
  Full-rank Example
  -----------------

  >>> X = np.array([[1.00,  0.60, 1.20, 3.90],
  ...               [1.00,  5.00, 4.00, 2.50],
  ...               [1.00,  1.00,-4.00,-5.50],
  ...               [1.00, -1.00,-2.00,-6.50],
  ...               [1.00, -4.20,-8.40,-4.80]])
  >>> y = np.array([3, 4, -1, -5, -1]).T

  Least-squares Fit:

  >>> beta,cov,ss = linreg(y, X)

  Results:

  >>> beta.T
  array([ 0.1517341 ,  0.9022158 , -0.8039499 ,  0.90558767])
  >>> np.allclose(ss, 0.28901734104046217)
  True
  >>> cov
  array([[ 0.4515896 , -0.15213552,  0.11721259, -0.0032113 ],
         [-0.15213552,  0.11347499, -0.08170984,  0.01441519],
         [ 0.11721259, -0.08170984,  0.08486762, -0.0297224 ],
         [-0.0032113 ,  0.01441519, -0.0297224 ,  0.0266895 ]])

  Rank Deficient Example
  ----------------------

  >>> X = np.array([[0.05,  0.05, 0.25, -0.25],
  ...               [0.25,  0.25, 0.05, -0.05],
  ...               [0.35,  0.35, 1.75, -1.75],
  ...               [1.75,  1.75, 0.35, -0.35],
  ...               [0.30, -0.30, 0.30,  0.30],
  ...               [0.40, -0.40, 0.40,  0.40]])
  >>> y = np.array([1, 2, 3, 4, 5, 6]).T

  Least-squares Fit (or lack thereof):

  >>> beta,cov,ss = linreg(y,X)

  Results (not ideal due to extremely ill-conditioned design)

  Parameters estimates are well off
  >>> np.allclose(beta.T, np.array([ 9. , -8. , -1.5, -4. ]))
  True

  Residual variance highly inflated -- should be ~0.83

  >>> np.allclose(ss, 5.9875)
  True

  Huge oscillating covariances, crazy negative variances
  >>> np.allclose(cov, np.array([[ -1.12589991e+15,  1.12589991e+15,  1.12589991e+15,  1.12589991e+15],
  ...                            [  1.12589991e+15, -1.12589991e+15, -1.12589991e+15, -1.12589991e+15],
  ...                            [  1.12589991e+15, -1.12589991e+15, -1.12589991e+15, -1.12589991e+15],
  ...                            [  1.12589991e+15, -1.12589991e+15, -1.12589991e+15, -1.12589991e+15]]))
  True
  '''
  y = np.asarray(y,dtype=float)
  X = np.asarray(X,dtype=float)

  #if y.ndim == 1:
  #  y = y[:,np.newaxis]

  # Add type specific means, if requested
  if add_mean:
    X = grand_mean(X)

  # n=observations and m=covariates
  n,m = X.shape

  # Fit using normal equations
  # FIXME: not a good idea for numerical stability when posed with
  # rank-deficient or ill-conditioned designs
  W  = inv(dot(X.T,X))
  b  = dot(W,dot(X.T,y))
  ss = ((y - dot(X,b))**2).sum()/(n-m)
  return b,W,ss


def logit(y, X, initial_beta=None, add_mean=False, max_iterations=50):
  '''
  Logistic/binomial regression

  Example outcome data taken from:

     http://luna.cas.usf.edu/~mbrannic/files/regression/Logistic.html

  Load and parse test data

  >>> from pkg_resources import resource_stream
  >>> D = []
  >>> indices = [0,1,2]
  >>> for line in resource_stream('glu','test/datasets/heart.dat'):
  ...   fields = line.split()
  ...   row = [ fields[i] for i in indices ]
  ...   row = map(float,row)
  ...   D.append(row)

  >>> D=np.asarray(D,dtype=float)

  Extract dependent variable (first column) and design array (following columns)

  >>> y,X = D[:,0].astype(int),D[:,1:]

  Compute the logistic model

  >>> L,b,W = logit(y,X,add_mean=True)

  Compute standard errors

  >>> stde = W.diagonal()**.5

  Compute Z-statistics

  >>> z = b.T/stde

  Show results

  >>> print 'obs =',len(D)
  obs = 20
  >>> print 'logL =',L
  logL = -9.41018298385
  >>> print '-2logL =',-2*L
  -2logL = 18.8203659677
  >>> print 'beta =',b.T
  beta = [[-6.36346994 -1.02410605  0.11904492]]
  >>> print 'stde =',stde
  stde = [ 3.21389766  1.17107845  0.05497905]
  >>> print 'Z =',z
  Z = [[-1.97998524 -0.87449825  2.16527797]]
  >>> print 'X2 =',z**2
  X2 = [[ 3.92034156  0.76474719  4.68842868]]

  Test model against empty model

  >>> print 'score test=%.6f df=%d' % GLogit(y,X,add_mean=True).score_test().test()
  score test=7.584906 df=2
  '''
  y = np.asarray(y,dtype=int)
  X = np.asarray(X,dtype=float)

  if y.ndim==1:
    y = y[:,np.newaxis]

  if y.shape[0] != X.shape[0]:
    raise ValueError('incompatible dimensions')

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
    b = np.zeros((m,1), dtype=float)
    b[0,0] = np.log(m1/(1-m1))

  # Initialize likelihood to -infinity
  L = -np.inf

  # Iterate until convergence
  for i in xrange(max_iterations):
    eta1 = dot(X,b)

    # Compute expected values and linear predictors
    e    = np.exp(eta1)
    mu0  = 1/(1+e)
    mu1  = e*mu0
    eta0 = np.log(mu0)

    # Compute likelihood
    newL = eta0.sum() + y.choose( (0,eta1) ).sum()

    d,L = newL-L,newL

    # Stop when likelihood converges
    if d < CONV:
      break

    # Form weights and quadrants of the expected information matrix
    w = mu1 * (1-mu1)

    # Update regression estimates
    z = (y - mu1)/w + eta1
    b,ss,W,resids,rank,s,vt = linear_least_squares(X,z,weights=w)

  else:
    raise LinAlgError('logit estimator failed to converge')

  # Return log-likelihood, regression estimates, and covariance matrix
  return L,b,W


class GLogit(object):
  '''
  Example trichotomous outcome data taken from 'Logistic Regression: A
  Self-Learning Text', David G. Kleinbaum, Springer-Verlag Telos, ISBN
  978-0387941424 based on dataset from Hill et al., 1995.  The study
  involves 288 women who had been diagnosed with endometrial cancer.

  Load and parse test data

  >>> from pkg_resources import resource_stream
  >>> D = []
  >>> indices = [4,5,3,6]
  >>> for line in resource_stream('glu','test/datasets/cancer.dat'):
  ...   fields = line.split()
  ...   row = [ fields[i] for i in indices ]
  ...   if '.' in row:
  ...     continue
  ...   D.append([float(f) for f in row])

  >>> D=np.asarray(D,dtype=float)

  Extract dependent variable (first column) and design matrix (following columns)

  >>> y,X = D[:,0].astype(int),D[:,1:]

  Compute the polytomous logit model for 3 outcomes

  >>> g = GLogit(y,X,add_mean=True,ref=0)
  >>> L,b,W = g.fit()

  Compute standard errors

  >>> stde = W.diagonal()**.5

  Compute Z-statistics

  >>> z = b.T/stde

  Show results

  >>> print 'obs =',len(D)
  obs = 286
  >>> print 'logL =',L
  logL = -247.202541257
  >>> print 'beta =',b.T
  beta = [[-1.88217996  0.98705921 -0.64389914  0.88946435 -1.20321646  0.28228558
    -0.10708623 -1.79131151]]
  >>> print 'ss =',stde
  ss = [ 0.40248124  0.41178981  0.34356069  0.52534808  0.31897576  0.32796594
    0.30673956  1.04647716]
  >>> print 'Z =',z
  Z = [[-4.67644149  2.39699764 -1.87419331  1.69309528 -3.77212504  0.86071615
    -0.34911126 -1.71175403]]
  >>> print 'X2 =',z**2
  X2 = [[ 21.86910498   5.74559769   3.51260056   2.86657163  14.22892731
      0.7408323    0.12187867   2.93010187]]

  Test model against empty model

  >>> print 'score test=%.6f df=%d' % g.score_test().test()
  score test=15.944228 df=6
  >>> print 'wald test=%.6f df=%d' % g.wald_test().test()
  wald test=13.942303 df=6
  >>> print 'lr test=%.6f df=%d' % g.lr_test().test()
  lr test=18.218406 df=6
  '''
  def __init__(self, y, X, ref=None, vars=None, add_mean=False, max_iterations=100):
    y = np.asarray(y,dtype=int)
    X = np.asarray(X,dtype=float)

    if y.ndim == 1:
      y = y[:,np.newaxis]

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    # Add type specific means, if requested
    if add_mean:
      X = grand_mean(X)

    categories = np.unique(y.ravel())

    if ref is None:
      ref = categories[0]
    else:
      if ref not in categories:
        raise ValueError('reference class for dependent variable not observed in data')
      categories = [ref] + [ r for r in categories if r!=ref ]

    if len(categories) < 2:
      raise ValueError('less than two dependent variable categories observed')

    # Recode y with ordinal categories from 0..k-1, where 0 is reference
    # FIXME: This should be a scalar loop
    y_ord = np.zeros_like(y)

    for i,v in enumerate(categories):
      y_ord[y==v] = i

    self.y     = y
    self.y_ord = y_ord
    self.X     = X
    self.vars  = vars
    self.categories = categories

    self.L    = None
    self.beta = None
    self.W    = None


  def fit(self, initial_beta=None, max_iterations=50):
    y = self.y_ord
    X = self.X
    k = len(self.categories)-1

    # n=observations and m=covariates
    n,m = X.shape

    # Initial estimate for mu1_0 and mu2_0
    mu = np.array( [(y==i+1).sum() for i in xrange(k) ], dtype=float)/n

    # Initial estimates of beta (b1,b2)
    if initial_beta is None:
      b = []
      for i in xrange(k):
        bi      = np.zeros( (m,1), dtype=float)
        bi[0,0] = np.log(mu[i]/(1-mu.sum()))
        b.append(bi)
    else:
      b = np.vsplit(initial_beta,k)

    # Initialize log likelihood to -infinity
    L = -np.inf

    # Iterate until convergence
    for iteration in xrange(max_iterations):
      eta = [ np.dot(X,bi) for bi in b ]

      # Compute expected values and linear predictors
      mu0  = 1/(1+np.exp(eta).sum(axis=0))
      eta0 = np.log(mu0)
      mu   = [ np.exp(eta[i])*mu0 for i in xrange(k) ]

      # Compute likelihood
      newL = eta0.sum() + y.choose( [0]+eta ).sum()

      d,L = newL-L,newL

      # Stop when likelihood converges
      if d < CONV:
        break

      # Form weights in block form
      w = self.update_weights(mu)
      w_min = abs(min(e.min() for row in w for e in row))

      if w_min < EPS:
        for row in w:
          for e in row:
            e[e<EPS] = 0

      #if w_min < EPS/100000:
      #  raise LinAlgError('glogit estimator failed due to extreme ill-conditioning (w_min=%e)' % w_min)

      # Compute upper Cholesky factor and inverse weights
      ww = block_cholesky(w,lower=False)
      uw = bdiag_copy(ww)
      iw = block_inverse_from_triangular(ww,lower=False)

      # Form augmented dependent vector and design matrix
      zz = self.dependent_vector(mu,eta,iw,uw)
      XX = self.design_matrix(uw)

      # Update model fit
      beta,ss,W,resids,rank,s,vt = linear_least_squares(XX,zz,cond=COND)
      b    = np.vsplit(beta,k)

    else:
      raise LinAlgError('glogit estimator failed to converge')

    # Return log-likelihood, regression estimates, and covariance matrix
    self.L       = L
    self.beta    = beta
    self.W       = W
    self.weights = w

    return L,beta,W


  def update_weights(self,mu):
    '''
    Form block matrix of block diagonal weights
    '''
    k = len(mu)
    w = [ [None]*k for i in xrange(k) ]
    for i in xrange(0,k):
      for j in xrange(i,k):
        if i == j:
          w[i][i] = mu[i]*(1-mu[i])
        else:
          w[i][j] = w[j][i] = -mu[i]*mu[j]
    return w


  def information_matrix(self,w):
    '''
    Form information matrix based on SVD using truncated singular values:
      X     = U*S*V'                   U is unitary, so that U'*U=U*U'=I;
     X'X    = V*S**2*V'                V is unitary, so that V'*V=V*V'=I;
    '''
    uw = block_cholesky(w,lower=False)
    XX = self.design_matrix(uw)

    # Fit linear model with y=1 in order to obtain s and vt from the SVD of XX
    # This is cheaper than computing the full SVD of X, since U is not required
    _,_,_,_,_,s,vt = linear_least_squares(XX,np.ones( (XX.shape[0],1) ))

    s[s<COND] = 0
    s = s**2

    info = np.dot(vt.T,(s*vt.T).T)

    return info


  def covariance_matrix(self,w):
    '''
    Form covariance matrix based on SVD using truncated singular values:
      X     = U*S*V'                   U is unitary, so U'*U=U*U'=I;
     X'X    = V*S**2*V'                V is unitary, so V'*V=V*V'=I;
      I     = V*S**2*V'*(V*S**-2*V')   S is diagonal positive semi-definite
    (X'X).I = V*S**-2*V'
    '''
    uw = block_cholesky(w,lower=False)
    XX = self.design_matrix(uw)

    # Fit linear model with y=1 in order to obtain the covariance matrix
    # This is cheaper than computing the full SVD of X, since U is not required
    _,_,W,_,_,_,_ = linear_least_squares(XX,np.ones( (XX.shape[0],1) ))

    return W


  def dependent_vector(self,mu,eta,iw,uw):
    '''
    Form new dependent vector by pre-conditioning the sum of the decorrelated
    scores and linear predictors by the Cholesky factor of the
    regression weights.
    '''
    k  = len(iw)
    y  = self.y_ord
    yy = [ ((y==i+1) - mu[i]).reshape(-1) for i in xrange(k) ]
    zz = [ [eta[i]+np.sum((iw[i][j]*yy[j]) for j in xrange(k)).reshape(-1,1)] for i in xrange(k) ]
    return np.bmat([ [np.sum(uw[i][j]*zz[j][0].T for j in xrange(i,k)).T] for i in xrange(k) ]).A


  def design_matrix(self,uw):
    '''
    Form augmented design matrix by pre-conditioning the sum of the
    decorrelated scores and linear predictors by the Cholesky factor of the
    regression weights.
    '''
    X  = self.X
    k  = len(uw)
    z  = np.zeros_like(X)
    XX = [ [z]*k for i in xrange(k) ]
    for i in xrange(k):
      for j in xrange(i,k):
        XX[i][j] = (uw[i][j]*X.T).T
    return bmat2d(XX)


  def score_test(self,parameters=None,indices=None):
    return GLogitScoreTest(self,parameters=parameters,indices=indices)


  def wald_test(self,parameters=None,indices=None):
    return GLogitWaldTest(self,parameters=parameters,indices=indices)


  def lr_test(self,parameters=None,indices=None):
    return GLogitLRTest(self,parameters=parameters,indices=indices)


class GLogitScoreTest(object):
  def __init__(self,model,parameters=None,indices=None):
    y = model.y
    X = model.X
    k = len(model.categories)-1
    n = X.shape[1]

    if parameters is None and indices is None:
      indices  = [ j*n+i for i in range(1,n)
                         for j in range(k) ]
    elif parameters is not None:
      pindices = [ j*n+i for i in parameters
                         for j in range(k) ]
      indices  = sorted(set(indices or []) | set(pindices))

    # Form null design matrix by finding the rows of the design matrix
    # that are to be scored and excluding those from the null model
    design_indices = set(i%n for i in indices)
    null_indices   = [ i for i in xrange(n) if i not in design_indices ]
    X_null = X[:,null_indices]

    # Fit null model
    self.null = GLogit(y,X_null)
    self.null.fit()

    # Augment null beta with zeros for all parameters to be tested
    null_indices = (i for i in xrange(n*k) if i not in indices)
    beta = np.zeros((X.shape[1]*k,1), dtype=float)
    for b,i in izip(self.null.beta.ravel(),null_indices):
      beta[i,0] = b

    self.beta    = beta
    self.model   = model
    self.indices = indices

    self.build_scores()


  def build_scores(self):
    y = self.model.y_ord
    X = self.model.X
    k = len(self.model.categories)-1

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    # Form linear predictors and expected values for the current data
    # using the estimates under the null
    b = np.vsplit(self.beta,k)

    eta = [ dot(X,bi) for bi in b ]

    # Compute expected values and linear predictors
    mu0  = 1/(1+np.sum(np.exp(etai) for etai in eta))
    #eta0 = np.log(mu0)
    mu   = [ np.exp(eta[i])*mu0 for i in xrange(k) ]

    # Form weights and blocks of the expected information matrix
    self.w = self.model.update_weights(mu)

    # Build information matrix and compute the inverse
    self.W = self.model.covariance_matrix(self.w)

    # Compute scores
    self.scores = np.bmat([ (y==i+1) - mu[i] for i in xrange(k) ]).A


  def test(self):
    scores  = self.scores
    indices = self.indices
    W       = self.W
    df      = len(indices)
    X       = self.model.X
    n       = X.shape[1]

    # Compute score vectors with the full design matrix
    indices = np.array(indices, dtype=int)
    s       = (X[:,indices%n]*scores[:,indices//n]).sum(axis=0)

    # Extract the covariance matrix of the score parameters
    V = W.take(indices,axis=0).take(indices,axis=1)

    # Compute score test statistic
    x2 = float(dot(s,dot(V,s.T)))

    return x2,df


class GLogitWaldTest(object):
  '''
  Perform a simplified Wald test on a GLogit object

  Tests linear hypotheses of the form Cp = 0 in using the Wald test, where C
  is a mxn matrix of linear contrasts, m is the number of contrasts, n is
  the number of parameters in the GLogit model.

  The test statistic is W = (Cp)' [C cov(p) C']^-1 (Cp).

  W distributed is approximately chi2 with m degrees of freedom.

  Contrasts are currently limited to one or more single parameters,
  specified as a list of indices from 0..n-1.  This limitation will be
  lifted soon.
  '''
  def __init__(self,model,parameters=None,indices=None):
    k = len(model.categories)-1
    n = model.X.shape[1]

    if parameters is None and indices is None:
      indices  = [ j*n+i for i in range(1,n)
                         for j in range(k) ]
    elif parameters is not None:
      pindices = [ j*n+i for i in parameters
                         for j in range(k) ]
      indices  = sorted(set(indices or []) | set(pindices))

    if model.beta is None:
      model.fit()

    self.model   = model
    self.indices = indices

  def test(self):
    indices = self.indices
    beta    = self.model.beta
    df      = len(indices)

    # FIXME: Use robust inverse of covariance matrix
    b  = beta[indices,:]
    w  = inv(self.model.covariance_matrix(self.model.weights)[indices,:][:,indices])
    x2 = float(dot(b.T,dot(w,b)))

    return x2,df


class GLogitLRTest(object):
  '''
  Compute the likelihood ratio test statistic of a given logit model and a
  list of parameter indices to test.  The null model removing those indices
  is fit and -2(l_null-l_model) is distributed as a chi-squared distribution
  with len(indices) degrees of freedom.
  '''
  def __init__(self,model,initial_beta=None,parameters=None,indices=None):
    y = model.y
    X = model.X
    k = len(model.categories)-1
    m,n = X.shape

    if parameters is None and indices is None:
      indices  = [ j*n+i for i in range(1,n)
                         for j in range(k) ]
    elif parameters is not None:
      pindices = [ j*n+i for i in parameters
                         for j in range(k) ]
      indices  = sorted(set(indices or []) | set(pindices))

    # Form null design matrix by finding the rows of the design matrix
    # that are to be scored and excluding those from the null model
    design_indices = set(i%n for i in indices)
    null_indices   = [ i for i in xrange(n) if i not in design_indices ]
    X_null = X[:,null_indices]

    # Fit null model
    self.null = GLogit(y,X_null)
    self.null.fit()

    self.model   = model
    self.indices = indices

  def test(self):
    return -2*(self.null.L-self.model.L),len(self.indices)


class Linear(object):
  '''
  Full-rank Example
  -----------------

  >>> X = np.array([[ 0.60, 1.20, 3.90],
  ...                [ 5.00, 4.00, 2.50],
  ...                [ 1.00,-4.00,-5.50],
  ...                [-1.00,-2.00,-6.50],
  ...                [-4.20,-8.40,-4.80]])
  >>> y = np.array([3, 4, -1, -5, -1]).T

  Least-squares Fit:

  >>> l=Linear(y,X,add_mean=True)
  >>> x=l.fit()

  Results:

  >>> np.allclose(l.beta.T, np.array([[ 0.1517341 ,  0.9022158 , -0.8039499 ,  0.90558767]]))
  True
  >>> np.allclose(l.ss,0.28901734104046128)
  True
  >>> np.allclose(l.W, np.array([[ 0.4515896 , -0.15213552,  0.11721259, -0.0032113 ],
  ...                          [-0.15213552,  0.11347499, -0.08170984,  0.01441519],
  ...                          [ 0.11721259, -0.08170984,  0.08486762, -0.0297224 ],
  ...                          [-0.0032113 ,  0.01441519, -0.0297224 ,  0.0266895 ]]))
  True
  >>> np.allclose(l.resids, 0.28901734104046128)
  True
  >>> np.allclose(l.L, 0.03207358773597957)
  True
  >>> np.allclose(l.score_test().test(), (3.977768,3))
  True
  >>> np.allclose(l.wald_test().test(), (178.920000,3))
  True
  >>> np.allclose(l.lr_test().test(), (25.962562,3))
  True

  Check results using R:

  >>> import rpy
  >>> rpy.set_default_mode(rpy.NO_CONVERSION)
  >>> data = rpy.r.data_frame(x1=X[:,0],x2=X[:,1],x3=X[:,2],y=y)
  >>> linear_model = rpy.r.lm(rpy.r("y ~ x1+x2+x3"), data=data)
  >>> null_model   = rpy.r.lm(rpy.r("y ~ 1"), data=data)
  >>> rpy.set_default_mode(rpy.BASIC_CONVERSION)
  >>> print ', '.join('%s=%.8f' % i for i in sorted(linear_model.as_py()['coefficients'].iteritems()))
  (Intercept)=0.15173410, x1=0.90221580, x2=-0.80394990, x3=0.90558767
  >>> summary=rpy.r.summary(linear_model)
  >>> print summary['cov.unscaled']
  [[ 0.4515896  -0.15213552  0.11721259 -0.0032113 ]
   [-0.15213552  0.11347499 -0.08170984  0.01441519]
   [ 0.11721259 -0.08170984  0.08486762 -0.0297224 ]
   [-0.0032113   0.01441519 -0.0297224   0.0266895 ]]
  >>> print summary['sigma']**2
  0.28901734104
  >>> print rpy.r.logLik(linear_model)
  0.032073587736
  >>> print -2*(rpy.r.logLik(null_model)-rpy.r.logLik(linear_model))
  25.9625615383

  Rank Deficient Example
  ----------------------

  >>> X = np.array([[0.05,  0.05, 0.25, -0.25],
  ...                [0.25,  0.25, 0.05, -0.05],
  ...                [0.35,  0.35, 1.75, -1.75],
  ...                [1.75,  1.75, 0.35, -0.35],
  ...                [0.30, -0.30, 0.30,  0.30],
  ...                [0.40, -0.40, 0.40,  0.40]])
  >>> y = np.array([1, 2, 3, 4, 5, 6]).T

  Least-squares Fit:

  >>> l=Linear(y,X,add_mean=True)
  >>> x=l.fit()

  Results:

  >>> np.allclose(l.beta.T, np.array([[ 1.18918919,  3.81711712, -2.31801802,  3.41711712,  2.71801802]]))
  True
  >>> np.allclose(l.ss, 0.19351351351351354)
  True
  >>> np.allclose(l.W, np.array([[ 0.67567568, -0.65315315,  0.29279279, -0.65315315, -0.29279279],
  ...                          [-0.65315315,  0.97165916, -0.44275526,  0.84665916,  0.56775526],
  ...                          [ 0.29279279, -0.44275526,  0.46715465, -0.56775526, -0.34215465],
  ...                          [-0.65315315,  0.84665916, -0.56775526,  0.97165916,  0.44275526],
  ...                          [-0.29279279,  0.56775526, -0.34215465,  0.44275526,  0.46715465]]))
  True
  >>> np.allclose(l.resids, 0.38702702702702707)
  True
  >>> np.allclose(l.L, -0.29057053820869672)
  True

  >>> np.allclose(l.score_test().test(), (4.889421,4))
  True
  >>> np.allclose(l.wald_test().test(), (-201.10990990991002, 4))
  True
  >>> np.allclose(l.lr_test().test(), (22.868770,4))
  True

  Check results using R (we do considerably better, since in R x4 is aliased
  and results are reported for a reduced parameter space without x4.  The
  quality of model fit is identical, though the interpretation is less
  obvious):

  >>> rpy.set_default_mode(rpy.NO_CONVERSION)
  >>> data = rpy.r.data_frame(x1=X[:,0],x2=X[:,1],x3=X[:,2],x4=X[:,3],y=y)
  >>> linear_model = rpy.r.lm(rpy.r("y ~ x1+x2+x3+x4"), data=data)
  >>> null_model   = rpy.r.lm(rpy.r("y ~ 1"), data=data)
  >>> rpy.set_default_mode(rpy.BASIC_CONVERSION)
  >>> print ', '.join('%s=%.8f' % i for i in sorted(linear_model.as_py()['coefficients'].iteritems()))
  (Intercept)=1.18918919, x1=6.53513514, x2=-5.03603604, x3=0.69909910, x4=nan
  >>> summary=rpy.r.summary(linear_model)
  >>> print summary['cov.unscaled']
  [[ 0.67567568 -0.94594595  0.58558559 -0.36036036]
   [-0.94594595  2.57432432 -1.81981982  0.2545045 ]
   [ 0.58558559 -1.81981982  1.61861862 -0.2012012 ]
   [-0.36036036  0.2545045  -0.2012012   0.5533033 ]]
  >>> print summary['sigma']**2
  0.193513513514
  >>> print rpy.r.logLik(linear_model)
  -0.290570538209
  >>> print -2*(rpy.r.logLik(null_model)-rpy.r.logLik(linear_model))
  22.8687697922
  '''
  def __init__(self, y, X, vars=None, add_mean=False):
    y = np.asarray(y,dtype=float)
    X = np.asarray(X,dtype=float)

    if y.ndim == 1:
      y = y[:,np.newaxis]

    assert y.shape[1] == 1
    assert y.shape[0] == X.shape[0]

    # Add grand mean, if requested
    if add_mean:
      X = grand_mean(X)

    self.y    = y
    self.X    = X
    self.vars = vars
    self.L    = None
    self.beta = None
    self.W    = None
    self.ss   = None
    self.rank = None

    self.singular_values        = None
    self.right_singlular_matrix = None

  def fit(self):
    y = self.y
    X = self.X

    # n=observations and m=covariates
    n,m = X.shape

    beta,ss,W,resids,rank,s,vt = linear_least_squares(X,y)

    L      = -(n/2.)*(1+np.log(2*pi)-np.log(n)+np.log(resids.sum()))

    self.L      = L
    self.beta   = beta
    self.W      = W
    self.ss     = ss.sum()
    self.resids = resids.sum()
    self.rank   = rank
    self.singular_values = s
    self.right_singlular_matrix = vt.T

    return L,beta,W,ss

  def p_values(self, phred=False):
    y     = self.y
    b     = self.beta.reshape(-1)
    stde  = (self.ss*self.W.diagonal())**0.5
    t     = b/stde
    n,m   = self.X.shape
    p     = 2*scipy.stats.distributions.t.cdf(-abs(t),n-m)

    if phred:
      p   = p.clip(1e-99,1)
      p   = (-10*np.log10(p)).astype(int).clip(0,99)

    return p

  def r2(self):
    ss_t  = np.var(self.y,ddof=1)
    r2    = 1 - self.ss/ss_t
    return r2

  def score_test(self,parameters=None,indices=None):
    return LinearScoreTest(self,parameters=parameters,indices=indices)

  def wald_test(self,parameters=None,indices=None):
    return LinearWaldTest(self,parameters=parameters,indices=indices)

  def lr_test(self,parameters=None,indices=None):
    return LinearLRTest(self,parameters=parameters,indices=indices)


class LinearScoreTest(object):
  def __init__(self,model,parameters=None,indices=None):
    y = model.y
    X = model.X
    n = X.shape[1]

    if parameters is None and indices is None:
      indices = range(1,n)
    else:
      indices = sorted(set(indices or []) | set(parameters or []))

    design_indices = set(indices)
    null_indices = [ i for i in xrange(n) if i not in design_indices ]
    X_null = X[:,null_indices]

    # Fit null model
    self.null = Linear(y,X_null)
    self.null.fit()

    # Augment null beta with zeros for all parameters to be tested
    beta = np.zeros((n,1), dtype=float)
    for b,i in izip(self.null.beta.ravel(),null_indices):
      beta[i,0] = b

    self.model   = model
    self.beta    = beta
    self.indices = indices
    self.scores  = model.y - dot(X,beta)


  def test(self):
    scores  = self.scores
    indices = self.indices
    W       = self.model.W
    df      = len(indices)
    X       = self.model.X
    ss      = self.null.ss

    # Compute score vectors with the full design matrix
    indices = np.asarray(indices, dtype=int)
    s = (X[:,indices]*scores).sum(axis=0)

    # Extract the covariance matrix of the score parameters
    V = W[indices,:][:,indices]/ss

    # Compute score test statistic
    x2 = float(dot(s,dot(V,s.T)))

    return x2,df


class LinearWaldTest(object):
  '''
  Perform a simplified Wald test on a Linear model

  Tests linear hypotheses of the form Cp = 0 in using the Wald test, where C
  is a mxn matrix of linear contrasts, m is the number of contrasts, n is
  the number of parameters in the Linear model.

  The test statistic is W = (Cp)' [C cov(p) C']^-1 (Cp).

  W distributed is approximately chi2 with m degrees of freedom.

  Contrasts are currently limited to one or more single parameters,
  specified as a list of indices from 0..n-1.  This limitation will be
  lifted soon.
  '''
  def __init__(self,model,parameters=None,indices=None):
    if parameters is None and indices is None:
      n = model.X.shape[1]
      indices = range(1,n)
    else:
      indices = sorted(set(indices or []) | set(parameters or []))

    self.model   = model
    self.indices = indices


  def test(self):
    indices = self.indices
    beta    = self.model.beta
    W       = self.model.W
    ss      = self.model.ss
    df      = len(indices)

    # FIXME: Use robust inverse of covariance matrix
    b  = beta[indices,:]
    w  = inv(W[indices,:][:,indices])/ss
    x2 = float(dot(b.T,dot(w,b)))

    return x2,df


class LinearLRTest(object):
  '''
  Compute the likelihood ratio test statistic of a given linear model and a
  list of parameter indices to test.  The null model removing those indices
  is fit and -2(l_null-l_model) is distributed as a chi-squared distribution
  with len(indices) degrees of freedom.
  '''
  def __init__(self,model,parameters=None,indices=None):
    y = model.y
    X = model.X
    n = X.shape[1]

    if parameters is None and indices is None:
      indices = range(1,n)
    else:
      indices = sorted(set(indices or []) | set(parameters or []))

    design_indices = set(indices)
    null_indices = [ i for i in xrange(n) if i not in design_indices ]
    X_null = X[:,null_indices]

    # Fit null model
    self.null = Linear(y,X_null)
    self.null.fit()

    self.model   = model
    self.indices = indices


  def test(self):
    return -2*(self.null.L-self.model.L),len(self.indices)


if __name__ == '__main__':
  import doctest
  doctest.testmod()
