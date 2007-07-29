# -*- coding: utf-8 -*-
'''
File:          logit.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from math         import pi
from itertools    import izip

# Import Python numerical libraries (from http://scipy.org/)
from numpy        import array,matrix,asmatrix,asarray,asarray_chkfinite,asanyarray,zeros,zeros_like,ones,\
                         where,exp,log,vsplit,unique,bmat,inf,sum,dot,transpose,diag,sqrt,conjugate
from scipy        import linalg
from scipy.linalg import lapack,calc_lwork,eig,eigh,svd,cholesky,norm,LinAlgError

from utils        import tally

CONV = 1e-8


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

  >>> X = matrix([[ 0.60, 1.20, 3.90],
  ...             [ 5.00, 4.00, 2.50],
  ...             [ 1.00,-4.00,-5.50],
  ...             [-1.00,-2.00,-6.50],
  ...             [-4.20,-8.40,-4.80]])
  >>> y = matrix([3, 4, -1, -5, -1]).T

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,cond=5e-4)

  Results:

  >>> beta.T
  array([[ 0.95333333, -0.84333333,  0.90666667]])
  >>> ss
  array([ 0.17])
  >>> cov
  array([[ 0.06222222, -0.04222222,  0.01333333],
         [-0.04222222,  0.05444444, -0.02888889],
         [ 0.01333333, -0.02888889,  0.02666667]])
  >>> resids
  array([ 0.34])

  Rank, singluar Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 15.,   6.,   3.])
  >>> vt
  array([[-0.33333333, -0.66666667, -0.66666667],
         [ 0.66666667,  0.33333333, -0.66666667],
         [ 0.66666667, -0.66666667,  0.33333333]])

  Rank Deficient Example
  ----------------------

  >>> X = matrix([[0.05,  0.05, 0.25, -0.25],
  ...             [0.25,  0.25, 0.05, -0.05],
  ...             [0.35,  0.35, 1.75, -1.75],
  ...             [1.75,  1.75, 0.35, -0.35],
  ...             [0.30, -0.30, 0.30,  0.30],
  ...             [0.40, -0.40, 0.40,  0.40]])
  >>> y = matrix([1, 2, 3, 4, 5, 6]).T

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,cond=5e-4)

  Results:

  >>> beta.T
  array([[ 4.96666667, -2.83333333,  4.56666667,  3.23333333]])
  >>> ss
  array([ 0.82666667])
  >>> cov
  array([[ 0.34027778, -0.15972222,  0.21527778,  0.28472222],
         [-0.15972222,  0.34027778, -0.28472222, -0.21527778],
         [ 0.21527778, -0.28472222,  0.34027778,  0.15972222],
         [ 0.28472222, -0.21527778,  0.15972222,  0.34027778]])
  >>> resids
  array([ 2.48])

  Rank, singluar Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([  3.00000000e+00,   2.00000000e+00,   1.00000000e+00,
           6.93889390e-18])

  >>> vt
  array([[-0.5, -0.5, -0.5,  0.5],
         [-0.5, -0.5,  0.5, -0.5],
         [ 0.5, -0.5,  0.5,  0.5],
         [ 0.5, -0.5, -0.5, -0.5]])

  Weighted Example (diagonal)
  ---------------------------

  >>> X = matrix([[-0.57, -1.28,  -0.39],
  ...             [-1.93,  1.08,  -0.31],
  ...             [ 2.30,  0.24,  -0.40],
  ...             [-0.02,  1.03,  -1.43]])
  >>> y = matrix([1.31,-4.01,5.56,3.22]).T
  >>> w = matrix([0.5,1,2,5])

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,weights=w)

  Results:

  >>> beta.T
  array([[ 2., -1., -3.]])
  >>> ss
  array([  7.88860905e-31])
  >>> cov
  array([[ 0.07494261,  0.05450618,  0.04577265],
         [ 0.05450618,  0.55082517,  0.39779859],
         [ 0.04577265,  0.39779859,  0.38118818]])
  >>> resids
  array([  7.88860905e-31])

  Rank, singluar Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 4.11371205,  3.81196479,  1.0665822 ])
  >>> vt
  array([[ 0.13309546,  0.61438512, -0.77769951],
         [-0.98717469,  0.15197524, -0.04888411],
         [ 0.0881574 ,  0.77423153,  0.62673265]])

  Weighted Example (general)
  ---------------------------

  >>> X = matrix([[-0.57, -1.28,  -0.39],
  ...             [-1.93,  1.08,  -0.31],
  ...             [ 2.30,  0.24,  -0.40],
  ...             [-0.02,  1.03,  -1.43]])
  >>> y = matrix([1.31,-4.01,5.56,3.22]).T
  >>> w = diag([0.5,1,2,5])

  Least-squares Fit:

  >>> beta,ss,cov,resids,rank,s,vt = linear_least_squares(X,y,weights=w)

  Results:

  >>> beta.T
  array([[ 2., -1., -3.]])
  >>> ss
  array([  7.88860905e-31])
  >>> cov
  array([[ 0.07494261,  0.05450618,  0.04577265],
         [ 0.05450618,  0.55082517,  0.39779859],
         [ 0.04577265,  0.39779859,  0.38118818]])
  >>> resids
  array([  7.88860905e-31])

  Rank, singluar Values, and right-hand singular vectors:

  >>> rank
  3
  >>> s
  array([ 4.11371205,  3.81196479,  1.0665822 ])
  >>> vt
  array([[ 0.13309546,  0.61438512, -0.77769951],
         [-0.98717469,  0.15197524, -0.04888411],
         [ 0.0881574 ,  0.77423153,  0.62673265]])
  '''
  a1,b1 = map(asarray_chkfinite,(a,b))

  if len(a1.shape) != 2:
    raise ValueError('incompatible design matrix dimensions')

  m,n = a1.shape

  if len(b1.shape)==2:
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
    k = len(weights.shape)
    if k==1 or (k==2 and 1 in weights.shape):
      weights = asarray(weights).reshape(-1)

      if weights.size != m:
        raise ValueError('incompatible weight vector dimensions')

      if not sqrtweights:
        weights = sqrt(weights)

      a1 = transpose(weights*transpose(a1))
      b1 = transpose(weights*transpose(b1))

    # Generalized covariance
    elif len(weights.shape) == 2:
      if weights.shape[0] != weights.shape[1] != m:
        raise ValueError('incompatible weight matrix dimensions')

      if not sqrtweights:
        weights = cholesky(weights,lower=False)

      a1 = dot(weights,a1)
      b1 = dot(weights,b1)

    else:
      raise ValueError('incompatible weight dimensions')

  gelss, = lapack.get_lapack_funcs(('gelss',),(a1,b1))

  if n>m:
    # extend b matrix as it will be filled with a larger solution
    b2 = zeros((n,nrhs), dtype=gelss.dtype)
    if len(b1.shape)==2:
      b2[:m,:] = b1
    else:
      b2[:m,0] = b1
    b1 = b2

  if gelss.module_name.startswith('flapack'):
    lwork = calc_lwork.gelss(gelss.prefix,m,n,nrhs)[1]
    v,x,s,rank,info = gelss(a1, b1, cond=cond, lwork=lwork,
                            overwrite_a=0, overwrite_b=0)
  else:
    raise NotImplementedError('gelss not available from %s' % gelss.module_name)

  if info>0:
    raise LinAlgError('SVD did not converge')
  elif info<0:
    raise ValueError('illegal value in %d-th argument of %s.gelss' % (-info.gelss.module_name))

  resids = array([],   dtype=x.dtype)
  cov    = array([[]], dtype=x.dtype)
  ss     = array([],   dtype=x.dtype)
  vt     = v[:min(m,n)]
  beta   = x[:n]

  if n<m and rank==n:
    resids = (x[n:]**2).sum(axis=0)
  else:
    resids = ((b-dot(a1,beta)).A**2).sum(axis=0)

  # Estimated residual variance
  if m>rank:
    ss = resids/(m-rank)

  # Truncated inverse squared singular values (vector)
  s2 = s**-2
  s2[rank:] = 0

  # Form covariance matrix based on SVD using truncated singular values:
  #   X     = U*S*V'                   U is unitary, so that U'*U=U*U'=I;
  #  X'X    = V*S**2*V'                V is unitary, so that V'*V=V*V'=I;
  #   I     = V*S**2*V'*(V*S**-2*V')   S is diagonal positive semi-definite
  # (X'X).I = V*S**-2*V'
  v   = transpose(vt)
  cov = dot(v,transpose(s2*v))

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
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = array([[4.,1.,1.],
  ...            [2.,4.,1.],
  ...            [0.,1.,4.]])
  >>> y = sqrtm(x)
  >>> y
  array([[ 1.97119712,  0.23914631,  0.23914631],
         [ 0.51131184,  1.95468751,  0.2226367 ],
         [-0.03301922,  0.25565592,  1.98770673]])
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True
  '''
  A = asanyarray(x)
  if len(A.shape)!=2:
    raise ValueError('Non-matrix input to matrix function')

  # Refactor into complex Schur form
  T,Z = linalg.schur(A)
  T,Z = linalg.rsf2csf(T,Z)
  n,n = T.shape

  # Compute upper triangular square root R of T a column at a tim
  R = diag(sqrt(diag(T)))
  for j in range(n):
    for i in range(j-1,-1,-1):
      k = slice(i+1,j)
      R[i,j] = (T[i,j] - dot(R[i,k],R[k,j]))/(R[i,i] + R[j,j])

  X = dot(dot(Z,R),conjugate(transpose(Z))).astype(A.dtype.char)

  # XXX: Should be a warning
  if 0:
    nzeig = (diag(T)==0).any()
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
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = array([[4.,1.,1.],
  ...            [2.,4.,1.],
  ...            [0.,1.,4.]])
  >>> y = sqrtm_svd(x)
  >>> y
  array([[ 1.94181569,  0.11454081,  0.37488366],
         [ 0.63982574,  1.94044655,  0.2169288 ],
         [-0.16206229,  0.28998723,  1.97233742]])
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  False
  '''
  u,s,vt = svd(x)
  return dot(u,transpose((s**0.5)*transpose(vt)))


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
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True

  Defective example:

  >>> x = array([[4.,1.,1.],
  ...            [2.,4.,1.],
  ...            [0.,1.,4.]])
  >>> y = sqrtm_eig(x)
  >>> y
  array([[ 0.76018647,  1.01358196,  0.50679098],
         [ 1.01358196,  3.08349342, -1.0563295 ],
         [ 0.50679098, -1.0563295 ,  2.06991146]])
  >>> error = norm(x-dot(y,y))
  >>> error > 1e-14
  True
  '''
  d,e = eig(x)
  d = (d**0.5).astype(float)
  return dot(e,transpose(d*e)).astype(float)


def sqrtm_symmetric(x,cond=1e-7):
  '''
  Compute the square root y of x such that y*y = x, where x is a symmetic
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
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True
  '''
  d,e = eigh(x)
  d[d<cond] = 0
  return dot(e,transpose((d**0.5)*e)).astype(float)


def sqrtm_symmetric2(x):
  '''
  Compute the square root y of x such that y*y = x, where x is a symmetic
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
  >>> error = norm(x-dot(y,y))
  >>> error < 1e-14
  True
  '''
  l=cholesky(x,lower=1)
  u,s,vt = svd(l)
  return dot(u,transpose(s*u))


def bdiag_to_mat(x):
  n = len(x)
  m = [ [None]*n for i in range(n) ]
  for i in range(n):
    for j in range(n):
      m[i][j] = diag(x[i][j])
  return bmat(m)


def block_cholesky(w):
  '''
  Efficiently compute the lower Cholesky decomposition of a block matrix
  composed entirely of diagonal matrices.  Input is given as a 2-dimensional
  list of diagonal vectors.  As the results are also in the form of a block
  matrix composed of diagonal vectors, output is also a 2-dimensional list
  of diagonal vectors.

  This is essentially the classical scalar Cholesky-Crout algorithm adapted
  to exploit the block structure of the input and the diagonal nature of
  each block element.

  The traditional algorithm costs n^3*m^3/3 FLOPS, where n is the number
  of row/column blocks and m is the diagonal length of each block.  This
  algorithm requires only n^3*m/3 FLOPS, a reduction over the usual
  algorithm of a factor of m^2 for suitably structured inputs.

  Example of a 2x2 block matrix of 2x2 diagonal matrices:

  >>> w11 = array([16.,36.])
  >>> w22 = array([64.,25.])
  >>> w12 = array([-4.,-16.])
  >>> w   = [[w11,w12],[w12,w22]]
  >>> x   = bdiag_to_mat(w)
  >>> x
  matrix([[ 16.,   0.,  -4.,   0.],
          [  0.,  36.,   0., -16.],
          [ -4.,   0.,  64.,   0.],
          [  0., -16.,   0.,  25.]])

  Perform a lower Cholesky decomposition on the expanded form:

  >>> l = matrix(cholesky(x.astype(float),lower=1))
  >>> l
  matrix([[ 4.        ,  0.        ,  0.        ,  0.        ],
          [ 0.        ,  6.        ,  0.        ,  0.        ],
          [-1.        ,  0.        ,  7.93725393,  0.        ],
          [ 0.        , -2.66666667,  0.        ,  4.22952585]])

  Verify x=l*l':

  >>> norm(dot(l,transpose(l))-x) < 1e-14
  True

  Now do the same thing using our algorithm, still using the packed format:

  >>> l2 = bdiag_to_mat(block_cholesky(w))
  >>> l2
  matrix([[ 4.        ,  0.        ,  0.        ,  0.        ],
          [ 0.        ,  6.        ,  0.        ,  0.        ],
          [-1.        ,  0.        ,  7.93725393,  0.        ],
          [ 0.        , -2.66666667,  0.        ,  4.22952585]])

  Verify the factorization matches the expanded version above and solves x=l*l':

  >>> norm(l2-l) < 1e-14
  True
  >>> norm(dot(l2,transpose(l2))-x) < 1e-14
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

  Elements of L can be defermined in a manner akin to the Cholesky-Crout
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
  c = [ [zeros(m)]*n for i in range(n) ]

  for j in range(n):
    for i in range(j,n):
      s = w[i][j] - sum(c[i][k]*c[j][k] for k in range(min(i,j)))
      if i==j:
        c[i][j] = s**0.5
      else:
        c[i][j] = s/c[j][j]

  return c


def grand_mean(X):
  return bmat([ ones((len(X),1)), X ])


def linreg(y, X, add_mean=False):
  '''
  Full-rank Example
  -----------------

  >>> X = matrix([[ 0.60, 1.20, 3.90],
  ...             [ 5.00, 4.00, 2.50],
  ...             [ 1.00,-4.00,-5.50],
  ...             [-1.00,-2.00,-6.50],
  ...             [-4.20,-8.40,-4.80]])
  >>> y = matrix([3, 4, -1, -5, -1]).T

  Least-squares Fit:

  >>> beta,cov,ss = linreg(y, X)

  Results:

  >>> beta.T
  matrix([[ 0.95333333, -0.84333333,  0.90666667]])
  >>> ss
  0.17
  >>> cov
  matrix([[ 0.06222222, -0.04222222,  0.01333333],
          [-0.04222222,  0.05444444, -0.02888889],
          [ 0.01333333, -0.02888889,  0.02666667]])

  Rank Deficient Example
  ----------------------

  >>> X = matrix([[0.05,  0.05, 0.25, -0.25],
  ...             [0.25,  0.25, 0.05, -0.05],
  ...             [0.35,  0.35, 1.75, -1.75],
  ...             [1.75,  1.75, 0.35, -0.35],
  ...             [0.30, -0.30, 0.30,  0.30],
  ...             [0.40, -0.40, 0.40,  0.40]])
  >>> y = matrix([1, 2, 3, 4, 5, 6]).T

  Least-squares Fit (or lack thereof):

  >>> beta,cov,ss = linreg(y,X)

  Results (not ideal due to extremely ill-conditioned design)

  Parameters estimates are well off
  >>> beta.T
  matrix([[ 10.125   ,  -8.59375 ,   0.421875,  -2.90625 ]])

  Residual variance highly inflated -- should be ~0.83
  >>> ss
  6.84427490234

  Huge oscillating covariances, crazy negative variances
  >>> cov
  matrix([[ -1.12589991e+15,   1.12589991e+15,   1.12589991e+15,
             1.12589991e+15],
          [  1.12589991e+15,  -1.12589991e+15,  -1.12589991e+15,
            -1.12589991e+15],
          [  1.12589991e+15,  -1.12589991e+15,  -1.12589991e+15,
            -1.12589991e+15],
          [  1.12589991e+15,  -1.12589991e+15,  -1.12589991e+15,
            -1.12589991e+15]])
  '''
  y = asmatrix(y)
  X = asmatrix(X)

  # Add type specific means, if requested
  if add_mean:
    X = grand_mean(X)

  # n=observations and m=covariates
  n,m = X.shape

  # Fit using normal equations
  # FIXME: not a good idea for numerical stability when posed with
  # rank-deficient or ill-conditioned designs
  W  = (X.T*X).I
  b  = W*X.T*y
  ss = ((y - X*b).A**2).sum()/(n-m)
  return b,W,ss


def logit(y, X, initial_beta=None, add_mean=False, max_iterations=50):
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
    z = y - mu1 + w*eta1.A
    b = W*X.T*z

  else:
    raise RuntimeError('logit estimator failed to converge')

  # Return log-likelihood, regression estimates, and covariance matrix
  return L,b,W


class GLogit(object):

  def __init__(self, y, X, ref=None, vars=None, add_mean=False, max_iterations=50):
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
        raise ValueError('reference class for dependent variable not observed in data')
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
    mu = array( [(y==i+1).sum() for i in range(k) ], dtype=float)/n

    # Initial estimates of beta (b1,b2)
    if initial_beta is None:
      b = []
      for i in range(k):
        bi = asmatrix(zeros(m), dtype=float).T
        bi[0,0] = log(mu[i]/(1-mu.sum()))
        b.append(bi)
    else:
      b = vsplit(initial_beta,k)

    # Initialize log likelihood to -infinity
    L = -inf

    # Iterate until convergence
    for iteration in range(max_iterations):
      eta = [ X*bi for bi in b ]

      # Compute expected values and linear predictors
      mu0  = 1/(1+sum(exp(etai) for etai in eta)).A
      eta0 = log(mu0)
      mu   = [ exp(eta[i]).A*mu0 for i in range(k) ]

      # Compute likelihood
      newL = eta0.sum() + y.choose( [0]+eta ).sum()

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
      raise RuntimeError('glogit estimator failed to converge')

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
    y = self.y_ord
    X = self.X
    k = len(mu)

    xz = []
    for i in range(k):
      # Form adjusted dependent variable: z_i = y_i-mu_i + eta
      z   = (y==i+1) - mu[i] + sum( w[i][j]*eta[j].A for j in range(k) )
      # Evaluate X'z_i
      xz += [[X.T*z]]

    # Obtain new betas: (X'*w*X*)^-1 * X'*w*z
    return W * bmat(xz)


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
      # that are to be scored and excluding those from the null model
      design_indices = tally(i%n for i in indices)
      assert design_indices
      assert all(m == k for m in design_indices.itervalues())

      null_indices = [ i for i in range(n) if i not in design_indices ]
      X_null = X[:,null_indices]

      # Fit null model
      self.null_L,beta_null,W = GLogit(y,X_null).fit(initial_beta=initial_beta)

      # Augment null beta with zeros for all parameters to be tested
      null_indices = (i for i in range(n*k) if i not in indices)
      beta = asmatrix(zeros((X.shape[1]*k,1), dtype=float))
      for b,i in izip(beta_null.A.ravel(),null_indices):
        beta[i,0] = b
    else:
      self.null_L = None

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

  def __init__(self, y, X, vars=None, add_mean=False):
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
      null = Linear(y,X_null)
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


def test_glogit():
  '''
  Example trichotomous outcome data taken from 'Logistic Regression: A
  Self-Learning Text', David G. Kleinbaum, Springer-Verlag Telos, ISBN
  978-0387941424 based on dataset from <need to find the book>.

  Load and parse test data

  >>> from pkg_resources import resource_stream
  >>> D = []
  >>> indices = [4,5,3,6]
  >>> for line in resource_stream(__name__,'test/cancer.dat'):
  ...   fields = line.split()
  ...   row = [ fields[i] for i in indices ]
  ...   if '.' in row:
  ...     continue
  ...   row = map(float,row)
  ...   D.append(row)

  >>> D=asmatrix(D,dtype=float)

  Extract dependent variable (first column) and design matrix (following columns)

  >>> y,X = D[:,0],D[:,1:]

  Compute the polytomous logit model for 3 outcomes

  >>> L,b,W = GLogit(y,X,ref=0).fit()

  Compute standard errors

  >>> stde = W.diagonal().A**.5

  Compute Z-statistics

  >>> z = b.T.A/stde

  Show results

  >>> print 'obs =',len(D)
  obs = 286
  >>> print 'logL =',L
  logL = -264.821392873
  >>> print '-2logL =',-2*L
  -2logL = 529.642785746
  >>> print 'beta =',b.T
  beta = [[-0.50851018 -1.41766704 -0.23349802 -0.56820429 -0.75411164 -2.43513171]]
  >>> print 'ss =',stde
  ss = [[ 0.22752862  0.29154207  0.47386178  0.22469103  0.25034798  1.03396066]]
  >>> print 'Z =',z
  Z = [[-2.2349284  -4.86265001 -0.49275555 -2.52882501 -3.01225372 -2.35514928]]
  >>> print 'X2 =',z**2
  X2 = [[  4.99490495  23.64536512   0.24280803   6.39495594   9.07367248
      5.54672813]]

  Test model against empty model

  >>> print 'score test =',GLogit(y,X).score_test().test()
  score test = (44.858794628046766, 4)
  '''


if __name__ == '__main__':
  import doctest
  doctest.testmod()
