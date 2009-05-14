# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Principle Components Analysis (PCA) to find large-scale correlations among samples based on genotype data'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   itertools                 import izip,islice

import numpy

try:
  import cvxopt, cvxopt.blas, cvxopt.lapack
except ImportError:
  cvxopt = None

from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.genoarray import count_genotypes, count_alleles_from_genocounts, \
                                        major_allele_from_allelecounts


def encode_loci(genos):
  genos = genos.as_ldat()

  for (locus,geno),model in izip(genos,genos.models):
    genocounts   = count_genotypes(geno)
    allelecounts = count_alleles_from_genocounts(model,genocounts)

    try:
      major,freq = major_allele_from_allelecounts(model,allelecounts)
    except ValueError:
      major = None

    # Uninformative
    if major is None:
      continue

    other = [ a for a,n in izip(model.alleles[1:],allelecounts[1:]) if a!=major and n ]

    # Monomorphic or non-bialleleic
    if len(other) != 1:
      continue

    other = other[0]

    # Form numeric encoding of informative data (0,.5,1)
    # NB: dicts are faster than lists here-- sigh.
    gmap = { model[major,major]:0.0,
             model[other,major]:0.5,
             model[other,other]:1.0 }
    inf  = [ gmap[g] for g in geno if g ]

    # Mean center and normalize scale
    s,n  = sum(inf),len(inf)
    avg  = s/n
    p    = (s+0.5)/(n+1)
    norm = numpy.sqrt(p*(1-p))

    # Form genotype to score mapping
    # NB: dicts are faster than lists here-- sigh.
    gmap = { model[None,None]  :  0.0,
             model[major,major]: (0.0-avg)/norm,
             model[other,major]: (0.5-avg)/norm,
             model[other,other]: (1.0-avg)/norm }

    yield [ gmap[g] for g in geno ]


def tracy_widom(evals):
  evals = numpy.asarray(evals).reshape(-1)
  l1 = evals[::-1].cumsum()[::-1]
  l2 = (evals[::-1]**2).cumsum()[::-1]
  n  = numpy.arange(len(evals),0,-1,dtype=float)

  oldf = numpy.seterr(all='ignore')

  try:
    s2 = n**2*l2/(l1**2)
    v  = n*(n+2)/(s2-n)
    l  = n*evals/l1
    vs = numpy.sqrt(v-1)
    ns = numpy.sqrt(n)
    mu = (vs+ns)**2/v
    ss = (vs+ns)/v*(1/vs+1/ns)**(1./3.)
    tw = (l-mu)/ss
  finally:
    numpy.seterr(**oldf)

  return tw,v


def pca_cov_numpy(cov,n=None):
  '''
  Let
     n be the number of subjects
     m be the number of loci
     data is an nxm matrix
     evec is a nx1 vector of eigenvalues
     eval is a nxm matrix of eigenvectors

  The time complexity of PCA using SVD is O(n*m**2)
  EIG:
  '''
  cov = numpy.asarray(cov)

  shape = cov.shape
  if len(shape) != 2 or shape[0] != shape[1]:
    raise ValueError('Invalid covariance matrix for PCA')

  values,vectors = numpy.linalg.eigh(cov)
  order = numpy.flipud(values.argsort())
  values = numpy.where(values < 0, 0, values)

  return values[order],vectors[:,order].T


def pca_sequence_numpy(genos,n=None,chunksize=500):
  '''
  Let
     n be the number of subjects
     m be the number of loci
     data is an nxm matrix
     evec is a nx1 vector of eigenvalues
     eval is a nxm matrix of eigenvectors

  The time complexity of PCA using symmetric eigenvalue decomposition is O(m*n**2)
  '''
  s,m = len(genos.samples),0
  cov = numpy.zeros( (s,s), dtype=float )

  data = encode_loci(genos)

  while 1:
    chunk = list(islice(data,chunksize))

    if not chunk:
      break

    m    += len(chunk)
    chunk = numpy.array(chunk, dtype=float)
    cov  += numpy.dot(chunk.T, chunk)

  return pca_cov_numpy(cov/m,n)


def pca_cov_cvxopt(cov,n=None):
  '''
  Let
     n be the number of subjects
     m be the number of loci
     data is an nxm matrix
     evec is a nx1 vector of eigenvalues
     eval is a nxm matrix of eigenvectors

  The time complexity of PCA using SVD is O(n*m**2)
  EIG:
  '''
  m = cov.size[0]
  n = n or m

  # Create result matrices
  values  = cvxopt.base.matrix(0.0, (1,n) )
  vectors = cvxopt.base.matrix(0.0, (m,n) )

  # Compute the top m eigenvalues and vectors
  cvxopt.lapack.syevr(cov,values,Z=vectors,jobz='V',
                      range='I',il=max(1,m-n+1),iu=m)

  # syevr returns eigenvalues and vectors in ascending order,
  # so reverse them
  values  = values[:,::-1]
  vectors = vectors[:,::-1].T

  return values,vectors


def pca_sequence_cvxopt(genos,n=None,chunksize=500):
  '''
  Let
     n be the number of subjects
     m be the number of loci
     data is an nxm matrix
     evec is a nx1 vector of eigenvalues
     eval is a nxm matrix of eigenvectors

  The time complexity of PCA using symmetric eigenvalue decomposition is O(m*n**2)
  '''
  s,m = len(genos.samples),0
  cov = cvxopt.base.matrix(0.0, (s,s))

  data = encode_loci(genos)

  # Read blocks of up to chunksize
  while 1:
    chunk = list(islice(data,chunksize))

    if not chunk:
      break

    m += len(chunk)

    # Compute cov += dot(C.T,C)
    chunk = cvxopt.base.matrix(chunk)
    cvxopt.blas.syrk(chunk,cov,beta=1.0)

  return pca_cov_cvxopt(cov/m,n)


def do_pca(genos,n=None):
  genos = genos.as_ldat()

  if cvxopt:
    return pca_sequence_cvxopt(genos,n)
  else:
    return pca_sequence_numpy(genos,n)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotypes'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output principle components (eigenvectors) to FILE (default is "-" for standard out)')
  parser.add_option('-O', '--vecout', dest='vecout', metavar='FILE',
                    help='Output eigenvalues and statistics to FILE')

  parser.add_option('--vectors', dest='vectors', metavar='N', type='int', default=10,
                    help='Output the top N eigenvectors.  Default=10')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  genos = load_genostream(args[0],format=options.informat,genorepr=options.ingenorepr,
                                  genome=options.loci,phenome=options.pedigree,
                                  transform=options).as_ldat()

  n = options.vectors or None

  # Compute on all N or else we can't get TW_STATS if options.vecout
  m = n if not options.vecout else None
  values,vectors = do_pca(genos,n=m)

  if options.vecout:
    # Compute TW Stats before truncating data
    twstats,n_eff=tracy_widom(values)

  if n:
    values  = values[:n]
    vectors = vectors[:n,:]
    if options.vecout:
      twstats = twstats[:n]
      n_eff   = n_eff[:n]

  if options.output:
    out = table_writer(options.output,hyphen=sys.stdout)
    out.writerow(['ID'] + [ 'EV%d' % (i+1) for i in range(len(values)) ])
    for i,sample in enumerate(genos.samples):
      out.writerow( [sample]+['%7.4f' % v for v in vectors[:,i]] )

  if options.vecout:
    out = table_writer(options.vecout)
    out.writerow(['N', 'EV','TW_STAT','N_EFF'])
    for i in xrange(len(values)):
      out.writerow(['%d' % (i+1), '%.4f' % values[i],
                                  '%.4f' % twstats[i],
                                  '%.2f' % n_eff[i]])


if __name__ == '__main__':
  main()
