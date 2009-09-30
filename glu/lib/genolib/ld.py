# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = '''\
Fast linkage disequilibrium (LD) estimation for allelic
SNP genotype data'''
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


__all__ = []

from   math        import log
from   itertools   import izip
from   collections import defaultdict

from   glu.lib.genolib.genoarray import GenotypeLookupError,GENO_ARRAY_VERSION

epsilon = 10e-10


def count_haplotypes_native(genos1, genos2):
  '''
  Count the various haplotype combinations and return a vector containing:
     c11 - haplotype counts for allele 1 by allele 1
     c12 - haplotype counts for allele 1 by allele 2
     c21 - haplotype counts for allele 2 by allele 1
     c22 - haplotype counts for allele 2 by allele 2
     dh  - double heterozygote haplotypes (uninformative)

  >>> from glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,build_model
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*1400)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos1 = model.genotypes*350
  >>> genos2 = model.genotypes*200+model.genotypes[1:]*200

  >>> count_haplotypes_native(genos1,genos2)
  (600, 200, 200, 600, 250)

  >>> def encode(genos): return [ model[g] for g in genos ]
  >>> genos1 = encode([('A','A'),('A','B'),('A','B'),('B','B'),('B','B'),('B','B'),('A','A')])
  >>> genos2 = encode([('A','A'),('A','B'),('A','A'),('A','A'),('A','B'),('B','B'),('B','B')])
  >>> for i in range(len(genos1)+1):
  ...   count_haplotypes_native(genos1[:i],genos2[:i])
  (0, 0, 0, 0, 0)
  (2, 0, 0, 0, 0)
  (2, 0, 0, 0, 1)
  (3, 0, 1, 0, 1)
  (3, 0, 3, 0, 1)
  (3, 0, 4, 1, 1)
  (3, 0, 4, 3, 1)
  (3, 2, 4, 3, 1)

  # X-linked test

  >>> model = build_model('AB',allow_hemizygote=True)
  >>> genos1  = encode([(None,'A'),(None,'B'),(None,'A'),('A','A'),('A','B'),('B','B'),('A','A')])
  >>> genos2  = encode([(None,'A'),(None,'B'),(None,'A'),('A','A'),('A','B'),('B','B'),('B','B')])
  >>> count_haplotypes_native(genos1[:i],genos2[:i])
  Traceback (most recent call last):
     ...
  ValueError: Hemizygote LD estimation is not currently supported
  '''
  if len(genos1) != len(genos2):
    raise ValueError('genos1 and genos2 must be of same length')

  model1       = genos1[0].model if genos1 else None
  model2       = genos2[0].model if genos2 else None
  diplo_counts = count_diplotypes(genos1, genos2)
  het1,het2    = find_heterozygotes(model1,model2,diplo_counts)
  indices      = list( (2*i+j,i,j,(g1,g2)) for i,g1 in enumerate(het1) if g1 for j,g2 in enumerate(het2) if g2)
  x            = [0]*5

  for g1,g2,n in diplo_counts:
    if not g1 or not g2 or not n:
      continue

    if g1.hemizygote() or g2.hemizygote():
      raise ValueError('Hemizygote LD estimation is not currently supported')

    if g1.heterozygote() and g2.heterozygote():
      x[4] += n
      continue

    # Homozygotes count twice, since they appear in only one class --
    # conversely, all other configurations appear in two.
    if g1.homozygote() and g2.homozygote():
      n *= 2

    # Sum the counts of each category of allele configurations
    for i,a,b,dip in indices:
      if (g1[a],g2[b]) == dip:
        x[i] += n

  return tuple(x)


def count_diplotypes(genos1, genos2):
  '''
  Return a list of diplotype frequencies and a sets of alleles from each locus

  >>> from glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,build_model
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*1400)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos1 = model.genotypes*350
  >>> genos2 = model.genotypes*200+model.genotypes[1:]*200

  >>> for g1,g2,n in sorted(count_diplotypes(genos1,genos2)):
  ...   print g1,g2,n
  (None, None) (None, None) 200
  (None, None) ('A', 'A') 50
  (None, None) ('A', 'B') 50
  (None, None) ('B', 'B') 50
  ('A', 'A') ('A', 'A') 250
  ('A', 'A') ('A', 'B') 50
  ('A', 'A') ('B', 'B') 50
  ('A', 'B') ('A', 'A') 50
  ('A', 'B') ('A', 'B') 250
  ('A', 'B') ('B', 'B') 50
  ('B', 'B') ('A', 'A') 50
  ('B', 'B') ('A', 'B') 50
  ('B', 'B') ('B', 'B') 250

  >>> genos1 = GenotypeArray(descr,genos1)
  >>> genos2 = GenotypeArray(descr,genos2)
  >>> for g1,g2,n in sorted(count_diplotypes(genos1,genos2)):
  ...   print g1,g2,n
  (None, None) (None, None) 200
  (None, None) ('A', 'A') 50
  (None, None) ('A', 'B') 50
  (None, None) ('B', 'B') 50
  ('A', 'A') ('A', 'A') 250
  ('A', 'A') ('A', 'B') 50
  ('A', 'A') ('B', 'B') 50
  ('A', 'B') ('A', 'A') 50
  ('A', 'B') ('A', 'B') 250
  ('A', 'B') ('B', 'B') 50
  ('B', 'B') ('A', 'A') 50
  ('B', 'B') ('A', 'B') 50
  ('B', 'B') ('B', 'B') 250
  '''
  diplo_counts = defaultdict(int)
  for g1,g2 in izip(genos1,genos2):
    diplo_counts[g1,g2] += 1
  return tuple( (g1,g2,n) for (g1,g2),n in diplo_counts.iteritems() )


def find_heterozygotes(model1,model2,diplo_counts):
  '''Return exemplar heterozygotes for each locus'''
  a1 = set()
  a2 = set()
  for g1,g2,n in diplo_counts:
    a1.update(g1)
    a2.update(g2)

  return find_hetz(model1,a1),find_hetz(model2,a2)


def find_hetz(model,alleles):
  '''Return the heterozygote genotype for the given alleles'''
  alleles = set(alleles)
  alleles.discard(None)
  het = sorted(alleles)

  if len(het) > 2:
    raise ValueError('Only biallelic loci are allowed')

  while len(het) < 2:
    het.append(None)

  het = tuple(het)
  if model is not None:
    try:
      het = model[het]
    except GenotypeLookupError:
      pass

  return het


def estimate_ld_native(c11,c12,c21,c22,dh):
  '''
  Compute r-squared (pair-wise) measure of linkage disequilibrium for genotypes at two loci

          Haplotype Counts
          Locus1  Locus2
      c11    A       A
      c12    A       B
      c21    B       A
      c22    B       B

      dh     A       A    - double heterozygote haplotypes
             B       B
                or
             A       B
             B       A

  >>> import numpy
  >>> ld=estimate_ld_native(4,0,0,4,1)
  >>> numpy.allclose(ld, (1.,1.))
  True

  >>> ld=estimate_ld_native(1,1,0,6,1)
  >>> numpy.allclose(ld, (0.58333333, 1))
  True
  '''

  # Bail out on monomorphic markers
  information = (c11+c12, c21+c22, c11+c21, c12+c22)
  if not dh and 0 in information:
    return 0.,0.

  TOLERANCE = 10e-9

  # Initial estimate
  n = c11 + c12 + c21 + c22 + 2*dh
  p = float(c11 + c12 + dh)/n
  q = float(c11 + c21 + dh)/n

  p11 = p*q
  p12 = p*(1-q)
  p21 = (1-p)*q
  p22 = (1-p)*(1-q)

  loglike = -999999999

  for i in xrange(100):
    oldloglike=loglike

    # Force estimates away from boundaries
    p11=max(epsilon, p11)
    p12=max(epsilon, p12)
    p21=max(epsilon, p21)
    p22=max(epsilon, p22)

    a = p11*p22 + p12*p21

    loglike = c11*log(p11) + c12*log(p12) + c21*log(p21) + c22*log(p22) + dh*log(a)

    if abs(loglike-oldloglike) < TOLERANCE:
      break

    nx1 = dh*p11*p22/a
    nx2 = dh*p12*p21/a

    p11 = (c11+nx1)/n
    p12 = (c12+nx2)/n
    p21 = (c21+nx2)/n
    p22 = (c22+nx1)/n

  d = p11*p22 - p12*p21

  if d > 0:
    dmax =  min( p*(1-q), (1-p)*q )
  else:
    dmax = -min( p*q, (1-p)*(1-q) )

  dprime = d/dmax
  r2     = d*d/(p*(1-p)*q*(1-q))

  return r2,dprime


def bound_ld_native(c11,c12,c21,c22,dh):
  # Hack to estimate maxd and r2max
  n = c11 + c12 + c21 + c22 + 2*dh
  p = float(c11 + c12 + dh)/n
  q = float(c11 + c21 + dh)/n

  if p and p > 0.5:
    p = 1-p
    c11,c12,c21,c22 = c21,c22,c11,c12
  if q and q > 0.5:
    q = 1-q
    c11,c12,c21,c22 = c12,c11,c22,c21
  if p > q:
    p,q=q,p
    c11,c12,c21,c22 = c22,c21,c12,c11

  # Obtain rough estimate of d ignoring double-heterozygotes
  n -= 2*dh
  d  = float(c11*c22 - c12*c21)/n/n

  # Distinguish coupling from repulsion:
  #   Magic constant -0.005 can be refined by interval arithmetic, exploiting
  #   the minimum MAF and uncertainty in the estimates of p and q
  if d > -0.005:
    dmax =  min( p*(1-q), (1-p)*q )
  else:
    dmax = -min( p*q, (1-p)*(1-q) )

  if p > 0:
    r2max = dmax*dmax / (p*(1-p)*q*(1-q))
  else:
    r2max = 1.0

  return r2max


try:
  if GENO_ARRAY_VERSION != 'C':
    raise ImportError('Using Python version')

  # Load the optimized C versions, if available
  from glu.lib.genolib._genoarray import count_haplotypes as count_haplotypes_fast, \
                                         estimate_ld      as estimate_ld_fast


  count_haplotypes = count_haplotypes_fast
  estimate_ld      = estimate_ld_fast

  def test_count_haplotypes():
    '''
    >>> import numpy
    >>> from glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,build_model
    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*1400)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos1 = model.genotypes*350
    >>> genos2 = model.genotypes*200+model.genotypes[1:]*200

    >>> count_haplotypes_native(genos1,genos2)
    (600, 200, 200, 600, 250)

    >>> count_haplotypes_fast(genos1,genos2)
    (600, 200, 200, 600, 250)

    >>> def encode(genos): return GenotypeArray(GenotypeArrayDescriptor([model]*len(genos)), genos)
    >>> genos1 = encode([('A','A'),('A','B'),('A','B'),('B','B'),('B','B'),('B','B'),('A','A')])
    >>> genos2 = encode([('A','A'),('A','B'),('A','A'),('A','A'),('A','B'),('B','B'),('B','B')])
    >>> count_haplotypes_fast(genos1,genos2)
    (3, 2, 4, 3, 1)

    >>> for i in range(len(genos1)+1):
    ...   count_haplotypes_fast(genos1[:i],genos2[:i])
    (0, 0, 0, 0, 0)
    (2, 0, 0, 0, 0)
    (2, 0, 0, 0, 1)
    (3, 0, 1, 0, 1)
    (3, 0, 3, 0, 1)
    (3, 0, 4, 1, 1)
    (3, 0, 4, 3, 1)
    (3, 2, 4, 3, 1)

    # X-linked test

    >>> model = build_model('AB',allow_hemizygote=True)
    >>> genos1  = encode([(None,'A'),(None,'B'),(None,'A'),('A','A'),('A','B'),('B','B'),('A','A')])
    >>> genos2  = encode([(None,'A'),(None,'B'),(None,'A'),('A','A'),('A','B'),('B','B'),('B','B')])
    >>> count_haplotypes_fast(genos1[:i],genos2[:i])
    Traceback (most recent call last):
       ...
    ValueError: Hemizygote LD estimation is not currently supported
    '''

  def test_estimate_ld():
    '''
    >>> import numpy
    >>> ld=estimate_ld_fast(4,0,0,4,1)
    >>> numpy.allclose(ld, (1.,1.))
    True

    >>> ld=estimate_ld_fast(1,1,0,6,1)
    >>> numpy.allclose(ld, (0.58333333, 1))
    True
    '''

except ImportError:
  # If not, fall back on the pure-Python version
  estimate_ld      = estimate_ld_native
  count_haplotypes = count_haplotypes_native

bound_ld = bound_ld_native


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
