# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Estimate genetic admixture proportions from series of assumed ancestral populations'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

import numpy as np

from   itertools                 import izip

from   glu.lib.fileutils         import table_writer, table_reader
from   glu.lib.progressbar       import progress_loop
from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.genoarray import genotype_count_matrix, genotype_indices, build_model, \
                                        GenotypeRepresentationError


# Set absolute and relative tolerances for solutions
EPSILON=np.finfo(float).eps
ABSTOL=1e-6
RELTOL=1e-9


def admixture_log_likelihood_python(f,x):
  '''
  >>> f = np.array([[0.25, 0.50, 0.25],
  ...               [0.50, 0.25, 1.00],
  ...               [0.50, 1.00, 1.00],
  ...               [0.25, 0.50, 0.50],
  ...               [0.75, 0.25, 0.25]])
  >>> x = np.array([0.25, 0.50, 0.25])

  >>> np.allclose(admixture_log_likelihood_python(f,x), -3.6150156523923)
  True
  '''
  return np.log(np.dot(f,x)).sum()


def admixture_log_likelihood_derivative_python(f,x):
  '''
  >>> f = np.array([[0.25, 0.50, 0.25],
  ...               [0.50, 0.25, 1.00],
  ...               [0.50, 1.00, 1.00],
  ...               [0.25, 0.50, 0.50],
  ...               [0.75, 0.25, 0.25]])
  >>> x = np.array([0.25, 0.50, 0.25])

  >>> d=admixture_log_likelihood_derivative_python(f,x)
  >>> np.allclose(d, [ 4.80952381,  4.78571429,  5.61904762])
  True
  '''
  return (f/np.dot(f,x)[:,np.newaxis]).sum(axis=0)


def individual_frequencies_python(populations,ind):
  mask    = ind>0
  indices = np.where(mask)[0]
  ind     = ind[mask]
  f       = populations[:,indices,ind].T
  return f


try:
  from glu.modules.struct._admix import individual_frequencies as individual_frequencies_c, \
                                        admixture_log_likelihood as admixture_log_likelihood_c, \
                                        admixture_log_likelihood_derivative as admixture_log_likelihood_derivative_c

  individual_frequencies = individual_frequencies_c
  admixture_log_likelihood = admixture_log_likelihood_c
  admixture_log_likelihood_derivative = admixture_log_likelihood_derivative_c

  def test_admixture_log_likelihood():
    '''
    >>> f = np.array([[0.25, 0.50, 0.25],
    ...               [0.50, 0.25, 1.00],
    ...               [0.50, 1.00, 1.00],
    ...               [0.25, 0.50, 0.50],
    ...               [0.75, 0.25, 0.25]])
    >>> x = np.array([0.25, 0.50, 0.25])

    >>> np.allclose(admixture_log_likelihood_c(f,x), -3.6150156523923)
    True

    >>> d=admixture_log_likelihood_derivative_c(f,x)
    >>> np.allclose(d, [ 4.80952381,  4.78571429,  5.61904762])
    True
    '''

except ImportError:
  individual_frequencies = individual_frequencies_python
  admixture_log_likelihood = admixture_log_likelihood_python
  admixture_log_likelihood_derivative = admixture_log_likelihood_derivative_python


def estimate_admixture_em(f,x0=None,iters=100):
  '''
  Problem: Maximize a likelihood to determine mixing proportions a series of
  events from a series of k Bernoulli distributions (a simplification of a
  binomial mixture problem)

  Let F be a (n x k) matrix of known frequencies of n events from k
  distributions.  We wish to maximize the following log-likelihood function
  over a (k x 1) vector x of proportions:

    lnL(x) = sum(ln(F*x))

    subject to sum(x) <= 1
               min(x) >= 0

  Note: n is typically 10,000-20,000
        k is typically 2..5

  The classical EM algorithm is used to estimate x.  It is _extremely_ slow
  to converge and is better suited for finding the neighborhood of the
  solution, rather than for precisely estimating the coefficients.  Thus we
  iterate a fixed number of times and do not bother checking for
  convergence.  Those estimates can then be used as a feasible starting
  point for algorithms with better local convergence properties.
  '''
  n,k = f.shape

  if x0 is None:
    x = np.ones(k)/k
  else:
    # Ensure we copy, since iterations update x in place
    x = np.array(x0)

  for i in xrange(iters):
    x *= admixture_log_likelihood_derivative(f,x)/n

  return x


def estimate_admixture_cvxopt(f, x0, maxiters=25, failover=True):
  '''
  Problem: Maximize a likelihood to determine mixing proportions a series of
  events from a series of k Bernoulli distributions (a simplification of a
  binomial mixture problem)

  Let F be a (n x k) matrix of known frequencies of n events from k
  distributions.  We wish to maximize the following log-likelihood function
  over a (k x 1) vector x of proportions:

    lnL(x) = sum(ln(F*x))

    subject to sum(x) <= 1
               min(x) >= 0

  Note: n is typically 10,000-20,000
        k is typically 2..5

  The CVXOPT NLP solver is used.  It is based on an interior point method
  that solves both the primary constrained problem and its dual by solving
  the KKT system of equations at every iteration.  The algorithm uses
  analytical first and second derivatives and is implemented in mostly
  Python code, so it is also relatively slow.

  Note: The solver will attempt to evaluate the likelihood and its
  derivatives at infeasible points and for the same values repeatedly.
  Also, it uses a custom matrix type that is generally less useful than the
  NumPy versions.
  '''
  # WARNING: This code mixes NumPy and CVXOPT data types.  Proceed with caution.
  from cvxopt import solvers, matrix

  n,k = f.shape
  x0  = matrix(x0)

  # Precompute cross-products of population frequencies (columns of f) for the Hessian
  ff = [ [ (f[:,i]*f[:,j])[:,np.newaxis] if i>=j else 0 for j in range(k) ]
                                                        for i in range(k) ]

  # Store last function values, since the solver seems to want to
  # re-evaluate them several times
  last   = []
  iters  = [0]

  # Build likelihood function
  def lnL(x=None, z=None):
    # Return number of non-linear constraints and initial parameters
    if x is None:
      return 0,x0

    # Do not check constraints or else optimization will often get "stuck"
    x = np.asarray(x, dtype=float).reshape(-1)

    # Check to see if we've just solved this case
    if last:
      d = (abs(x-last[0])).sum()
      if d==0:
        last_x,last_z,last_l,last_df,last_h = last
        if z is None:
          return last_l,last_df
        elif z[0]==last_z:
          return last_l,last_df,last_h

    iters[0] += 1

    # Compute mixture probabilities and log-likelihood
    l  = -admixture_log_likelihood(f,x)

    if not np.isfinite(l):
      return None

    # Compute derivatives
    df = -admixture_log_likelihood_derivative(f,x)
    df = matrix(df, tc='d').T

    if z is None:
      last[:] = [x+0,None,l,df,None]
      return l,df

    # Compute Hessian, if requested
    u2 = (np.dot(f,x)**2)[:,np.newaxis]
    h  = [ [ z[0]*(ff[i][j]/u2).sum() if i>=j else 0 for i in range(k) ]
                                                     for j in range(k) ]
    h  = matrix(h, tc='d')

    last[:] = [x+0,z[0],l,df,h]
    return l,df,h

  # Set up constraint matrices
  #   k inequality constraints for x[i]>=0
  G = matrix([ [0.]*i + [-1.] + [0]*(k-i-1) for i in range(k) ]).T
  h = matrix([0.]*k)

  #   1 equality constraint for sum(x)==1
  A = matrix(np.ones(k)).T
  b = matrix(np.ones(1))

  # Set solver options
  solvers.options['show_progress'] = False
  solvers.options['maxiters']      = maxiters
  solvers.options['abstol']        = ABSTOL
  solvers.options['reltol']        = RELTOL

  # Run solver
  sol = solvers.cp(lnL, G, h, A=A, b=b)

  # Return results (parameters, number of iterations, final log-likelihood)
  x = np.asarray(sol['x'])
  l = lnL(x)[0]

  # Allow algorithm to fail and re-try problem with SQP (without
  # recursive failover)
  if sol['status'] != 'optimal' and failover:
    x2,l2,iters2 = estimate_admixture_sqp(f, x0, failover=False)
    iters[0] += iters2

    if l2<l:
      x,l=x2,l2

  return x,l,iters[0]


def estimate_admixture_sqp(f, x0, failover=True):
  '''
  Problem: Maximize a likelihood to determine mixing proportions a series of events from a
  series of k Bernoulli distributions (a simplification of a binomial mixture problem)

  Let F be a (n x k) matrix of known frequencies of n events from k
  distributions.  We wish to maximize the following log-likelihood function
  over a (k x 1) vector x of proportions:

    lnL(x) = sum(ln(F*x))

    subject to sum(x) <= 1
               min(x) >= 0

  Note: n is typically 10,000-20,000
        k is typically 2..5

  Estimation is by the Sequential Quadratic Programming algorithm as
  implemented in SciPy, which is based on the Sequential Least SQuares
  Programming optimization algorithm (SLSQP), originally developed by Dieter
  Kraft.  See http://www.netlib.org/toms/733
  '''
  from scipy.optimize import fmin_slsqp

  n,k = f.shape
  x0 = np.array(x0)

  # Store last function values, since the solver seems to want to
  # re-evaluate them several times
  iters = [0]

  # Build likelihood function
  def lnL(x):
    # Do not check constraints or else optimization will often get "stuck"
    iters[0] += 1

    l = -admixture_log_likelihood(f,x)

    if not np.isfinite(l):
      return np.inf

    return l

  def dlnL(x):
    return -admixture_log_likelihood_derivative(f,x)

  # Define equality constraints and derivatives
  # Weight norms since that seems to push estimates away from bounds (should
  # be much larger than the norm of the log likelihood)
  w = 1e4**k
  eq = np.ones(1)
  def eqcons(x):
    eq[0] = w*(x.sum()-1)
    return eq

  deq = np.ones( (1,k) )*w
  def fprime_eqcons(x):
    return deq

  x,fx,its,imode,smode = fmin_slsqp(lnL, fprime=dlnL, x0=x0, bounds=[(0,1)]*k,
                                    f_eqcons=eqcons, fprime_eqcons=fprime_eqcons,
                                    iter=25, full_output=True, iprint=-1, acc=ABSTOL/10)

  # Sometimes x is returned as a list...
  x = np.asarray(x)

  # Allow algorithm to fail and re-try problem with CVXOPT (without
  # recursive failover)
  if imode != 0 and failover:
    # Unless iteration limit was exceeded, start from the initial parameters
    x1 = x if imode == 9 else x0

    x2,fx2,iters2 = estimate_admixture_cvxopt(f, x1, failover=False)
    iters[0] += iters2

    if fx2<fx:
      x=x2

  x = np.clip(x,0,1)
  return x,fx,iters[0]


def classify_ancestry(labels,x,threshold):
  '''
  An individual is considered of a given ancestry based on the supplied
  labels and estimated admixture coefficients if their coefficient is
  greater than a given threshold.

  Otherwise, an individual who has no single estimated admixture coefficient
  that meets the specified threshold then one of two behaviors result.  If
  only one population group exceeds 1-threshold then the ancestry is deemed
  'ADMIXED' for that population.  Otherwise, a list of populations with
  estimated admixture above 1-threshold is returned.
  '''
  popset = set()

  cmax = -1
  for pop,coeff in izip(labels,x):
    if coeff >= 1-threshold:
      popset.add(pop)
      cmax = max(cmax,coeff)

  if len(popset)==1 and cmax < threshold:
    ipop = 'ADMIXED %s' % popset.pop()
  else:
    ipop = ','.join(sorted(popset))

  return ipop


def compute_frequencies(freq_model,sample_count,models,geno_counts):
  # Set missing genotypes to zero
  geno_counts[:,0] = 0

  if freq_model.upper() == 'GENO':
    # Set each genotype to be observed at least once
    np.clip(geno_counts,1,1e300,out=geno_counts)
    geno_counts[:,0] = 0

    # Compute frequencies
    n = geno_counts.sum(axis=1)[:,np.newaxis]
    geno_freqs = geno_counts/n

  elif freq_model.upper() == 'HWP':
    geno_freqs = np.zeros(geno_counts.shape, dtype=float)

    for i,model in enumerate(models):
      n = 2*geno_counts[i].sum()

      if not n or len(model.alleles)!=3:
        p = 1/sample_count
      else:
        a,b  =  model.alleles[1:3]
        inds = (model[a,a].index,
                model[a,b].index,
                model[b,b].index)

        hom1 = geno_counts[i,inds[0]]
        hets = geno_counts[i,inds[1]]

        p    = (2*hom1+hets)/n

      q = 1-p

      geno_freqs[i,inds[0]] =   p*p
      geno_freqs[i,inds[1]] = 2*p*q
      geno_freqs[i,inds[2]] =   q*q
  else:
    raise ValueError('Invalid genotype likelihood model specified: %s' % model)

  return geno_freqs


def load_from_genos(options):
  # Load samples to test
  sys.stderr.write('Loading %s...\n' % options.test_genotypes)
  test = load_genostream(options.test_genotypes,format=options.informat,genorepr=options.ingenorepr,
                         genome=options.loci,phenome=options.pedigree,transform=options,
                         hyphen=sys.stdin).as_sdat()

  # Initialize masks
  loci = test.loci
  locusset = set(test.loci)

  # Load source populations, align loci, and compute frequencies
  pops = []
  for filename in options.pop_genotypes:
    sys.stderr.write('Loading %s...\n' % filename)
    genos = load_genostream(filename,format=options.informat,genorepr=options.ingenorepr,
                                genome=test.genome,phenome=options.pedigree,
                                transform=options,includeloci=locusset,orderloci=loci)

    # Count genotypes
    pop_loci,samples,geno_counts = genotype_count_matrix(genos)

    # Update masks
    locusset &= set(pop_loci)
    loci = [ l for l in loci if l in locusset ]
    mask = np.array([ l in locusset for l in pop_loci ],dtype=bool)

    geno_counts = geno_counts[mask]
    geno_freqs  = compute_frequencies(options.model,len(genos.samples),genos.models,geno_counts)

    # Append to list of source populations
    pops.append( (pop_loci,geno_freqs) )

  # Perform final mask of individual data
  test = test.transformed(includeloci=loci)

  # Perform final mask of frequency data
  for i,(pop_loci,geno_freqs) in enumerate(pops):
    mask = np.array([ (l in locusset) for l in pop_loci ], dtype=bool)
    pops[i] = geno_freqs[mask]

  pops = np.array(pops, dtype=float)

  return test,pops


RC={'A':'T','C':'G','T':'A','G':'C'}


def load_from_freqs(options):
  # Load samples to test
  sys.stderr.write('Loading %s...\n' % options.test_genotypes)
  test = load_genostream(options.test_genotypes,format=options.informat,genorepr=options.ingenorepr,
                         genome=options.loci,phenome=options.pedigree,transform=options,
                         hyphen=sys.stdin).as_sdat()

  # Initialize masks
  locusmap = dict( (l,i) for i,l in enumerate(test.loci) )

  # Load source populations, align loci, and compute frequencies
  pops = []
  sys.stderr.write('Loading %s...\n' % options.pop_genotypes[0])
  fdata = table_reader(options.pop_genotypes[0],want_header=True)

  header = next(fdata)
  pops   = header[2:]

  print '!!! pops=',pops
  print '!!! starting with %d loci' % len(locusmap)

  n      = len(pops)

  seen    = set()
  freqmap = {}
  for row in fdata:
    lname = row[0]

    if lname not in locusmap:
      continue

    if lname in seen:
      freqmap.pop(lname,None)
      continue

    seen.add(lname)

    alleles = tuple(row[1].split(','))

    if len(alleles)!=2:
      continue

    model = test.models[ locusmap[lname ] ]

    if alleles not in model:
      try:
        model = build_model(alleles,base=model)
      except GenotypeRepresentationError:
        alleles = RC[alleles[0]],RC[alleles[1]]
        if alleles not in model:
          try:
            model = build_model(alleles,base=model)
          except GenotypeRepresentationError:
            print '!!! incompatible model: alleles=%s, model=%s' % (alleles,model.alleles[1:])
            continue

    freqmap[lname] = freqs = np.zeros( (n,4), dtype=float )

    q   = np.clip(np.array(row[2:],dtype=float),0.005,1-0.005)
    p   = 1-q

    a,b = alleles
    aa  = model[a,a].index
    ab  = model[a,b].index
    bb  = model[b,b].index

    freqs[:,aa] =   p*p
    freqs[:,ab] = 2*p*q
    freqs[:,bb] =   q*q

  # Perform final mask of individual data
  loci = [ l for l in test.loci if l in freqmap ]

  print '!!! found %d loci' % len(freqmap)

  test = test.transformed(includeloci=loci)

  locusmap  = dict( (l,i) for i,l in enumerate(loci) )
  pop_freqs = np.zeros( (n,len(loci),4), dtype=float)

  for i,lname in enumerate(loci):
    pop_freqs[:,i,:] = freqmap.pop(lname)

  print '!!! found %d loci' % len(loci)

  return test,pop_freqs


def build_labels(options):
  k = len(options.pop_genotypes)

  labels = []
  for label in options.labels:
    labels.extend( l.strip() for l in label.split(',') )

  if '' in labels:
    raise ValueError('Blank population label specified')

  if len(labels) != len(set(labels)):
    raise ValueError('Duplicate population label specified')

  if len(labels) > k:
    raise ValueError('Too many population labels specified')

  while len(labels) < k:
    labels.append('POP%d' % (len(labels)+1))

  return labels


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('test_genotypes',            help='Input genotypes to test')
  parser.add_argument('pop_genotypes',  nargs='+', help='Population reference genotypes')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('--labels', metavar='LABELS', action='append', default=[],
                    help='Population labels (specify one per population separated with commas)')
  parser.add_argument('--model', metavar='MODEL', default='HWP',
                    help='Model for genotype frequencies.  HWP to assume Hardy-Weinberg proportions, '
                         'otherwise GENO to fit genotypes based on frequency.  (Default=HWP)')
  parser.add_argument('-t', '--threshold', metavar='N', type=float, default=0.80,
                    help='Imputed ancestry threshold (default=0.80)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='output table file name')
  parser.add_argument('-P', '--progress', action='store_true',
                    help='Show analysis progress bar, if possible')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()

  labels    = build_labels(options)

  if len(options.pop_genotypes)>1:
    test,pops = load_from_genos(options)
  else:
    test,pops = load_from_freqs(options)

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE']+labels+['IMPUTED_ANCESTRY'])

  if options.progress and test.samples:
    test = progress_loop(test, length=len(test.samples), units='samples')

  for sample,genos in test:
    # Compute genotype frequencies
    f      = individual_frequencies(pops,genotype_indices(genos))

    # Find feasible starting values
    x0     = estimate_admixture_em(f,iters=10)

    # Estimate admixture
    x,l,it = estimate_admixture_sqp(f, x0)

    ipop   = classify_ancestry(labels, x, options.threshold)

    out.writerow([sample]+['%.4f' % a for a in x] + [ipop])


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  main()
