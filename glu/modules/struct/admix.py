# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Estimate admixture proportions of a series of samples with a series of assumed ancestral populations'
__copyright__ = 'Copyright (c) 2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

import numpy as np
import scipy.optimize

from   itertools                 import izip

from   glu.lib.fileutils         import table_writer
from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.genoarray import genotype_count_matrix, genotype_indices


def log_like(genos,pops):
  k = len(pops)

  # Build indices into freuqnecy data from genotypes
  ind  = np.asarray(genotype_indices(genos),dtype=int)

  # Pre-generate genotypes frequencies
  mask = ind>0
  indices = np.arange(len(ind))
  pop_freq = [ pop[indices,ind][mask] for pop in pops ]

  def log_like_fn(mix):
    mix = np.asarray(mix)

    # Check bounds
    if mix.min() < -0.0 or mix.max() > 1:
      return np.inf

    s = mix.sum()
    if s-1>1e-6 or s<0:
      return np.inf

    if (~np.isfinite(mix)).any():
      return np.inf

    # Augment parameters
    if len(mix)+1==k and s<1:
      mix = mix.tolist() + [1-s]

    # Compute weighted mixture of likelihoods per locus
    l = sum(m*pop for m,pop in izip(mix,pop_freq))

    # Mask missing genotypes, take the natural log, and sum over loci
    l = -np.log(l).sum()

    # Return log-likelihood
    return l

  # Return likelihood function now that precomputations are done
  return log_like_fn


def progress_bar(samples, sample_count):
  try:
    from glu.lib.progressbar import progress_loop
  except ImportError:
    return samples

  update_interval = max(1,min(sample_count//100,250))

  return progress_loop(samples, length=sample_count, units='samples', update_interval=update_interval)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] test_genotypes pop1_genotypes pop2_genotypes [pop3_genotypes...]'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True)

  parser.add_option('--label', dest='labels', metavar='LABEL', action='append', default=[],
                    help='Population label (specify one per population)')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='output table file name')
  parser.add_option('-P', '--progress', dest='progress', action='store_true',
                    help='Show analysis progress bar, if possible')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) < 2:
    parser.print_help()
    sys.exit(2)

  k = len(args)-1

  # Build population labels
  labels = [ label.strip() for label in options.labels ]

  if '' in labels:
    raise ValueError('Blank population label specified')

  if len(labels) != len(set(labels)):
    raise ValueError('Duplicate population label specified')

  if len(labels) > len(args)-1:
    raise ValueError('Too many population labels specified')

  while len(labels) < k:
    labels.append('POP%d' % (len(labels)+1))

  # Load samples to test
  sys.stderr.write('Loading %s...\n' % args[0])
  test = load_genostream(args[0],format=options.informat,genorepr=options.ingenorepr,
                                 genome=options.loci,phenome=options.pedigree,
                                 transform=options, hyphen=sys.stdin).as_sdat()

  # Initialize masks
  loci = test.loci
  locusset = set(test.loci)

  # Load source populations, align loci, and compute frequencies
  pops = []
  for arg in args[1:]:
    sys.stderr.write('Loading %s...\n' % arg)
    genos = load_genostream(arg,format=options.informat,genorepr=options.ingenorepr,
                                genome=test.genome,phenome=options.pedigree,
                                transform=options,includeloci=locusset,orderloci=loci)

    # Count genotypes
    pop_loci,samples,geno_counts = genotype_count_matrix(genos)

    # Update masks
    locusset &= set(pop_loci)
    loci = [ l for l in loci if l in locusset ]
    mask = np.array([ l in locusset for l in pop_loci ],dtype=bool)
    geno_counts = geno_counts[mask]

    # Set each genotype to be observed at least once
    np.clip(geno_counts,1,1e300,out=geno_counts)

    # Set missing genotypes to zero
    geno_counts[:,0] = 0

    # Compute frequencies
    n = geno_counts.sum(axis=1)[:,np.newaxis]
    geno_freqs = geno_counts/n

    # Append to list of source populations
    pops.append( (pop_loci,geno_freqs) )

  # Perform final mask of individual data
  test = test.transformed(includeloci=loci)

  # Perform final mask of frequency data
  for i,(pop_loci,geno_freqs) in enumerate(pops):
    mask = np.array([ (l in locusset) for l in pop_loci ],dtype=bool)
    pops[i] = geno_freqs[mask]

  # Initialize grid search parameters (0..1 in 0.1 increments for k-1 parameters)
  # N.B.: Half of the values will result in domain errors, but those will be
  #       kicked out very quickly
  grid = (slice(0,1.01,0.1),)*(k-1)

  if options.progress and test.samples:
    test = progress_bar(test, len(test.samples))

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE']+labels)
  for sample,genos in test:
    # Build likelihood function
    f = log_like(genos,pops)

    # Perform a grid search
    mix = scipy.optimize.brute(f, grid, finish=None)

    # Refine by some rounds of the Nelder-Mead downhill simplex algorithm
    mix = scipy.optimize.fmin(f, mix, disp=0)

    # Augment admixture estimates and ensure they conform to bounds
    mix = np.asarray(mix.tolist()+[1-mix.sum()],dtype=float)
    np.clip(mix,0,1,out=mix)

    # Write output
    out.writerow( [sample] + ['%.4f' % abs(a) for a in mix] )


if __name__=='__main__':
  main()
