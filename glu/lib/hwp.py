# -*- coding: utf-8 -*-
'''
File:          hwp.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from scipy         import stats

from glu.lib.utils import tally


def count_genos(genos):
  '''
  Count non-hemizygous non-missing genotypes

  @param   genos: genotypes
  @type    genos: sequence
  @return       : Count of observed genotypes for homozygote 1,
                  heterozygotes, and homozygote 2
  @rtype        : tuple
  '''

  f = tally(genos)

  hom1 = hom2 = het = 0

  for g,n in f.iteritems():
    if not g:
      continue

    a1,a2 = g

    if not a1 or not a2:
      continue
    elif a1!=a2:
      het = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)
  return hom1,het,hom2


def hwp_exact_biallelic(hom1_count, het_count, hom2_count):
  '''
  Exact SNP test for deviations from Hardy-Weinberg proportions.  Based on
  'A Note on Exact Tests of Hardy-Weinberg Equilibrium', Wigginton JE,
  Cutler DJ and Abecasis GR; Am J Hum Genet (2005) 76: 887-93

  @param  hom1_count: Count of observed homogygote 1
  @type   hom1_count: int
  @param   het_count: Count of observed heterogygote
  @type    het_count: int
  @param  hom1_count: Count of observed homogygote 2
  @type   hom2_count: int
  @return           : Exact p-value for deviation (2-sided) from Hardy-Weinberg Proportions (HWP)
  @rtype            : float
  '''

  # Computer the number of rare and common alleles
  rare   = 2*min(hom1_count,hom2_count)+het_count
  common = 2*max(hom1_count,hom2_count)+het_count

  if not rare:
    return 1.

  # Compute the expected number of heterogygotes under HWP
  hets = rare*common/(rare+common)

  # Account for rounding error on the number of hets, if the
  # parity of rare and hets do not match
  if rare%2 != hets%2:
    hets += 1

  # Initialize the expected number of rare and common homogygotes under HWP
  hom_r = (rare-hets)/2
  hom_c = (common-hets)/2

  # Initialize heterozygote probability vector, such that once filled in
  # P(hets|observed counts) = probs[hets/2]/sum(probs)
  probs = [0]*(rare/2+1)

  # Set P(expected hets)=1, since the remaining probabilities will be
  # computed relative to it
  probs[hets/2] = 1.0

  # Fill in relative probabilities for less than the expected hets
  for i,h in enumerate(xrange(hets,1,-2)):
    probs[h/2-1] = probs[h/2]*h*(h-1) / (4*(hom_r+i+1)*(hom_c+i+1))

  # Fill in relative probabilities for greater than the expected hets
  for i,h in enumerate(xrange(hets,rare-1,2)):
    probs[h/2+1] = probs[h/2]*4*(hom_r-i)*(hom_c-i) / ((h+1)*(h+2))

  # Compute the pvalue by summing the probabilities <= to that of the
  # observed number of heterogygotes and normalize by the total
  p_obs = probs[het_count/2]
  pvalue = sum(p for p in probs if p <= p_obs)/sum(probs)

  return pvalue


def hwp_chisq_biallelic(hom1_count, het_count, hom2_count):
  '''
  Return the asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes

  @param  hom1_count: Count of observed homogygote 1
  @type   hom1_count: int
  @param   het_count: Count of observed heterogygote
  @type    het_count: int
  @param  hom1_count: Count of observed homogygote 2
  @type   hom2_count: int
  @return           : asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes
  @rtype            : float

  >>> hom1_count, het_count, hom2_count = 15,36,20
  >>> hwp_chisq_biallelic(hom1_count, het_count, hom2_count)
  0.87188388159827424
  '''
  n = hom1_count + het_count + hom2_count

  if not n:
    return 1.0

  p = float(2*hom1_count+het_count)/(2*n)
  q = float(2*hom2_count+het_count)/(2*n)

  def score(o,e):
    return (o-e)**2/e if e>0 else 0.

  xx = (score(hom1_count,   n*p*p)
     +  score( het_count, 2*n*p*q)
     +  score(hom2_count,   n*q*q))

  return float(stats.distributions.chi2.sf(xx,1))


def hwp_biallelic(genos):
  '''
  @param   genos: genotypes
  @type    genos: sequence
  @return       : asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes
  @rtype        : float

  '''
  return hwp_biallelic_counts(*count_genos(genos))


def hwp_biallelic_counts(hom1_count,het_count,hom2_count,exact_threshold=8000):
  '''
  Return the asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes

  @param       hom1_count: Count of observed homogygote 1
  @type        hom1_count: int
  @param        het_count: Count of observed heterogygote
  @type         het_count: int
  @param       hom1_count: Count of observed homogygote 2
  @type        hom2_count: int
  @param  exact_threshold: threshold for the exact test
  @type   exact_threshold: int
  @return                : asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes
  @rtype                 : float

  >>> hom1_count, het_count, hom2_count = 15,36,20
  >>> hwp_chisq_biallelic(hom1_count, het_count, hom2_count)
  0.87188388159827424
  '''
  # Only use the exact test when there are less than 8000 rare alleles
  # otherwise, use the asymptotic test
  if 2*min(hom1_count,hom2_count)+het_count < exact_threshold:
    p = hwp_exact_biallelic(hom1_count, het_count, hom2_count)
  else:
    p = hwp_chisq_biallelic(hom1_count, het_count, hom2_count)

  return p


def _test():
  import doctest
  return doctest.testmod()


if __name__=='__main__':
  _test()
