# -*- coding: utf-8 -*-
'''
File:          hwp.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


from scipy         import stats

from glu.lib.utils import izip_exact


HWP_EXACT_THRESHOLD=8000


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

  >>> hwp_exact_biallelic(3, 100, 5)
  6.1950944758140156e-21
  >>> hwp_exact_biallelic(14, 57, 50)
  0.84227975657079257
  >>> hwp_exact_biallelic(32, 31, 51)
  2.5637734657550698e-06
  >>> hwp_exact_biallelic(3, 47, 5)
  1.1629848615126043e-07
  >>> hwp_exact_biallelic(32, 150, 55)
  2.4853416055530676e-05
  >>> hwp_exact_biallelic(7, 122, 32)
  6.278449684759094e-13
  >>> hwp_exact_biallelic(3, 99, 14)
  4.3090717622326841e-16
  >>> hwp_exact_biallelic(13, 146, 54)
  4.0674031361063006e-10
  >>> hwp_exact_biallelic(100, 177, 57)
  0.18195180192910754
  >>> hwp_exact_biallelic(57, 184, 155)
  0.83102796343705576
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


def biallelic_counts(model,counts):
  '''
  Convert genotype counts (derived from genolib.genoarray.count_genotypes)
  into a tuple of homozygote1, heterozygote, and homozygote2 counts.
  '''
  if len(model.alleles) > 3:
    raise ValueError('Biallelic locus required')

  # Set to uninformative values
  homs = [0,0]
  hets = [0]

  for g,c in izip_exact(model.genotypes,counts):
    # Ignore hemizygotes and missing
    if g.homozygote():
      homs.append(c)
    elif g.heterozygote():
      hets.append(c)

  # Take final counts from the ends of each list
  hom1_count,hom2_count = min(homs[-2],homs[-1]),max(homs[-2],homs[-1])
  return hom1_count,hets[-1],hom2_count


def hwp_biallelic(model,counts):
  '''
  @param   genos: genotypes
  @type    genos: sequence
  @return       : asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes
  @rtype        : float
  '''
  return hwp_biallelic_counts(*hwp_biallelic_counts(model,counts))


def hwp_biallelic_counts(hom1_count,het_count,hom2_count,exact_threshold=None):
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

  >>> hwp_biallelic_counts(15,36,20,0)
  0.87188388159827424
  >>> hwp_biallelic_counts(15,36,10)
  0.19960651078273423
  >>> hwp_biallelic_counts(15,36,10,0)
  0.14135568108811056
  '''
  if exact_threshold is None:
    exact_threshold = HWP_EXACT_THRESHOLD

  # Only use the exact test when there are less than a fixed number of rare
  # alleles otherwise, use the asymptotic test
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
