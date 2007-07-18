from scipy import stats

from utils import tally


def count_genos(genos):
  '''
  Estimate allele and genotype frequencies
  Missing alleles are coded as ' '
  '''
  f = tally(g for g in genos if g and ' ' not in g)

  hom1 = hom2 = het = 0

  for g,n in f.iteritems():
    if g[0] != g[1]:
      het = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)
  return hom1,het,hom2


def hwp_exact_biallelic(hom1_count, het_count, hom2_count):
  '''
  hwp_exact_biallelic(count, het_count, hom2_count):

  Exact SNP test for deviations from Hardy-Weinberg proportions.  Based on
  'A Note on Exact Tests of Hardy-Weinberg Equilibrium', Wigginton JE,
  Cutler DJ and Abecasis GR; Am J Hum Genet (2005) 76: 887-93

  Input: Count of observed homogygote 1, count of observed heterozygotes,
         count of observed homogyhote 2.
  Output: Exact p-value for deviation (2-sided) from Hardy-Weinberg
          Proportions (HWP)
  Complexity: time and space O(min(hom1_count,hom2_count)+het_count)
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

  # Fill in relative probabilities fore greater than the expected hets
  for i,h in enumerate(xrange(hets,rare-1,2)):
    probs[h/2+1] = probs[h/2]*4*(hom_r-i)*(hom_c-i) / ((h+1)*(h+2))

  # Compute the pvalue by summing the probabilities <= to that of the
  # observed number of heterogygotes and normalize by the total
  p_obs = probs[het_count/2]
  pvalue = sum(p for p in probs if p <= p_obs)/sum(probs)

  return pvalue


def hwp_chisq_biallelic(hom1_count, het_count, hom2_count):
  '''Return the asymptotic Hardy-Weinberg Chi-squared value and p-value for the given genotypes'''

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

  return stats.distributions.chi2.sf(xx,1)


def hwp_biallelic(genos):
  return hwp_biallelic_counts(*count_genos(genos))


def hwp_biallelic_counts(hom1_count,het_count,hom2_count,exact_threshold=8000):

  # Only use the exact test when there are less than 8000 rare alleles
  # otherwise, use the asymptotic test
  if 2*min(hom1_count,hom2_count)+het_count < exact_threshold:
    p = hwp_exact_biallelic(hom1_count, het_count, hom2_count)
  else:
    p = hwp_chisq_biallelic(hom1_count, het_count, hom2_count)

  return p
