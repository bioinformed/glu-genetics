# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'miscellaneous statistical functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from math import modf


def mean(seq):
  '''
  Find the mean value for the input sequence

  @param seq: a sequence of numbers
  @type  seq: sequence such as list
  @return:    mean for the sequence and ValueError for 0 length sequences
  @rtype:     float

  >>> mean([5])
  5.0
  >>> mean(range(10))
  4.5
  >>> mean([])
  Traceback (most recent call last):
       ...
  ValueError: Input sequence may not be empty
  '''
  seq = iter(seq)
  total = 0
  i = None
  for i,v in enumerate(seq):
    total += v

  if i is None:
    raise ValueError('Input sequence may not be empty')

  return float(total)/(i+1)


def median(seq, presorted=False):
  '''
  Return the median values of the data.  If the data are already sorted,
  then the presorted parameter should be set to avoid resorting the data.
  Otherwise, a copy will be made and then sorted.

  @param      data: sequence of data elements that can be meaningfully ordered
  @type       data: if presorted, sequence.  If not presorted, iterable
  @param presorted: Indicator if the input data sequence is supplied in sorted order
                    (default=False)
  @type  presorted: bool
  @return         : median value
  '''
  return quantile(seq, 0.5, presorted)


def quantile(data, k, presorted=False):
  '''
  Return the k-th weighted quantile of values in data.  If the data are
  already sorted, then the presorted parameter should be set to avoid
  resorting the data.  Otherwise, a copy will be made and then sorted.

  @param      data: sequence of data elements that can be meaningfully ordered
  @type       data: if presorted, sequence.  If not presorted, iterable
  @param         k: the percentile value in the range [0..1]
  @type          k: float
  @param presorted: Indicator if the input data sequence is supplied in sorted order
                    (default=False)
  @type  presorted: bool
  @return         : k-th weighted quantile of data

  Quantiles of data=0..16 from 0..1 every 0.5
  >>> observed = [ quantile(range(17), i/20.) for i in range(21) ]
  >>> expected = [0, 0.8, 1.6, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 8,
  ...             8.8, 9.6, 10.4, 11.2, 12, 12.8, 13.6, 14.4, 15.2, 16]
  >>> # Verify that the sum of the errors is less than 10^-20
  >>> assert sum( abs(o-e) for o,e in zip(observed,expected) ) < 10e-20
  '''
  if not (0 <= k <= 1):
    raise ValueError('Quantile value %f out of range [0..1]' % k)

  if not presorted:
    data = sorted(data)

  if not data:
    raise ValueError('Input sequence may not be empty')

  n = len(data) - 1
  f,w = modf(n*k)
  w = int(w)

  if f < 1e-10:
    q = data[w]
  else:
    q = (1-f)*data[w] + f*data[w+1]

  return q


def wilson_score_interval(x,n,alpha=0.05):
  '''
  wilson_score_interval(x,n,alpha=0.05) -> frequency, lower limit, upper limit

  Compute the binomial proportion confidence interval for a proportion x/n
  for power alpha or, equivalently, a confidence interval of 1-alpha.

  The Wilson interval is an improvement (the actual coverage probability is
  closer to the nominal value) over the normal approximation interval and
  was first developed by Edwin Bidwell Wilson (1927).

  For more information, see:

    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

    Wilson, E. B. (1927). "Probable inference, the law of succession, and
    statistical inference".  Journal of the American Statistical Association
    22: 209–212.  JSTOR 2276774.

    Agresti, Alan; Coull, Brent A. (1998). "Approximate is better than
    'exact' for interval estimation of binomial proportions".  The American
    Statistician 52: 119–126.  doi:10.2307/2685469.  JSTOR 2685469.
    MR1628435.

  >>> wilson_score_interval(0,40)
  (0.0, 0.0, 0.08762160119728668)
  >>> wilson_score_interval(1,40)
  (0.025, 0.0012823323596887644, 0.12881368963474096)
  >>> wilson_score_interval(39,40)
  (0.975, 0.87118631036525906, 0.9987176676403112)
  >>> wilson_score_interval(40,40)
  (1.0, 0.91237839880271343, 1.0)
  '''
  import scipy.stats
  from   math import log

  p = x/n
  z = -scipy.stats.distributions.norm.ppf(alpha/2)

  if x==1:
    p_lower =   - log(1-alpha)/n
  else:
    p_lower = (p + z*z/n/2 - z*( (p*(1-p) + z*z/4/n)/n )**0.5)/(1+z*z/n)

  if x==n-1:
    p_upper = 1 + log(1-alpha)/n
  else:
    p_upper = (p + z*z/n/2 + z*( (p*(1-p) + z*z/4/n)/n )**0.5)/(1+z*z/n)


  p_lower = max(0.,p_lower)
  p_upper = min(1.,p_upper)

  return p,p_lower,p_upper


def jeffreys_interval(x,n,alpha=0.05):
  '''
  jeffreys_interval(x,n,alpha=0.05) -> frequency, lower limit, upper limit

  Compute the binomial proportion confidence interval for a proportion x/n
  for power alpha or, equivalently, a confidence interval of 1-alpha.

  The 'Jeffreys interval' is the Bayesian credible interval obtained when
  using the non-informative Jeffreys prior for the binomial proportion p.
  The Jeffreys prior for this problem is a Beta distribution with parameters
  (1/2, 1/2).  After observing x successes in n trials, the posterior
  distribution for p is a Beta distribution with parameters (x + 1/2, n – x
  + 1/2).

  When x ≠ 0 and x ≠ n, the Jeffreys interval is taken to be the 100(1 – α)%
  equal-tailed posterior probability interval, i.e.  the α/2 and 1 – α/2
  quantiles of a Beta distribution with parameters (x + 1/2, n – x + 1/2).
  In order to avoid the coverage probability tending to zero when p → 0 or
  1, when x = 0 the upper limit is calculated as before but the lower limit
  is set to 0, and when x = n the lower limit is calculated as before but
  the upper limit is set to 1.

  For more information, see:

    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval

    Brown, Lawrence D.; Cai, T. Tony; DasGupta, Anirban (2001). "Interval
    Estimation for a Binomial Proportion".  Statistical Science 16 (2):
    101–133.

  >>> jeffreys_interval(0,40)
  (0.0, 0.0, 0.060497975214079638)
  >>> jeffreys_interval(1,40)
  (0.025, 0.0012823323596887644, 0.1109493804784687)
  >>> jeffreys_interval(39,40)
  (0.975, 0.8890506195215313, 0.9987176676403112)
  >>> jeffreys_interval(40,40)
  (1.0, 0.93950202478592038, 1.0)
  '''
  import scipy.stats
  from   math import log

  p    = x/n
  beta = scipy.stats.distributions.beta(x+0.5,n-x+0.5)

  if x==0:
    p_lower = 0.0
  elif x==1:
    p_lower =   - log(1-alpha)/n
  else:
    p_lower = beta.ppf(  alpha/2)

  if x==n:
    p_upper = 1.0
  elif x==n-1:
    p_upper = 1 + log(1-alpha)/n
  else:
    p_upper = beta.ppf(1-alpha/2)

  p_lower = max(0.,p_lower)
  p_upper = min(1.,p_upper)

  return p,p_lower,p_upper


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
