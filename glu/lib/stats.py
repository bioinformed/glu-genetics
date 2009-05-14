# -*- coding: utf-8 -*-

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
  @param         k: the perventile value in the range [0..1]
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


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
