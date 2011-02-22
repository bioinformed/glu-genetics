# -*- coding: utf-8 -*-

__abstract__  = 'efficiently merge sorted iterables'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   operator  import itemgetter
from   itertools import tee,izip,imap,count


def imerge(its, key=None):
  '''
  Generator to efficiently merge sorted iterables with optional key
  function.  Input sequence must be in native order if no key function is
  not specified, or in key order if a key function is specified.

  Equivalent, but much more efficient than:

    iter(sorted(*chain(*its), key=key))

  so long as the following invariant holds:

    for it in its:
      it = list(it)
      assert it == sorted(it,key=key)

  @param  its: list of sorted input sequences
  @type   its: sequence
  @param  key: optional key comparison function
  @type   key: binary function
  @return:     sequence of sorted values
  @rtype:      iterator

  >>> list(imerge([range(5),range(5,11)]))
  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  >>> list(imerge([reversed(range(5)),reversed(range(5,11))],key=lambda x: -x))
  [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
  '''
  if key is None:
    return merge(*its)

  counter = count()
  its_d   = []

  for it in its:
    it1,it2 = tee(it)
    its_d.append( izip(imap(key,it1), counter, it2) )

  merged = merge(*its_d)

  return imap(itemgetter(2), merged)


try:
  from heapq import merge

except ImportError:
  from heapq import heapify, heappop, heapreplace

  def merge(*iterables):
    '''Merge multiple sorted inputs into a single sorted output.

    Similar to sorted(itertools.chain(*iterables)) but returns a generator,
    does not pull the data into memory all at once, and assumes that each of
    the input streams is already sorted (smallest to largest).

    Based on heapq.merge from Python 2.7.

    Equivalent, but much more efficient than:

      iter(sorted(*chain(*its)))

    so long as the following invariant holds:

      for it in its:
        it = list(it)
        assert it == sorted(it)

    @param its: list of sorted input sequences
    @type  its: sequence
    @return:    sequence of sorted values
    @rtype:     iterator

    >>> list(merge(range(5),range(5,11)))
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    >>> list(merge([1,3,5,7], [0,2,4,8], [5,10,15,20], [], [25]))
    [0, 1, 2, 3, 4, 5, 5, 7, 8, 10, 15, 20, 25]

    '''
    _heappop, _heapreplace, _StopIteration = heappop, heapreplace, StopIteration

    h = []
    h_append = h.append

    for itnum,it in enumerate(map(iter, iterables)):
      try:
        next = it.next
        h_append([next(), itnum, next])
      except _StopIteration:
        pass

    heapify(h)

    while 1:
      try:
        while 1:
          v, itnum, next = s = h[0]   # raises IndexError when h is empty
          yield v
          s[0] = next()               # raises StopIteration when exhausted
          _heapreplace(h, s)          # restore heap condition
      except _StopIteration:
        _heappop(h)                   # remove empty iterator
      except IndexError:
        return


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
