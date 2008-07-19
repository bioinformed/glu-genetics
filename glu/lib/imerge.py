# -*- coding: utf-8 -*-

__abstract__  = 'efficiently merge sorted iterables'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import heapq
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
    return _imerge(its)

  counter = count()
  its_d = []
  for it in its:
    it1,it2 = tee(it)
    its_d.append( izip(imap(key,it1), counter, it2) )

  merged = _imerge(its_d)

  return imap(itemgetter(2), merged)


def _imerge(its):
  '''
  Generator to efficiently merge sorted iterables using native Python sort
  order.

  Based on a non-micro optimized version of Raymond Hettinger's recipe at
  http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/491285

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

  >>> list(_imerge([range(5),range(5,11)]))
  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  '''
  pqueue = []
  for i in imap(iter, its):
    try:
      pqueue.append((i.next(), i.next))
    except StopIteration:
      pass

  heapq.heapify(pqueue)

  while pqueue:
    val, it = pqueue[0]
    yield val
    try:
      heapq.heapreplace(pqueue, (it(), it))
    except StopIteration:
      heapq.heappop(pqueue)


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
