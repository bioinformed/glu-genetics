# -*- coding: utf-8 -*-
'''
File:          utils.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      generic utility functions

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import array
import warnings
import collections

from   functools import update_wrapper
from   itertools import izip, count, chain

from   fileutils import hyphen,namefile,autofile


def tally(seq):
  '''
  tally(sequence) -> { item:count,... }

  Returns a dictionary of values mapped to the number of times each
  item appears in the input sequence.
  '''
  d = collections.defaultdict(int)
  for item in seq:
    d[item] += 1
  return d


def ilen(seq):
  '''
  ilen(seq) -> length of iterable seq

  >>> ilen(xrange(10))
  10
  >>> ilen(a for a in range(10) for b in range(5))
  50
  >>> ilen(iter(range(100)))
  100
  '''
  return sum(1 for x in seq)


def pair_generator(items):
  '''Generator for distinct pairs of items'''
  if not isinstance(items, (list,tuple)):
    items = list(items)
  n = len(items)
  for i in xrange(n):
    for j in xrange(0,i):
      yield items[i],items[j]


def percent(a,b):
  if not b:
    return 0.
  return float(a)/b*100


def xenumerate(a, b=None):
  "enumerate([start,] iterable)"
  if b is None:
    start,iterable = 0,a
  else:
    start,iterable = a,b
  return izip(count(start), iterable)


def pick(sequence, indices):
  '''
  pick(sequence, indices) -> new sequence with elements chosen from the specified indices

  Supports strings, arrays, and sequences with contructors that take a single
  iterable argument.

  >>> pick([1,2,3], [0,2])
  [1, 3]
  >>> pick(range(10), range(0,10,2))
  [0, 2, 4, 6, 8]
  >>> pick('abcdefg', [1,3,5])
  'bdf'
  '''
  if isinstance(sequence,array.array):
    return type(sequence)(sequence.typecode, (sequence[i] for i in indices) )
  elif isinstance(sequence, str):
    return ''.join( (sequence[i] for i in indices) )
  else:
    return type(sequence)( (sequence[i] for i in indices) )


def peekfirst(sequence):
  '''
  peekfirst(sequence) -> first, sequence

  Returns the first element of the sequence and new iterator of the original sequence

  >>> peekfirst( [1,2,3] )
  (1, [1, 2, 3])
  >>> peekfirst( (1,2,3) )
  (1, (1, 2, 3))
  >>> f,s = peekfirst( iter([1,2,3]) )
  >>> f,list(s)
  (1, [1, 2, 3])
  '''
  try:
    return sequence[0],sequence
  except TypeError:
    pass

  sequence = iter(sequence)

  try:
    first = sequence.next()
  except StopIteration:
    raise IndexError('Cannot obtain first element of zero length sequence')

  sequence = chain([first],sequence)
  return first,sequence


def groups(sequence, keyfunc=None):
  '''
  Returns an iterator to the first element of each group of equal valued
  items (or keyfunc(item) if keyfunc is specified).

  This generator is merely a simplified version of itertools.groupby, such that:

  list(groups(seq,keyfunc)) == [ key for key,subiter in groupby(sequence,keyfunc) ]

  >>> list(groups( [1,1,2,2,3,3,1] ))
  [1, 2, 3, 1]
  '''
  if keyfunc is None:
    last = object()
    for item in sequence:
      if item != last:
        last = item
        yield item
  else:
    last = object()
    for item in sequence:
      key = keyfunc(item)
      if key != last:
        last = key
        yield item


def unique(sequence, keyfunc=None):
  '''
  Returns an iterator to the first occurance of each item in a sequence.  If
  keyfunc is specified, then the first occurance of each keyfunc(item) is
  returned.

  >>> list(unique( [1,1,2,2,3,3,1] ))
  [1, 2, 3]
  '''
  seen = set()
  if keyfunc is None:
    for item in sequence:
      if item not in seen:
        seen.add(item)
        yield item
  else:
    for item in sequence:
      key = keyfunc(item)
      if key not in seen:
        seen.add(key)
        yield item


def deprecated(func):
  '''
  This is a decorator which can be used to mark functions as deprecated. It
  will result in a warning being emitted when the function is used.

  >>> @deprecated
  ... def old_func(): return 1+1
  >>> old_func()    # doctest:+SKIP

  (Doctest skipped because warnings are not displayed properly)
  '''
  def deprecated_wrapper(*args, **kwargs):
    warnings.warn('Call to deprecated function %s.' % func.__name__,
                  category=DeprecationWarning, stacklevel=2)
    return func(*args, **kwargs)
  return update_wrapper(deprecated_wrapper, func)


def deprecated_by(msg):
  '''
  This is a decorator which can be used to mark functions as deprecated
  along with a message indicating by what. It will result in a warning being
  emitted when the function is used.

  >>> @deprecated_by('new_shiny_func')
  ... def old_func(): return 1+1
  >>> old_func()    # doctest:+SKIP

  (Doctest skipped because warnings are not displayed properly)
  '''
  def deprecated(func):
    def deprecated_wrapper(*args, **kwargs):
      warnings.warn('Call to function %s which has been deprecated by %s.' %
                    (func.__name__,msg), category=DeprecationWarning, stacklevel=2)
      return func(*args, **kwargs)
    return update_wrapper(deprecated_wrapper, func)

  return deprecated


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
