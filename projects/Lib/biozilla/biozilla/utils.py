# -*- coding: utf-8 -*-
'''
File:          utils.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      Biozilla generic utility functions

Requires:      Python 2.5

Revision:      $Id: utils.py 510 2007-02-27 17:23:28Z jacobske $
'''

__version__ = '0.99'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import array

from   itertools import izip, count, chain

from   fileutils import hyphen,namefile,autofile


try:
  any=any

except NameError:
  def any(it):
    '''Return True if bool(x) is True for any x in the iterable.'''
    for i in it:
      if i:
        return True
    return False


try:
  all=all

except NameError:
  def all(it):
    '''Return True if bool(x) is True for all x in the iterable.'''
    for i in it:
      if not i:
        return False
    return True


try:
  tally=tally

except NameError:
  import collections

  # Python 2.5 has a convenient defaultdict data structure
  if hasattr(collections,'defaultdict'):
    def tally(seq):
      '''
         tally(sequence) -> { item:count,... }

         Returns a dictionary of values mapped to the number of times each
         item appears in the input sequence.
      '''
      d = collections.defaultdict(int)
      for item in seq:
        d[item] += 1
      return dict(d)

  # Python 2.4 does not, so fall back to the old and ugly method
  else:
    def tally(seq):
      '''
         tally(sequence) -> { item:count,... }

         Returns a dictionary of values mapped to the number of times each
         item appears in the input sequence.
      '''
      d = {}
      for item in seq:
        d[item] = d.setdefault(item,0) + 1
      return d


def ilen(seq):
  '''ilen(seq) -> length of iterable seq'''
  count = 0
  for x in seq:
    count += 1
  return count


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


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
