# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'generic utility functions'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


__all__ = ['as_set','is_str','tally','ilen','pair_generator','percent','xenumerate','pick',
           'peekfirst','groups','unique','izip_exact','deprecated','deprecated_by',
           'gcdisabled']

import gc
import array
import warnings
import collections

from   functools import update_wrapper
from   itertools import izip, count, chain


def as_set(items):
  '''
  Return items as a set if it not already a set or dict.

  @param  items: sequence, set, or mapping
  @rtype       : set or dict

  Examples:

  >>> s = set([1,2,3])
  >>> as_set(s) is s
  True

  >>> d = dict(a=1, b=2)
  >>> as_set(d) is d
  True

  >>> l = [1,2,3]
  >>> as_set(l) == set(l)
  True

  >>> t = (1,2,3)
  >>> as_set(t) == set(t)
  True

  >>> as_set(iter(t)) == set(t)
  True
  '''
  if isinstance(items, (dict,set)):
    return items
  return set(items)


def is_str(s):
  '''
  Return whether the input is a string

  @param  s: object
  @return  : indicate if the input is a string
  @rtype   : bool

  >>> is_str(True)
  False

  >>> is_str(None)
  False

  >>> is_str('abc')
  True

  >>> is_str(u'abc')
  True
  '''
  return isinstance(s, basestring)


def tally(seq):
  '''
  tally(sequence) -> { item:count,... }

  Returns a dictionary of values mapped to the number of times each
  item appears in the input sequence.

  @param  seq: input sequence
  @type   seq: sequence
  @return    : dictionary indicating the count for each key item
  @rtype     : dict

  >>> tally(['A','B','A','C','A','C','D'])
  defaultdict(<type 'int'>, {'A': 3, 'C': 2, 'B': 1, 'D': 1})
  >>> tally(('A','B','A','C','A','C','D'))
  defaultdict(<type 'int'>, {'A': 3, 'C': 2, 'B': 1, 'D': 1})
  '''
  d = collections.defaultdict(int)
  for item in seq:
    d[item] += 1
  return d


def ilen(seq):
  '''
  ilen(seq) -> length of iterable seq

  @param  seq: input sequence
  @type   seq: sequence
  @return    : length of the seq
  @rtype     : int

  >>> ilen(xrange(10))
  10
  >>> ilen(a for a in range(10) for b in range(5))
  50
  >>> ilen(iter(range(100)))
  100
  '''
  return sum(1 for x in seq)


def pair_generator(items):
  '''
  Generator for distinct pairs of items

  @param items: sequence of items
  @type  items: sequence
  @return     : distinct pairs of items
  @rtype      : generator of tuples

  >>> pairs = pair_generator(['A','B','C'])
  >>> for item1,item2 in pairs:
  ...   print item1,item2
  B A
  C A
  C B
  '''
  if not isinstance(items, (list,tuple)):
    items = list(items)
  n = len(items)
  for i in xrange(n):
    for j in xrange(0,i):
      yield items[i],items[j]


def percent(a,b):
  '''
  Return the percentage of two numbers
  '''
  if not b:
    return 0.
  return float(a)/b*100


def xenumerate(a, b=None):
  '''
  Return an enumerate object with a start position
  enumerate([start,] iterable)
  '''
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

  @param  sequence: input sequence
  @type   sequence: sequence
  @param   indices: sequence of indices
  @type    indices: sequence of int
  @return         : sequence with elements at the specified indices
  @rtype          : sequence

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

  @param  sequence: input sequence
  @type   sequence: sequence
  @return         : sequence with first element and the original sequence
  @rtype          : sequence

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

  @param  sequence: input sequence
  @type   sequence: sequence
  @param   keyfunc: converting function for grouping
  @type    keyfunc: function
  @return         : new sequence
  @rtype          : sequence

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
  Generator to produce the unique first occurance of each item in a
  sequence.  If keyfunc is specified, then the first occurance of each
  keyfunc(item) is returned, otherwise the entire item is taken as a key.
  Ordering is stable, since result elements will always appear in the order
  they first first appear in the input sequence.  Keys must be hashable.

  @param  sequence: input sequence
  @type   sequence: sequence
  @param   keyfunc: key function function
  @type    keyfunc: callable
  @return         : new sequence
  @rtype          : generator

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


# Helpers for izip_exact

class LengthMismatch(Exception):
  pass

def _izip_exact_thow():
  raise LengthMismatch
  yield None # unreachable

def _zip_exact_check(rest):
  for i in rest:
    try:
      i.next()
    except LengthMismatch:
      pass
    else:
      raise LengthMismatch
  return
  yield None # unreachable


def izip_exact(*iterables):
  '''
  izip_exact(iter1 [,iter2 [...]]) --> iterator object

  Return an iterator whose .next() method returns a tuple where the i-th
  element comes from the i-th iterable argument.  The .next() method
  continues until the shortest iterable in the argument sequence is
  exhausted.  If all sequences are exhausted a StopIteration exception is
  raised, otherwise a LengthMismatch exception is raised.

  Works like itertools.izip(), but throws a LengthMismatch exception if any
  iterable's length differs.

  Zero length iterable lists return iter(iterables).  Unit length iterables
  return iter(iterables[0]).  Otherwise an izip object is returned.

  If the len() of all iterables is available and match, then this function
  returns izip(*iterables).  In the length of any iterable cannot be
  determined, then sentinel objects are appended to each iterable and an
  izip object of these augmented iterables is returned.  If the length a
  proper subset of iterables can be determined, the second strategy is
  employed, even if known lengths do not match.  This is in anticipation of
  iterable side-effects that may ultimately balance the lengths.

  Inspired largely by Peter Otten's zip_exc, modified to check lengths.
  (http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/497006)

  >>> list(izip_exact())
  []
  >>> list(izip_exact([]))
  []
  >>> list(izip_exact((), (), ()))
  []

  >>> list(izip_exact("abc", range(3)))
  [('a', 0), ('b', 1), ('c', 2)]

  >>> list(izip_exact("", range(3)))
  Traceback (most recent call last):
       ...
  LengthMismatch

  >>> list(izip_exact(range(3), ()))
  Traceback (most recent call last):
       ...
  LengthMismatch

  >>> list(izip_exact(range(3), range(2), range(4)))
  Traceback (most recent call last):
       ...
  LengthMismatch

  >>> items = izip_exact(iter(range(3)), range(2), range(4))
  >>> items.next()
  (0, 0, 0)
  >>> items.next()
  (1, 1, 1)
  >>> items.next()
  Traceback (most recent call last):
       ...
  LengthMismatch
  '''
  if not iterables:
    return iter(iterables)
  elif len(iterables) == 1:
    return iter(iterables[0])

  first = iterables[0]
  rest  = iterables[1:]

  try:
    n = len(first)
    if all( len(i)==n for i in rest ):
      return izip(*iterables)
  except (TypeError,AttributeError):
    pass

  # Must use sentinel objects to enforce length equality
  rest  = [chain(i, _izip_exact_thow()) for i in rest]
  first = chain(first, _zip_exact_check(rest))
  return izip(*[first] + rest)


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


class gcdisabled(object):
  '''
  Conext manager to temporarily disable Python's cyclic garbage collector.
  The primary use is to avoid thrashing while allocating large numbers of
  non-cyclic objects due to an overly aggressive garbage collector behavior.

  Will disable GC if it is enabled upon entry and renable upon exit:

  >>> gc.isenabled()
  True
  >>> with gcdisabled():
  ...   print gc.isenabled()
  False
  >>> print gc.isenabled()
  True

  Will not reenable if GC was disabled upon entry:

  >>> gc.disable()
  >>> gc.isenabled()
  False
  >>> with gcdisabled():
  ...   gc.isenabled()
  False
  >>> gc.isenabled()
  False
  '''
  def __init__(self):
    self.isenabled = gc.isenabled()

  def __enter__(self):
    gc.disable()

  def __exit__(self, type, value, traceback):
    if self.isenabled:
      gc.enable()


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
