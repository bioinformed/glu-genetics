# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'generic utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


__all__ = ['as_set','is_str','tally','ilen','pair_generator','percent','xenumerate','pick',
           'peekfirst','groups','unique','izip_exact','LengthMismatch','deprecated','deprecated_by',
           'gcdisabled','chunk']

import gc

from   collections      import deque
from   itertools        import izip, count, chain, islice, repeat, imap

import glu.lib.pycompat
from   glu.lib.pycompat import *

__all__ += glu.lib.pycompat.__all__


def as_set(items):
  '''
  Return items as a set if it not already a set or dict.  An exception is
  made if items is a bare string, where a set with a single item is returned
  (rather than the set of characters in items).

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

  >>> as_set('ab') == set(['ab'])
  True
  '''
  if is_str(items):
    items = [items]
  elif isinstance(items, (dict,set)):
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


def constant_function(value):
  '''
  constant_function(value) -> lambda: value

  Return a trivial function that takes no arguments and returns only value.
  This is semantically equivalent to "lambda : value", except roughly twice
  as fast (at the expense of being obscure (not that lambda isn't obscure
  these days)).

  Based on: http://code.activestate.com/recipes/502215/

  Why bother?  The primary use-case envisioned is to provide a more
  efficient default values to collections.defaultdict.

  >>> c = constant_function(1)
  >>> c()
  1
  >>> c()
  1

  >>> c = constant_function('foo')
  >>> c()
  'foo'
  >>> c()
  'foo'

  >>> from collections import defaultdict
  >>> d = defaultdict(c)
  >>> d[99]
  'foo'
  '''
  return repeat(value).next


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
  Counter({'A': 3, 'C': 2, 'B': 1, 'D': 1})
  >>> tally(('A','B','A','C','A','C','D'))
  Counter({'A': 3, 'C': 2, 'B': 1, 'D': 1})
  '''
  return Counter(seq)


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


def consume(seq, n=None):
  '''
  consume(seq, n=None)

  Advance the iterator n-steps ahead. If n is none, consume entirely.

  Based on itertools recipes found in the Python documentation.  See:
  http://docs.python.org/library/itertools.html

  @param seq: input
  @type  seq: sequence
  @param   n: count of steps or None
  @type    n: int or NoneType

  >>> s = iter('abcdefg')
  >>> consume(s,3)
  >>> print ''.join(s)
  defg

  >>> s = iter('abcdefg')
  >>> consume(s)
  >>> print ''.join(s)
  <BLANKLINE>
  '''
  from collections import deque

  # The technique uses objects that consume iterators at C speed.
  if n is None:
    # feed the entire iterator into a zero-length deque
    deque(seq, maxlen=0)
  else:
    # advance to the empty slice starting at position n
    next(islice(seq, n, n), None)


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

  Supports strings, arrays, and sequences with constructors that take a single
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
  import array

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


def groups(sequence, key=None):
  '''
  Returns an iterator to the first element of each group of equal valued
  items (or key(item) if key is specified).

  This generator is merely a simplified version of itertools.groupby, such that:

  list(groups(seq,key)) == [ key for key,subiter in groupby(sequence,key) ]

  @param sequence: input sequence
  @type  sequence: sequence
  @param      key: converting function for grouping
  @type       key: function
  @return        : new sequence
  @rtype         : sequence

  >>> list(groups( [1,1,2,2,3,3,1] ))
  [1, 2, 3, 1]
  '''
  if key is None:
    last = object()
    for item in sequence:
      if item != last:
        last = item
        yield item
  else:
    last = object()
    for item in sequence:
      k = key(item)
      if k != last:
        last = k
        yield item


def unique(sequence, key=None):
  '''
  Generator to produce the unique first occurrence of each item in a
  sequence.  If key is specified, then the first occurrence of each
  key(item) is returned, otherwise the entire item is taken as a key.
  Ordering is stable, since result elements will always appear in the order
  they first first appear in the input sequence.  Keys must be hashable.

  @param  sequence: input sequence
  @type   sequence: sequence
  @param       key: key function function
  @type        key: callable
  @return         : new sequence
  @rtype          : generator

  >>> list(unique( [1,1,2,2,3,3,1] ))
  [1, 2, 3]
  '''
  seen = set()
  if key is None:
    for item in sequence:
      if item not in seen:
        seen.add(item)
        yield item
  else:
    for item in sequence:
      k = key(item)
      if k not in seen:
        seen.add(k)
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


class OrderError(ValueError): pass


def sort_almost_sorted(sequence,key=None,windowsize=1000,stable=True):
  '''
  sort_almost_sorted(sequence,key=None,windowsize=1000)

  Sorts an almost sorted sequence of items provided that all misordered
  items are within windowsize elements of the correct location in the final
  sequence.  Returns a generator that yields the correctly sorted sequence
  or raises OrderError if the correct sequence cannot be constructed with
  the given window size.

  If a key function is provided or the stable argument is True, then the
  resulting sort order is stable, with items returned in proper sort order
  and equal keys returned in the same relative order.  Otherwise, the
  resulting ordering is not guaranteed to be stable.

  Test window sizes:

    >>> list(sort_almost_sorted([1,2,4,5,3],windowsize=3))
    [1, 2, 3, 4, 5]

    >>> list(sort_almost_sorted([1,2,5,6,4],windowsize=1))
    Traceback (most recent call last):
         ...
    OrderError: Misordered keys beyond window size

    >>> list(sort_almost_sorted([2,3,4,5,6,7,9,10,1],windowsize=10))
    [1, 2, 3, 4, 5, 6, 7, 9, 10]

    >>> list(sort_almost_sorted([2,3,4,5,6,7,9,10,1],windowsize=8))
    Traceback (most recent call last):
    ...
    OrderError: Misordered keys beyond window size

  Test stability:

    # key=int(x), so all keys will compare equal and values are not in the
    # natural ordering.
    >>> list(sort_almost_sorted([1.7,1.5,1.6,1.4], key=lambda x: int(x)))
    [1.7, 1.5, 1.6, 1.4]

  Test key function:

    >>> list(sort_almost_sorted([1,2,3,4], key=lambda x: -x))
    [4, 3, 2, 1]
  '''
  from operator import itemgetter

  # Key sorts are always stable, since we already pay the price for a
  # decorated sequence.
  if key is not None:
    decorated = ( (key(item),i,item) for i,item in enumerate(sequence) )
    ordered   = _sort_almost_sorted(decorated,windowsize)
    return imap(itemgetter(2),ordered)

  # Otherwise, use a similar method as above to ensure stability
  elif stable:
    decorated = ( (item,i) for i,item in enumerate(sequence) )
    ordered   = _sort_almost_sorted(decorated,windowsize)
    return imap(itemgetter(0),ordered)

  # Unstable, undecorated sort
  else:
    return _sort_almost_sorted(sequence,windowsize)


def _sort_almost_sorted(sequence,windowsize):
  '''
  Internal function.  See sort_almost_sorted
  '''
  from heapq import heapify, heappushpop, heappop

  # STAGE 1: Fill initial window and heapify
  it   = iter(sequence)
  heap = list(islice(it,windowsize))

  # Establish invariant len(heap)>0
  if not heap:
    return

  heapify(heap)

  # STAGE 2: Slide window until end of sequence
  last = heap[0]

  # Loop invariants:
  #   len(heap)==c, where c is a constant 0 < c <= windowsize
  #   all(e>=last for e in heap)
  for item in it:
    item = heappushpop(heap, item)

    if item<last:
      raise OrderError('Misordered keys beyond window size')

    last = item

    yield item

  # Release reference to last item, since it is no longer needed
  del last

  # STAGE 3: Drain window, no need to check sort order of remaining elements
  #          since remaining elements must be within windowsize distance
  while heap:
    yield heappop(heap)


def check_sorted_iter(sequence,key=None):
  '''
  check_sorted_iter(sequence,key=None) -> sequence

  Returns a generator that yields all elements of sequence, provided that
  the elements are sorted in non-descending order.  Otherwise an OrderError
  exception is raised.

  >>> list(check_sorted_iter([1,2,3]))
  [1, 2, 3]

  >>> list(check_sorted_iter([3,2,1]))
  Traceback (most recent call last):
       ...
  OrderError: Invalid sort order

  >>> list(check_sorted_iter([1.7,1.5,1.6,1.4], key=lambda x: int(x)))
  [1.7, 1.5, 1.6, 1.4]
  '''
  it = iter(sequence)

  try:
    last = next(it)
  except StopIteration:
    return

  if key is None:
    yield last

    for item in it:
      if item<last:
        raise OrderError('Invalid sort order')

      last = item

      yield item
  else:
    yield last
    last = key(last)

    for item in it:
      current = key(item)

      if current<last:
        raise OrderError('Invalid sort order')

      last = current

      yield item


def parallel_coordinate_iterator(sequences, key):
  '''
  Iterate over several ordered sequences and return a list of values at each
  distinct key value.

  >>> from operator import itemgetter
  >>> seq1 = [(1,'A'),(2,'B'),        (4,'D')]
  >>> seq2 = [        (2,'b'),                (5,'E')]
  >>> seq3 = [(1,'W'),(2,'X'),(3,'Y'),        (5,'Z')]
  >>> seqs = (seq1,seq2,seq3)
  >>> p    = parallel_coordinate_iterator(seqs,itemgetter(0))
  >>> for row in p:
  ...   print row
  [(1, 'A'), None, (1, 'W')]
  [(2, 'B'), (2, 'b'), (2, 'X')]
  [None, None, (3, 'Y')]
  [(4, 'D'), None, None]
  [None, (5, 'E'), (5, 'Z')]
  '''
  from heapq     import merge
  from itertools import groupby
  from operator  import itemgetter

  n = len(sequences)

  # Must be a function or else i will be unbound across sequences and will
  # vary with outer iteration
  def _decorate(seq,index):
    return ((key(x),index,x) for x in seq)

  iters = [ _decorate(s,i) for i,s in enumerate(sequences) ]

  for _,items in groupby(merge(*iters),itemgetter(0)):
    row = [None]*n

    for k,i,value in items:
      if row[i] is not None:
        row = None
        break
      row[i] = value

    if row:
      yield row


def deprecated(func):
  '''
  This is a decorator which can be used to mark functions as deprecated. It
  will result in a warning being emitted when the function is used.

  >>> @deprecated
  ... def old_func(): return 1+1
  >>> old_func()    # doctest:+SKIP

  (Doctest skipped because warnings are not displayed properly)
  '''
  from functools import update_wrapper

  def deprecated_wrapper(*args, **kwargs):
    import warnings
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
  from functools import update_wrapper

  def deprecated(func):
    def deprecated_wrapper(*args, **kwargs):
      import warnings
      warnings.warn('Call to function %s which has been deprecated by %s.' %
                    (func.__name__,msg), category=DeprecationWarning, stacklevel=2)
      return func(*args, **kwargs)
    return update_wrapper(deprecated_wrapper, func)

  return deprecated


class gcdisabled(object):
  '''
  Context manager to temporarily disable Python's cyclic garbage collector.
  The primary use is to avoid thrashing while allocating large numbers of
  non-cyclic objects due to an overly aggressive garbage collector behavior.

  Will disable GC if it is enabled upon entry and re-enable upon exit:

  >>> gc.isenabled()
  True
  >>> with gcdisabled():
  ...   print gc.isenabled()
  False
  >>> print gc.isenabled()
  True

  Will not re-enable if GC was disabled upon entry:

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


_chunk_nothing=object()

def chunk(iterable, size, pad=_chunk_nothing):
  '''
  Generator that breaks an iterable into tuples of chunks (tuples) of a
  given size.  The optional keyword parameter 'pad' if the last block is
  padded with a given value or returned as-is if smaller than a full
  chunk.

  >>> list(chunk(xrange(7), 3))
  [(0, 1, 2), (3, 4, 5), (6,)]
  >>> list(chunk(xrange(7), 3, pad=None))
  [(0, 1, 2), (3, 4, 5), (6, None, None)]
  >>> list(chunk(reversed(xrange(7)), 3, pad=-1))
  [(6, 5, 4), (3, 2, 1), (0, -1, -1)]
  '''
  iterable = iter(iterable)
  while True:
    block = tuple(islice(iterable,size))

    if not block:
      break

    if pad is not _chunk_nothing and len(block) < size:
      block = tuple(chain(block, repeat(pad, size-len(block))))

    yield block


def chunk_iter(iterable, size, pad=_chunk_nothing):
  for c in chunk(iterable, size, pad):
    for item in c:
      yield item


def enum(typename, field_names):
  '''
  Create a new enumeration type

  Based on Recipe 577024: "Yet another 'enum' for Python" from
  http://code.activestate.com/recipes/577024-yet-another-enum-for-python/ by
  Gabriel Genellina.

  >>> STATE = enum('STATE', 'GET_QUIZ, GET_VERSE, TEACH')
  >>> STATE.GET_QUIZ
  0
  >>> STATE.GET_VERSE
  1
  >>> STATE.TEACH
  2
  >>> STATE.GET_VERSE = 8
  Traceback (most recent call last):
    ...
  AttributeError: 'STATE' object attribute 'GET_VERSE' is read-only
  >>> del STATE.GET_VERSE
  Traceback (most recent call last):
    ...
  AttributeError: 'STATE' object attribute 'GET_VERSE' is read-only
  '''
  if isinstance(field_names, str):
      field_names = field_names.replace(',', ' ').split()
  d = dict((reversed(nv) for nv in enumerate(field_names)), __slots__ = ())
  return type(typename, (object,), d)()


class iter_queue(deque):

  '''Iterator object that is also a queue that allows one to look-ahead in
  the stream using the peek() method.  Calling peek() when the iterator is
  exhausted results in a ValueError exception.

  Values returned by peek are accessible from the other using a
  deque-interface, from which the iterator inherits.  Each call to peek()
  looks ahead one value and that value is available as the right-most
  element in the deque.  There is no look-behind, as all buffered values
  returned from the iterator are removed after they are "seen" via the
  iterator protocol.  Thus the len(it) returns the number of items in the
  look-ahead and next(it) reduces the size of the look-ahead by one, if any
  are present.

  >>> iq = iter_queue(xrange(10))
  >>> iq.peek()
  0
  >>> iq.peek()
  1
  >>> iq[1]
  1
  >>> next(iq)
  0
  >>> iq[0]
  1
  >>> iq.peek(0)
  1
  >>> iq.peek()
  2
  >>> next(iq)
  1
  >>> iq.peek(4)
  6
  >>> list(iq)
  [2, 3, 4, 5, 6, 7, 8, 9]
  '''
  def __init__(self, iterable):
    self.iterable = iter(iterable)

  def __iter__(self):
    return self

  def next(self):
    try:
      return self.popleft()
    except IndexError:
      return next(self.iterable)

  def peek(self, i=None):
    '''Return next or i'th value but do not advance this iterator'''
    try:
      if i is not None:
        # Index is already buffered
        if i<len(self):
          value = self[i]

        # We must peek further
        else:
          for _ in xrange(i+1-len(self)):
            value = next(self.iterable)
            self.append(value)
      else:
        value = next(self.iterable)
        self.append(value)

    except StopIteration:
      raise ValueError('No more values')

    return value


def generator_thread(source, maxsize=1000, chunksize=1000):
  '''
  Run a generator in a seperate thread


  >>> list(generator_thread(range(10)))
  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  '''
  import threading
  import Queue

  def sendto_queue(source,q):
    put = q.put
    if chunksize<=1:
      for item in source:
        put(item)
    else:
      source = iter(source)
      items  = True
      while items:
        items = list(islice(source, chunksize))
        put(items)

    put(StopIteration)

  def genfrom_queue(q):
    _StopIteration = StopIteration
    get = q.get

    if chunksize<=1:
      while True:
        item = get()
        if item is _StopIteration:
          break
        yield item
    else:
      while True:
        items = get()
        if items is _StopIteration:
          break
        for item in items:
          yield item

  q = Queue.Queue(maxsize)

  thr = threading.Thread(target=sendto_queue, args=(source,q))
  thr.setDaemon(True)
  thr.start()

  return genfrom_queue(q)


def lazy_property(fn):
  '''
  Return a lazily computed property.

  This property will compute a read-only result once, but only when
  requested.  All subsequent requests will used the stored result of the
  first computation.
  '''
  attr_name = '_lazy_' + fn.__name__
  @property
  def _lazy_property(self):
      if not hasattr(self, attr_name):
          setattr(self, attr_name, fn(self))
      return getattr(self, attr_name)
  return _lazy_property


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
