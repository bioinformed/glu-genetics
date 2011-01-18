# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Various edit distance algorithms'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import numpy as np

from libc.stdio    cimport sprintf
from libc.stdlib   cimport malloc, realloc, free
from cpython       cimport PyString_Size, PyErr_NoMemory


cdef inline Py_ssize_t min2(Py_ssize_t a, Py_ssize_t b):
  if a < b:
    return a
  else:
    return b


cdef inline Py_ssize_t min3(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c):
  return min2(min2(a,b),c)


cdef inline Py_ssize_t min4(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c, Py_ssize_t d):
  return min2(min2(a,b),min2(c,d))


def hamming_distance(s1, s2):
  '''
  Calculate the Hamming distance between two sequences.

  This distance is the number of substitutions needed to transform the first
  sequence into the second.  Both sequences are required to be of equal
  length.

  See: http://en.wikipedia.org/wiki/Hamming_distance

       Hamming, Richard W. (1950), "Error detecting and error correcting
       codes", Bell System Technical Journal 26 (2): 147–160,
       http://www.ece.rutgers.edu/~bushnell/dsdwebsite/hamming.pdf

  This implementation requires O(N) time and O(1) space, where N is the
  length of the input sequences.

  >>> hamming_distance('abc', 'abc')
  0
  >>> hamming_distance('abb', 'abc')
  1
  >>> hamming_distance('abc', 'def')
  3
  >>> hamming_distance('a', 'ab')
  Traceback (most recent call last):
  ...
  LengthMismatch
  '''
  cdef Py_ssize_t d = 0
  cdef Py_ssize_t l1 = PyString_Size(s1)
  cdef Py_ssize_t l2 = PyString_Size(s2)

  if l1!=l2:
    raise ValueError('Length mismatch')

  cdef char *ss1 = s1
  cdef char *ss2 = s2

  for i in range(l1):
    if ss1[i]!=ss2[i]:
      d += 1

  return d


cdef reduce_match(s1, s2):
  '''
  Reduce matching problem size by removing any common prefix and suffix from
  s1 and s2 and returning the prefix size to adjust the edit sequence

  >>> reduce_match('xxx','yyy')
  ('xxx', 'yyy', 0)
  >>> reduce_match('abc','ac')
  ('b', '', 1)
  >>> reduce_match('abcdefg','abcxxx')
  ('defg', 'xxx', 3)
  >>> reduce_match('abcxxx','defxxx')
  ('abc', 'def', 0)
  >>> reduce_match('abc','abc')
  ('', '', 3)
  '''
  cdef Py_ssize_t l1 = PyString_Size(s1)
  cdef Py_ssize_t l2 = PyString_Size(s2)

  if l1!=l2:
    raise ValueError('Length mismatch')

  cdef char *ss1 = s1
  cdef char *ss2 = s2

  cdef Py_ssize_t  i = 0, j = l1-1

  for i in range(l1):
    if ss1[i]!=ss2[i]:
      break

  for j in range(l1,i,-1):
    if ss1[j-1]!=ss2[j-1]:
      break

  return s1[i:j],s2[i:j],i


def levenshtein_distance(s1, s2):
  '''
  Calculate a minimum number of edit operations to transform s1 into s2
  based on the Levenshtein distance.

  This distance is the number of insertions, deletions, and substitutions
  needed to transform the first sequence into the second.

  See: http://en.wikipedia.org/wiki/Levenshtein_distance

       В.И. Левенштейн (1965). "Двоичные коды с исправлением выпадений, вставок и
       замещений символов".  Доклады Академий Наук СCCP 163 (4): 845–8.  Appeared
       in English as: Levenshtein VI (1966).  "Binary codes capable of correcting
       deletions, insertions, and reversals".  Soviet Physics Doklady 10:
       707–10.  http://www.scribd.com/full/18654513?access_key=key-10o99fv9kcwswflz9uas

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and O(min(N,M)) space, where N and M are the lengths
  of the two sequences.  See the following for more information on these
  concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
  cost[i,j] = |
              |     | cost[i-1,j-1]      if c1==c2     (match)
              | min | cost[i-1,j-1] + 1  if c1!=c2     (substitution)
              |     | cost[i,  j-1] + 1                (insert c2)
              |     | cost[i-1,j  ] + 1                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j].  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> levenshtein_distance('ac', 'abc')
  1
  >>> levenshtein_distance('ba', 'ab')
  2
  >>> levenshtein_distance('eat', 'seat')
  1
  >>> levenshtein_distance('abcdefghijk','cdefghijklm')
  4
  '''
  # Reduce problem size by removing any common prefix or suffix
  s1,s2,_ = reduce_match(s1,s2)

  # Minimize storage by forcing s2 to the shorter sequence
  if len(s1) < len(s2):
    s1,s2 = s2,s1

  # Fast-path: If s2 is empty, return distance of len(s1)
  if not s2:
    return len(s1)

  cdef Py_ssize_t n = PyString_Size(s1)
  cdef Py_ssize_t m = PyString_Size(s2)

  cdef Py_ssize_t  i,j,cost,match,insert,delete

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters

  cdef Py_ssize_t *current  = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  # Initialize previous costs to zero, since the values are not used and are
  # going to be overwritten
  for i in range(m+1):
    current[i]  = i
    previous[i] = 0

  cdef char *ss1 = s1
  cdef char *ss2 = s2
  cdef char c1,c2

  # For each location in s1
  for i in range(n):
    c1 = ss1[i]

    # Swap current and previous cost row storage to recycle the old
    # 'previous' as the new 'current'
    previous,current = current,previous

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    # Update minimum match, insertion, and deletion costs for each
    # location in s2
    for j in range(m):
      c2 = ss2[j]

      if c1==c2:
        cost = 0
      else:
        cost = 1

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match        = previous[j]   + cost

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert       = current[j]    + 1

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete       = previous[j+1] + 1

      # Take minimum cost operation
      current[j+1] = min3(match, insert, delete)

  # Return the minimum edit cost for both complete sequences.  i.e. the cost
  # transforming s1[0:i+1]->s2[0:j+1], which is current[-1].
  distance = current[m]

  free(current)
  free(previous)

  return distance


def damerau_levenshtein_distance(s1, s2):
  '''
  Calculate a minimum number of edit operations to transform s1 into s2
  based on the Damerau-Levenshtein distance.

  This distance is the number of insertions, deletions, substitutions, and
  transpositions needed to transform the first sequence into the second.
  Transpositions are exchanges of two *consecutive* characters and may not
  overlap with other transpositions.

  The operations required to incrementally transform s1 into s2 are returned
  as a sequence of operations, represented by tuples of the form:

    Substitution:  EditOp(op='S', pos=position, old=old,  new=new)
    Insertion:     EditOp(op='I', pos=position, old=None, new=new)
    Deletion:      EditOp(op='D', pos=position, old=old,  new=None)
    Transposition: EditOp(op='T', pos=position, old=old,  new=new)

  When compress=False, substitution, insertion, and deletions operations are
  always performed a single character or element at a time.  When
  compress=True, adjacent operations are combined into multi-character or
  multi-element operations, when possible.

  See: http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance

       В.И. Левенштейн (1965). "Двоичные коды с исправлением выпадений, вставок и
       замещений символов".  Доклады Академий Наук СCCP 163 (4): 845–8.  Appeared
       in English as: Levenshtein VI (1966).  "Binary codes capable of correcting
       deletions, insertions, and reversals".  Soviet Physics Doklady 10:
       707–10.  http://www.scribd.com/full/18654513?access_key=key-10o99fv9kcwswflz9uas

       Damerau F (1964). "A technique for computer detection and correction
       of spelling errors".  Communications of the ACM 7 (3):171-6.
       http://www.cis.uni-muenchen.de/~heller/SuchMasch/apcadg/literatur/data/damerau_distance.pdf

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and O(min(N,M)) space, where N and M are the lengths
  of the two sequences.  See the following for more information on these
  concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
              |
              |     | cost[i-1,j-1]      if c1==c2     (match)
  cost[i,j] = |     | cost[i-1,j-1] + 1  if c1!=c2     (substitution)
              |     | cost[i,  j-1] + 1                (insert c2)
              | min | cost[i-1,j  ] + 1                (delete c1)
              |     | cost[i-2,j-2] + 1  if i>1, j>1,  (transpose)
              |     |                       s1[i-2]==c2,
              |     |                       s2[j-2]==c1

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j].  This implementation saves space by only storing the last three
  cost rows at any given time (cost[i-2], cost[i-1], and cost[i]).

  >>> damerau_levenshtein_distance('ba', 'ab')
  1
  >>> damerau_levenshtein_distance('ba', 'abc')
  2
  >>> damerau_levenshtein_distance('fee', 'deed')
  2
  >>> damerau_levenshtein_distance('eat', 'seat')
  1
  '''
  # Reduce problem size by removing any common prefix or suffix
  s1,s2,_ = reduce_match(s1,s2)

  print '!!! reduced =',s2,s2

  # Minimize storage by forcing s2 to the shorter sequence
  if len(s1) < len(s2):
    s1,s2 = s2,s1

  # Fast-path: If s2 is empty, return distance of len(s1)
  if not s1:
    return len(s2)

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters

  cdef Py_ssize_t n = PyString_Size(s1)
  cdef Py_ssize_t m = PyString_Size(s2)

  cdef Py_ssize_t  i,j,cost,match,insert,delete,trans

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters

  cdef Py_ssize_t *current   = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous1 = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous2 = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  # Initialize previous two cost rows to zero, since the values are not used
  # and are going to be overwritten

  for i in range(m+1):
    current[i]   = i
    previous1[i] = 0
    previous2[i] = 0

  cdef char *ss1 = s1
  cdef char *ss2 = s2
  cdef char c1,c2

  # For each location in s1
  for i in range(n):
    c1 = ss1[i]

    # Swap current and previous cost row storage to recycle the old
    # 'previous2' as the new 'current'
    previous2,previous1,current = previous1,current,previous2

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    # Update minimum match, insertion, and deletion costs for each
    # location in s2
    for j in range(m):
      c2 = ss2[i]

      if c1==c2:
        cost = 0
      else:
        cost = 1

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match  = previous1[j]   + cost

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert = current[j]     + 1

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete = previous1[j+1] + 1

      # Transpose: transform s1[0:i-1]->s2[0:j-1] + 1,
      #            if s1[i]==s2[j-1] and s1[i-1]==s2[j]
      trans  = match
      if i and j and c1==s2[j-1] and s1[i-1]==c2:
        trans = previous2[j-1]+1

      # Take minimum cost operation
      current[j+1] = min4(match, insert, delete, trans)

  # Return the minimum edit cost for both complete sequences.  i.e. the cost
  # transforming s1[0:i+1]->s2[0:j+1], which is current[-1].

  distance = current[m]

  free(current)
  free(previous1)
  free(previous2)

  return distance
