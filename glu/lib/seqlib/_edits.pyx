# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Various edit distance algorithms'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

# C imports
cimport cython
cimport numpy as np

from    libc.stdio          cimport sprintf
from    libc.limits         cimport INT_MIN
from    libc.stdlib         cimport malloc, realloc, calloc, free
from    cpython             cimport PyString_Size, PyErr_NoMemory

# Python imports
import  numpy as np

from    glu.lib.utils        import namedtuple
from    glu.lib.seqlib.edits import cigar_to_string, cigar_alignment


np.import_array()


# Alignment operation data structures
EditOp  = namedtuple('EditOp',  'op pos old new')

CigarOp = namedtuple('CigarOp', 'op count')


cdef inline Py_ssize_t min2(Py_ssize_t a, Py_ssize_t b):
  if a < b:
    return a
  else:
    return b


cdef inline Py_ssize_t min3(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c):
  return min2(min2(a,b),c)


cdef inline Py_ssize_t min4(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c, Py_ssize_t d):
  return min2(min2(a,b),min2(c,d))


cdef inline Py_ssize_t max2(Py_ssize_t a, Py_ssize_t b):
  if a >= b:
    return a
  else:
    return b


cdef inline Py_ssize_t max3(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c):
  return max2(max2(a,b),c)


cdef inline Py_ssize_t max4(Py_ssize_t a, Py_ssize_t b, Py_ssize_t c, Py_ssize_t d):
  return max2(max2(a,b),max2(c,d))


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
  ValueError: Length mismatch
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
  cdef Py_ssize_t l  = min2(l1,l2)

  cdef char *ss1 = s1
  cdef char *ss2 = s2

  cdef Py_ssize_t  i = 0, j = l1-1

  for i in range(l):
    if ss1[i]!=ss2[i]:
      break

  for j in range(l-i):
    if ss1[l1-j-1]!=ss2[l2-j-1]:
      break

  return s1[i:l1-j],s2[i:l2-j],i


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

  cdef Py_ssize_t *current  = <Py_ssize_t *>calloc( (m+1), sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous = <Py_ssize_t *>calloc( (m+1), sizeof(Py_ssize_t))

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
  cdef char *ss1    = s1
  cdef char *ss2    = s2

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters

  cdef Py_ssize_t *current   = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous1 = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))
  cdef Py_ssize_t *previous2 = <Py_ssize_t *>malloc( (m+1)*sizeof(Py_ssize_t))

  cdef Py_ssize_t  i,j,cost,match,insert,delete,trans
  cdef char c1,c2

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  # Initialize previous two cost rows to zero, since the values are not used
  # and are going to be overwritten

  for i in range(m+1):
    current[i]   = i
    previous1[i] = 0
    previous2[i] = 0

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
      c2 = ss2[j]

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
      if i and j and c1==ss2[j-1] and ss1[i-1]==c2:
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


cdef inline _op_str(char op):
  if op=='S':
    opstr = 'S'
  elif op=='N':
    opstr = 'N'
  elif op=='M':
    opstr = 'M'
  elif op=='=':
    opstr = '='
  elif op=='X':
    opstr = 'X'
  elif op=='I':
    opstr = 'I'
  elif op=='D':
    opstr = 'D'
  elif op=='T':
    opstr = 'T'
  else:
    raise ValueError('Unknown op=%s' % op)
  return opstr


cdef _roll_cigar(Py_ssize_t n, Py_ssize_t m, Py_ssize_t i, Py_ssize_t j, char *edits, int anchor_left):
  '''
  Compute the sequence of edits required to transform sequence s1 to s2
  using the operations encoded in the supplied matrix of edit operations.
  Used internally for edit matrices generated by levenshtein_sequence and
  damerau_levenshtein_sequence
  '''
  cdef Py_ssize_t l,count
  cdef char op,last_op

  cdef list cigar = []

  last_op = 0
  count   = 0

  while i>=0 and j>=0:
    op = edits[m*i+j]

    if op==0:
      break

    if last_op==op:
      count += 1
    else:
      if count:
        cigar.append( CigarOp(_op_str(last_op),count) )

      last_op = op
      count   = 1

    if op=='M' or op=='=' or op=='X':
      i -= 1
      j -= 1
    elif op=='D':
      i -= 1
    elif op=='I':
      j -= 1
    elif op=='T':
      i -= 2
      j -= 2
    else:
      raise ValueError('Invalid edit operation')

  if count:
    cigar.append( CigarOp(_op_str(last_op),count) )

  if anchor_left:
    if j>=0:
      cigar.append( CigarOp('I', j+1) )
      j = -1

    if i>=0:
      cigar.append( CigarOp('D', i+1) )
      i = -1

  cigar.reverse()

  return i+1,j+1,cigar


@cython.boundscheck(False)
def smith_waterman(s1, s2, match_score=10, mismatch_score=-9, gap_score=-12,
                           anchor_left=False,anchor_right=False,local_align=None):
  '''
  Align s1 to s2 using the Smith-Waterman algorithm for local gapped
  alignment.  An alignment score and sequence of alignment operations are returned.

  The operations to align s1 to s2 are returned as a sequence, represented
  by extended CIGAR (Compact Idiosyncratic Gapped Alignment Report)
  operations of the form:

    Match:         CigarOp(op='=', count=n)
    Mismatch:      CigarOp(op='X', count=n)
    Insertion:     CigarOp(op='I', count=n)
    Deletion:      CigarOp(op='D', count=n)

  Match operations are inclusive of matching and mismatch characters

  See: http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

       Smith, Temple F.; and Waterman, Michael S. (1981). "Identification
       of Common Molecular Subsequences". Journal of Molecular Biology 147: 195–197.
       http://gel.ym.edu.tw/~chc/AB_papers/03.pdf

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and space, where N and M are the lengths of the two
  sequences.  See the following for more information on these concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
  cost[i,j] = |
              |     | 0
              |     | cost[i-1,j-1] + m   if c1==c2     (match: perfect)
              | min | cost[i-1,j-1] + mm  if c1!=c2     (match: substitution)
              |     | cost[i,  j-1] + g                (insert c2)
              |     | cost[i-1,j  ] + g                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> s1,s2='b','abc'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2)
  >>> p1
  slice(0, 1, None)
  >>> p2
  slice(1, 2, None)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abc','b'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2)
  >>> p1
  slice(1, 2, None)
  >>> p2
  slice(0, 1, None)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abcbd','acd'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2,match_score=20)
  >>> p1
  slice(0, 5, None)
  >>> p2
  slice(0, 3, None)
  >>> score
  36
  >>> cigar_to_string(cigar)
  '1=1D1=1D1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'abcbd'
  '.-.-.'

  >>> s2,s1='abcbd','acd'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2,match_score=20)
  >>> p1
  slice(0, 3, None)
  >>> p2
  slice(0, 5, None)
  >>> score
  36
  >>> cigar_to_string(cigar)
  '1=1I1=1I1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'a-c-d'
  '.b.b.'

  >>> s1,s2='abcbd','beb'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2)
  >>> p1
  slice(1, 4, None)
  >>> p2
  slice(0, 3, None)
  >>> score
  11
  >>> cigar_to_string(cigar)
  '1=1X1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'bcb'
  '.e.'

  >>> s1,s2='abcdefghijklmnop','pqrstuvqw'
  >>> p1,p2,score,cigar = smith_waterman(s1,s2)
  >>> p1
  slice(15, 16, None)
  >>> p2
  slice(0, 1, None)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'p'
  '.'
  '''
  cdef int _local_align, _anchor_left, _anchor_right

  _local_align  = 1 if local_align  else 0
  _anchor_left  = 1 if anchor_left  else 0
  _anchor_right = 1 if anchor_right else 0

  if local_align is not None and not local_align:
    _local_align = _anchor_left = _anchor_right = 1
  elif _anchor_left or _anchor_right:
    _local_align  = 0

  # Prepare storage for the edit costs and operations
  # Allocate an empty character matrix to track the best edits at each step
  # in order to reconstruct an optimal sequence at the end

  cdef Py_ssize_t n    = PyString_Size(s1)
  cdef Py_ssize_t m    = PyString_Size(s2)
  cdef char *ss1       = s1
  cdef char *ss2       = s2

  cdef np.ndarray[int, ndim=1, mode="c"]    match_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"] mismatch_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]      gap_scores = np.empty(n, dtype=np.int32)

  match_scores[:]      = match_score
  mismatch_scores[:]   = mismatch_score
  gap_scores[:]        = gap_score

  cdef char *edits     = <char*>malloc(n*m*sizeof(char))

  cdef int  *curr_cost =  <int*>malloc((m+1)*sizeof(int))
  cdef int  *prev_cost =  <int*>malloc((m+1)*sizeof(int))

  cdef Py_ssize_t i,j,start_i,end_i,start_j,end_j,score,mcost,match,insert,delete
  cdef char c1,c2,op

  end_i = end_j = score = mcost = 0

  for j in range(m+1):
    curr_cost[j] = j if _anchor_left else 0

  for i in range(n):
    c1 = ss1[i]

    curr_cost,prev_cost = prev_cost,curr_cost

    curr_cost[0] = i+1 if _anchor_left else 0

    for j in range(m):
      c2 = ss2[j]

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      match  = prev_cost[j]

      if c1==c2:
        match += match_scores[i]
      else:
        match += mismatch_scores[i]

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert = curr_cost[j] + gap_scores[i]

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete = prev_cost[j+1] + gap_scores[i]

      # Take best costing operation
      if not _anchor_left:
        curr_cost[j+1] = mcost = max4(0, match, insert, delete)
      else:
        curr_cost[j+1] = mcost = max3(match, insert, delete)

      if mcost>score:
         score = mcost
         end_i = i
         end_j = j

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if not _anchor_left and mcost==0:
        op = 0
      elif mcost==match and c1==c2:
        op = '='
      elif mcost==match:
        op = 'X'
      elif mcost==insert:
        op = 'I'
      else:
        op = 'D'

      edits[m*i+j]=op

  # Build and return a minimal edit sequence using the saved operations
  if _anchor_right:
    end_i = n-1
    end_j = m-1
    score = mcost

  start_i,start_j,cigar = _roll_cigar(n,m,end_i,end_j,edits, _anchor_left)

  free(prev_cost)
  free(curr_cost)

  free(edits)

  return slice(start_i,end_i+1),slice(start_j,end_j+1),score,cigar


@cython.boundscheck(False)
def smith_waterman_gotoh(s1, s2, match_score=10, mismatch_score=-9,
                                 gap_open_score=-15, gap_extend_score=-6,
                                 anchor_left=False,anchor_right=False,local_align=None):
  '''
  Align s1 to s2 using the Smith-Waterman algorithm for local gapped
  alignment.  An alignment score and sequence of alignment operations are returned.

  The operations to align s1 to s2 are returned as a sequence, represented
  by extended CIGAR (Compact Idiosyncratic Gapped Alignment Report)
  operations of the form:

    Match:         CigarOp(op='=', count=n)
    Mismatch:      CigarOp(op='X', count=n)
    Insertion:     CigarOp(op='I', count=n)
    Deletion:      CigarOp(op='D', count=n)

  Match operations are inclusive of matching and mismatch characters

  See: http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

       Smith, Temple F.; and Waterman, Michael S. (1981). "Identification
       of Common Molecular Subsequences". Journal of Molecular Biology 147: 195–197.
       http://gel.ym.edu.tw/~chc/AB_papers/03.pdf

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and space, where N and M are the lengths of the two
  sequences.  See the following for more information on these concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
  cost[i,j] = |
              |     | 0
              |     | cost[i-1,j-1] + m   if c1==c2     (match: perfect)
              | min | cost[i-1,j-1] + mm  if c1!=c2     (match: substitution)
              |     | cost[i,  j-1] + g                (insert c2)
              |     | cost[i-1,j  ] + g                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> s1,s2='b','abc'
  >>> p1,p2,score,cigar = smith_waterman_gotoh(s1,s2)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abc','b'
  >>> p1,p2,score,cigar = smith_waterman_gotoh(s1,s2)
  >>> p1
  slice(1, 2, None)
  >>> p2
  slice(0, 1, None)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abbcbbd','acd'
  >>> p1,p2,score,cigar = smith_waterman_gotoh(s1,s2,match_score=30)
  >>> p1
  slice(0, 7, None)
  >>> p2
  slice(0, 3, None)
  >>> score
  48
  >>> cigar_to_string(cigar)
  '1=2D1=2D1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'abbcbbd'
  '.--.--.'

  >>> s1,s2='cc','accct'
  >>> p1,p2,score,cigar = smith_waterman_gotoh(s1,s2,match_score=1,mismatch_score=-1,gap_open_score=-4,gap_extend_score=-1,local_align=False)
  >>> p1
  slice(0, 2, None)
  >>> p2
  slice(0, 5, None)
  >>> score
  -6
  >>> cigar_to_string(cigar)
  '3I1=1X'
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  '---cc'
  'acc.t'

  >>> s1='AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
  >>> s2='GCTGGTGCGACACAT'
  >>> p1,p2,score,cigar = smith_waterman_gotoh(s1,s2,mismatch_score=-20)
  >>> p1
  slice(46, 53, None)
  >>> p2
  slice(6, 14, None)
  >>> score
  55
  >>> cigar_to_string(cigar)
  '2=1I5='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'GC-ACACA'
  '..G.....'
  '''
  # Fall back to the standard Smith Waterman when the overhead of the Gotoh
  # scoring is not needed
  #if gap_open_score==gap_extend_score:
  #  return smith_waterman(s1, s2, match_score=match_score, mismatch_score=mismatch_score,
  #                                gap_score=gap_open_score)

  cdef int _local_align, _anchor_left, _anchor_right

  _local_align  = 1 if local_align  else 0
  _anchor_left  = 1 if anchor_left  else 0
  _anchor_right = 1 if anchor_right else 0

  if local_align is not None and not local_align:
    _local_align = _anchor_left = _anchor_right = 1
  elif _anchor_left or _anchor_right:
    _local_align  = 0

  cdef Py_ssize_t n    = PyString_Size(s1)
  cdef Py_ssize_t m    = PyString_Size(s2)
  cdef char *ss1       = s1
  cdef char *ss2       = s2

  cdef np.ndarray[int, ndim=1, mode="c"]      match_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   mismatch_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   gap_open_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"] gap_extend_scores = np.empty(n, dtype=np.int32)

  match_scores[:]      = match_score
  mismatch_scores[:]   = mismatch_score
  gap_open_scores[:]   = gap_open_score
  gap_extend_scores[:] = gap_extend_score

  cdef char *edits     = <char*>malloc(n*m*sizeof(char))

  cdef int  *curr_cost =  <int*>malloc((m+1)*sizeof(int))
  cdef int  *prev_cost =  <int*>malloc((m+1)*sizeof(int))
  cdef int  *curr_gap1 =  <int*>malloc((m+1)*sizeof(int))
  cdef int  *curr_gap2 =  <int*>malloc((m+1)*sizeof(int))
  cdef int  *prev_gap2 =  <int*>malloc((m+1)*sizeof(int))

  cdef int local_match, local_mismatch, local_gap_open, local_gap_extend

  cdef Py_ssize_t i,j,start_i,end_i,start_j,end_j,score,mcost,match,insert,delete
  cdef char c1,c2,op

  end_i = end_j = score = mcost = 0

  if _anchor_left:
    curr_cost[0] = curr_gap1[0] = curr_gap2[0] = 0
    curr_cost[1] = curr_gap1[1] = curr_gap2[1] = gap_open_scores[0]
    for j in range(1,m):
      curr_cost[j+1] = curr_gap1[j+1] = curr_gap2[j+1] = curr_cost[j] + gap_extend_scores[0]
  else:
    for j in range(m+1):
      curr_cost[j] = curr_gap1[j] = curr_gap2[j] = 0

  for i in range(n):
    local_match      = match_scores[i]
    local_mismatch   = mismatch_scores[i]
    local_gap_open   = gap_open_scores[i]
    local_gap_extend = gap_extend_scores[i]

    c1 = ss1[i]

    curr_cost,prev_cost = prev_cost,curr_cost
    curr_gap2,prev_gap2 = prev_gap2,curr_gap2

    curr_cost[0] = prev_cost[0]+local_gap_extend if _anchor_left else 0
    curr_gap1[0] = curr_gap1[0]+local_gap_extend if _anchor_left else 0
    curr_gap2[0] = prev_gap2[0]+local_gap_extend if _anchor_left else 0

    for j in range(m):
      c2 = ss2[j]

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      match = prev_cost[j]

      if c1==c2:
        match += local_match
      else:
        match += local_mismatch

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert = curr_gap1[j+1] = max2(curr_gap1[j] + local_gap_extend,
                                     curr_cost[j] + local_gap_open)

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete = curr_gap2[j+1] = max2(prev_gap2[j+1] + local_gap_extend,
                                     prev_cost[j+1] + local_gap_open)

      # Take best costing operation
      if not _anchor_left:
        curr_cost[j+1] = mcost = max4(0, match, insert, delete)
      else:
        curr_cost[j+1] = mcost = max3(match, insert, delete)

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if mcost>score:
         score = mcost
         end_i = i
         end_j = j

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if not _anchor_left and mcost==0:
        op = 0
      elif mcost==match and c1==c2:
        op = '='
      elif mcost==match:
        op = 'X'
      elif mcost==insert:
        op = 'I'
      else:
        op = 'D'

      edits[m*i+j]=op

  # Build and return a minimal edit sequence using the saved operations
  if _anchor_right:
    end_i = n-1
    end_j = m-1
    score = mcost

  start_i,start_j,cigar = _roll_cigar(n,m,end_i,end_j,edits, _anchor_left)

  free(prev_gap2)
  free(curr_gap2)
  free(curr_gap1)
  free(prev_cost)
  free(curr_cost)

  free(edits)

  return slice(start_i,end_i+1),slice(start_j,end_j+1),score,cigar


@cython.boundscheck(False)
def smith_waterman_gotoh2_align(s1, s2, match_score=10, mismatch_score=-9,
                                        gap_open_score=-15, gap_extend_score=-6):
  '''
  Align s1 to s2 using the Smith-Waterman algorithm for local gapped
  alignment.  An alignment score and sequence of alignment operations are returned.

  The operations to align s1 to s2 are returned as a sequence, represented
  by extended CIGAR (Compact Idiosyncratic Gapped Alignment Report)
  operations of the form:

    Match:         CigarOp(op='=', count=n)
    Mismatch:      CigarOp(op='X', count=n)
    Insertion:     CigarOp(op='I', count=n)
    Deletion:      CigarOp(op='D', count=n)

  Match operations are inclusive of matching and mismatch characters

  See: http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

       Smith, Temple F.; and Waterman, Michael S. (1981). "Identification
       of Common Molecular Subsequences". Journal of Molecular Biology 147: 195–197.
       http://gel.ym.edu.tw/~chc/AB_papers/03.pdf

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and space, where N and M are the lengths of the two
  sequences.  See the following for more information on these concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
  cost[i,j] = |
              |     | 0
              |     | cost[i-1,j-1] + m   if c1==c2     (match: perfect)
              | min | cost[i-1,j-1] + mm  if c1!=c2     (match: substitution)
              |     | cost[i,  j-1] + g                (insert c2)
              |     | cost[i-1,j  ] + g                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> s1,s2='b','abc'
  >>> p1,p2,score,cigar = smith_waterman_gotoh2_align(s1,s2)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abc','b'
  >>> p1,p2,score,cigar = smith_waterman_gotoh2_align(s1,s2)
  >>> p1
  slice(1, 2, None)
  >>> p2
  slice(0, 1, None)
  >>> score
  10
  >>> cigar_to_string(cigar)
  '1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'b'
  '.'

  >>> s1,s2='abbcbbd','acd'
  >>> p1,p2,score,cigar = smith_waterman_gotoh2_align(s1,s2,match_score=30)
  >>> p1
  slice(0, 7, None)
  >>> p2
  slice(0, 3, None)
  >>> score
  48
  >>> cigar_to_string(cigar)
  '1=2D1=2D1='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'abbcbbd'
  '.--.--.'

  >>> s1='AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
  >>> s2='GCTGGTGCGACACAT'
  >>> p1,p2,score,cigar = smith_waterman_gotoh2_align(s1,s2,mismatch_score=-20)
  >>> p1
  slice(46, 53, None)
  >>> p2
  slice(6, 14, None)
  >>> score
  55
  >>> cigar_to_string(cigar)
  '2=1I5='
  >>> a1,a2 = cigar_alignment(s1[p1],s2[p2],cigar)
  >>> print "'%s'\\n'%s'" % (a1,a2) # doctest: +NORMALIZE_WHITESPACE
  'GC-ACACA'
  '..G.....'
  '''
  cdef Py_ssize_t n    = PyString_Size(s1)
  cdef Py_ssize_t m    = PyString_Size(s2)
  cdef char *ss1       = s1
  cdef char *ss2       = s2

  cdef np.ndarray[int, ndim=1, mode="c"]      match_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   mismatch_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   gap_open_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"] gap_extend_scores = np.empty(n, dtype=np.int32)

  match_scores[:]      = match_score
  mismatch_scores[:]   = mismatch_score
  gap_open_scores[:]   = gap_open_score
  gap_extend_scores[:] = gap_extend_score

  cdef char *edits     = <char*>malloc(n*m*sizeof(char))
  cdef int  *S         =  <int*>malloc(m*sizeof(int))
  cdef int  *coldel    =  <int*>malloc(m*sizeof(int))
  cdef int  NINF       = INT_MIN+100000
  cdef int  Sup, Sleft, Sdiag, Stemp1, Stemp2, Sij, Smax, inscost, delcost

  cdef int local_match, local_mismatch, local_gap_open, local_gap_extend

  cdef unsigned int i,j,start_i,end_i,start_j,end_j
  cdef char c1,c2,op,*editpos=edits

  Sij = Smax = end_i = end_j = 0

  for j in range(m):
     coldel[j] = NINF
     S[j]      = 0

  for i in range(n):
    local_match      = match_scores[i]
    local_mismatch   = mismatch_scores[i]
    local_gap_open   = gap_open_scores[i]
    local_gap_extend = gap_extend_scores[i]

    Sleft   = 0
    Sdiag   = 0
    inscost = NINF
    c1      = ss1[i]

    for j in range(m):
      c2  = ss2[j]
      Sup = S[j]

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delcost  = coldel[j] + local_gap_extend
      Stemp1   = Sup       + local_gap_open
      inscost +=             local_gap_extend
      Stemp2   = Sleft     + local_gap_open

      if delcost < Stemp1:
        delcost = Stemp1
      if inscost < Stemp2:
        inscost = Stemp2

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      if c1==c2:
        Sij = Sdiag + local_match
        op  = '='
      else:
        Sij = Sdiag + local_mismatch
        op  = 'X'

      # Take best costing operation
      if Sij < delcost:
        Sij = delcost
        op  = 'D'

      if Sij < inscost:
        Sij = inscost
        op  = 'I'

      if Sij < 0:
        Sij = 0
        op  = 0

      if Sij > Smax:
        Smax  = Sij
        end_i = i
        end_j = j

      Sdiag        = Sup
      S[j] = Sleft = Sij
      coldel[j]    = delcost
      editpos[0]   = op
      editpos     += 1

  start_i,start_j,cigar = _roll_cigar(n,m,end_i,end_j,edits,0)

  free(coldel)
  free(S)
  free(edits)

  return slice(<int>start_i,<int>end_i+1),slice(<int>start_j,<int>end_j+1),<int>Smax,cigar


@cython.boundscheck(False)
def smith_waterman_gotoh2_score(s1, s2, match_score=10, mismatch_score=-9,
                                        gap_open_score=-15, gap_extend_score=-6):
  '''
  Align s1 to s2 using the Smith-Waterman algorithm for local gapped
  alignment.  An alignment score and sequence of alignment operations are returned.

  The operations to align s1 to s2 are returned as a sequence, represented
  by extended CIGAR (Compact Idiosyncratic Gapped Alignment Report)
  operations of the form:

    Match:         CigarOp(op='=', count=n)
    Mismatch:      CigarOp(op='X', count=n)
    Insertion:     CigarOp(op='I', count=n)
    Deletion:      CigarOp(op='D', count=n)

  Match operations are inclusive of matching and mismatch characters

  See: http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm

       Smith, Temple F.; and Waterman, Michael S. (1981). "Identification
       of Common Molecular Subsequences". Journal of Molecular Biology 147: 195–197.
       http://gel.ym.edu.tw/~chc/AB_papers/03.pdf

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and space, where N and M are the lengths of the two
  sequences.  See the following for more information on these concepts:

    http://en.wikipedia.org/wiki/Dynamic_programming
    http://en.wikipedia.org/wiki/Big_O_notation

  The the cost to transform s1[:i]->s2[:j] is based on the following
  recurrence:

              |  i                       if i>=0, j=0  (delete s1[:i])
              |  j                       if i =0, j>0  (insert s2[:j])
  cost[i,j] = |
              |     | 0
              |     | cost[i-1,j-1] + m   if c1==c2     (match: perfect)
              | min | cost[i-1,j-1] + mm  if c1!=c2     (match: substitution)
              |     | cost[i,  j-1] + g                (insert c2)
              |     | cost[i-1,j  ] + g                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> s1,s2='b','abc'
  >>> smith_waterman_gotoh2_score(s1,s2)
  10

  >>> s1,s2='abc','b'
  >>> smith_waterman_gotoh2_score(s1,s2)
  10

  >>> s1,s2='abbcbbd','acd'
  >>> smith_waterman_gotoh2_score(s1,s2,match_score=30)
  48

  >>> s1='AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
  >>> s2='GCTGGTGCGACACAT'
  >>> smith_waterman_gotoh2_score(s1,s2,mismatch_score=-20)
  55
  '''
  cdef Py_ssize_t n    = PyString_Size(s1)
  cdef Py_ssize_t m    = PyString_Size(s2)
  cdef char *ss1       = s1
  cdef char *ss2       = s2

  cdef np.ndarray[int, ndim=1, mode="c"]      match_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   mismatch_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"]   gap_open_scores = np.empty(n, dtype=np.int32)
  cdef np.ndarray[int, ndim=1, mode="c"] gap_extend_scores = np.empty(n, dtype=np.int32)

  match_scores[:]      = match_score
  mismatch_scores[:]   = mismatch_score
  gap_open_scores[:]   = gap_open_score
  gap_extend_scores[:] = gap_extend_score

  cdef int  *S         =  <int*>malloc(m*sizeof(int))
  cdef int  *coldel    =  <int*>malloc(m*sizeof(int))
  cdef int  NINF       = INT_MIN+100000
  cdef int  Sup, Sleft, Sdiag, Stemp1, Stemp2, Sij, Smax, inscost, delcost

  cdef int local_match, local_mismatch, local_gap_open, local_gap_extend

  cdef unsigned int i,j,end_i,end_j
  cdef char c1,c2

  Sij = end_i = end_j = Smax = 0

  for j in range(m):
     coldel[j] = NINF
     S[j]      = 0

  for i in range(n):
    local_match      = match_scores[i]
    local_mismatch   = mismatch_scores[i]
    local_gap_open   = gap_open_scores[i]
    local_gap_extend = gap_extend_scores[i]

    Sleft   = 0
    Sdiag   = 0
    inscost = NINF
    c1      = ss1[i]

    for j in range(m):
      c2  = ss2[j]
      Sup = S[j]

      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      delcost  = coldel[j] + local_gap_extend
      Stemp1   = Sup       + local_gap_open
      inscost +=             local_gap_extend
      Stemp2   = Sleft     + local_gap_open

      if delcost < Stemp1:
        delcost = Stemp1
      if inscost < Stemp2:
        inscost = Stemp2

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      if c1==c2:
        Sij = Sdiag + local_match
      else:
        Sij = Sdiag + local_mismatch

      # Take best costing operation
      if Sij < delcost:
        Sij = delcost

      if Sij < inscost:
        Sij = inscost

      if Sij < 0:
        Sij = 0

      if Sij > Smax:
        end_i = i
        end_j = j
        Smax  = Sij

      S[j] = Sleft = Sij
      Sdiag        = Sup
      coldel[j]    = delcost

  free(coldel)
  free(S)

  return <int>end_i,<int>end_j,<int>Smax
