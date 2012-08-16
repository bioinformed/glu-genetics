# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Various edit distance algorithms'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   itertools         import izip, groupby

import numpy as np

from   glu.lib.utils     import izip_exact, namedtuple


EditOp  = namedtuple('EditOp',  'op pos old new')
CigarOp = namedtuple('CigarOp', 'op count')


def _compress_edit_sequence(seq):
  '''
  Compress an edit sequence by merging adjacent character substitutions,
  insertion and deletion operations into single multi-character operations.

  >>> Op = EditOp
  >>> list(_compress_edit_sequence([Op('S',5,'A','G'),Op('S',6,'A','T'),Op('S',7,'A','G')]))
  [EditOp(op='S', pos=5, old='AAA', new='GTG')]
  >>> list(_compress_edit_sequence([Op('S',5,'A','G'),Op('S',7,'A','T'),Op('S',9,'A','G')]))
  [EditOp(op='S', pos=5, old='A', new='G'), EditOp(op='S', pos=7, old='A', new='T'), EditOp(op='S', pos=9, old='A', new='G')]
  >>> list(_compress_edit_sequence([Op('D',5,'A',None),Op('D',5,'B',None),Op('D',5,'C',None)]))
  [EditOp(op='D', pos=5, old='ABC', new=None)]
  >>> list(_compress_edit_sequence([Op('D',5,'A',None),Op('D',6,'B',None),Op('D',8,'C',None)]))
  [EditOp(op='D', pos=5, old='A', new=None), EditOp(op='D', pos=6, old='B', new=None), EditOp(op='D', pos=8, old='C', new=None)]
  >>> list(_compress_edit_sequence([Op('I',5,None,'A'),Op('I',6,None,'B'),Op('I',7,None,'C')]))
  [EditOp(op='I', pos=5, old=None, new='ABC')]
  >>> list(_compress_edit_sequence([Op('I',5,None,'A'),Op('I',7,None,'B'),Op('I',8,None,'C')]))
  [EditOp(op='I', pos=5, old=None, new='A'), EditOp(op='I', pos=7, old=None, new='BC')]
  '''
  last = None

  for edit in seq:
    if not last:
      last = edit
    elif edit.op==last.op=='D' and edit.pos==last.pos:
      last = EditOp(last.op,last.pos,last.old+edit.old,None)
    elif edit.op==last.op=='S' and edit.pos==last.pos+len(last.old):
      last = EditOp(last.op,last.pos,last.old+edit.old,last.new+edit.new)
    elif edit.op==last.op=='I' and edit.pos==last.pos+len(last.new):
      last = EditOp(last.op,last.pos,None,last.new+edit.new)
    else:
      yield last
      last = edit

  if last:
    yield last


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
  return sum(ch1 != ch2 for ch1, ch2 in izip_exact(s1, s2))


def hamming_sequence(s1, s2, compress=True):
  '''
  Calculate a minimum sequence of edit operations based on the
  Hamming distance between sequences.

  This distance is based on the number of substitutions needed to transform
  the first sequence into the second.  Both sequences are required to be of
  equal length.

  The operations required to incrementally transform s1 into s2 are returned
  as a sequence of operations, represented by tuples of the form:

    Substitution:  EditOp(op='S', pos=position, old=old, new=new)

  When compress=False, operations are always performed a single character or
  element at a time.  When compress=True, adjacent operations are combined
  into multi-character or multi-element operations, when possible.

  See: http://en.wikipedia.org/wiki/Hamming_distance

       Hamming, Richard W. (1950), "Error detecting and error correcting
       codes", Bell System Technical Journal 26 (2): 147–160,
       http://www.ece.rutgers.edu/~bushnell/dsdwebsite/hamming.pdf

  This implementation requires O(N) time and space, where N is the length of
  the input sequences.

  >>> hamming_sequence('abc', 'abc')
  []
  >>> hamming_sequence('abb', 'abc')
  [EditOp(op='S', pos=2, old='b', new='c')]
  >>> hamming_sequence('abc', 'def')
  [EditOp(op='S', pos=0, old='abc', new='def')]
  >>> hamming_sequence('abc', 'def', compress=False)
  [EditOp(op='S', pos=0, old='a', new='d'), EditOp(op='S', pos=1, old='b', new='e'), EditOp(op='S', pos=2, old='c', new='f')]
  >>> hamming_sequence('a', 'ab')
  Traceback (most recent call last):
  ...
  LengthMismatch
  '''
  seq = [ EditOp('S',i,c1,c2) for i,(c1,c2) in enumerate(izip_exact(s1, s2)) if c1!=c2 ]
  if compress:
    seq = list(_compress_edit_sequence(seq))
  return seq


def _edit_sequence(s1, s2, edits, offset=0):
  '''
  Compute the sequence of edits required to transform sequence s1 to s2
  using the operations encoded in the supplied matrix of edit operations.
  Used internally for edit matrices generated by levenshtein_sequence and
  damerau_levenshtein_sequence
  '''
  i   = len(s1)-1
  j   = len(s2)-1
  seq = []

  while i>=0 and j>=0:
    c1 = s1[i]
    c2 = s2[j]
    op = edits[i,j]

    if op=='M':
      if c1!=c2:
        seq.append( EditOp('S',offset+j,c1,c2) )
      i -= 1
      j -= 1
    elif op=='I':
      seq.append( EditOp('I',offset+j,None,c2) )
      j -= 1
    elif op=='D':
      seq.append( EditOp('D',offset+j,c1,None) )
      i -= 1
    elif op=='T':
      seq.append( EditOp('T',offset+j-1,s1[i-1:i+1],s2[j-1:j+1]) )
      i -= 2
      j -= 2
    else:
      raise ValueError('Invalid edit operation')

  while j>=0:
    seq.append( EditOp('I',offset+j,None,s2[j]) )
    j -= 1

  while i>=0:
    seq.append( EditOp('D',offset,s1[i],None) )
    i -= 1

  seq.reverse()
  return seq


def edit_alignment(s1,s2,edit_seq):
  a1 = []
  a2 = []

  for edit in edit_seq:
    old = edit.old or ('-'*len(edit.new))
    new = edit.new or ('-'*len(edit.old))
    a1.extend(old)
    a2.extend(new)

  return ''.join(a1),''.join(a2)


def reduce_match(s1,s2):
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
  i=j=0
  n = min(len(s1),len(s2))

  while i<n and s1[i]==s2[i]:
    i+=1

  s1 = s1[i:][::-1]
  s2 = s2[i:][::-1]

  n -= i
  while j<n and s1[j]==s2[j]:
    j+=1

  s1 = s1[j:][::-1]
  s2 = s2[j:][::-1]

  return s1,s2,i


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

  It works with arbitrary sequences too:
  >>> levenshtein_distance('abcd', ['b', 'a', 'c', 'd', 'e'])
  3
  '''
  # Reduce problem size by removing any common prefix or suffix
  s1,s2,_ = reduce_match(s1,s2)

  # Minimize storage by forcing s2 to the shorter sequence
  if len(s1) < len(s2):
    s1,s2 = s2,s1

  # Fast-path: If s2 is empty, return distance of len(s1)
  if not s2:
    return len(s1)

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters
  m        = len(s2)

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  current  = range(m+1)

  # Initialize previous costs to zero, since the values are not used and are
  # going to be overwritten
  previous = [0]*(m+1)

  # For each location in s1
  for i,c1 in enumerate(s1):
    # Swap current and previous cost row storage to recycle the old
    # 'previous' as the new 'current'
    previous,current = current,previous

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    # Update minimum match, insertion, and deletion costs for each
    # location in s2
    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match        = previous[j]   + (c1 != c2)

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert       = current[j]    + 1

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete       = previous[j+1] + 1

      # Take minimum cost operation
      current[j+1] = min(match, insert, delete)

  # Return the minimum edit cost for both complete sequences.  i.e. the cost
  # transforming s1[0:i+1]->s2[0:j+1], which is current[-1].
  return current[-1]


def levenshtein_sequence(s1, s2, compress=True):
  '''
  Calculate a minimum number of edit operations to transform s1 into s2
  based on the Levenshtein distance.

  This sequence is comprised of the minimum number insertions, deletions,
  and substitutions needed to transform the first sequence into the second.

  The operations required to incrementally transform s1 into s2 are returned
  as a sequence of operations, represented by tuples of the form:

    Substitution:  EditOp(op='S', pos=position, old=old,  new=new)
    Insertion:     EditOp(op='I', pos=position, old=None, new=new)
    Deletion:      EditOp(op='D', pos=position, old=old,  new=None)

  When compress=False, substitution, insertion, and deletions operations are
  always performed a single character or element at a time.  When
  compress=True, adjacent operations are combined into multi-character or
  multi-element operations, when possible.

  See: http://en.wikipedia.org/wiki/Levenshtein_distance

       В.И. Левенштейн (1965). "Двоичные коды с исправлением выпадений, вставок и
       замещений символов".  Доклады Академий Наук СCCP 163 (4): 845–8.  Appeared
       in English as: Levenshtein VI (1966).  "Binary codes capable of correcting
       deletions, insertions, and reversals".  Soviet Physics Doklady 10:
       707–10.  http://www.scribd.com/full/18654513?access_key=key-10o99fv9kcwswflz9uas

  This function is related to levenshtein_sequence since the length of the
  uncompressed edit sequence is equal to the minimum edit distance:

    len(levenshtein_sequence(s1, s2, False)) == levenshtein_distance(s1,s2)

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
              |     | cost[i-1,j-1]      if c1==c2     (match)
              | min | cost[i-1,j-1] + 1  if c1!=c2     (substitution)
              |     | cost[i,  j-1] + 1                (insert c2)
              |     | cost[i-1,j  ] + 1                (delete c1)

  where c1=s1[i-1], c2=s2[j-1].  The resulting minimum edit distance is then
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last two
  cost rows at any given time (cost[i-1], and cost[i]).

  >>> levenshtein_sequence('', '')
  []
  >>> levenshtein_sequence('abc', 'abc')
  []
  >>> levenshtein_sequence('', 'abc', compress=False)
  [EditOp(op='I', pos=0, old=None, new='a'), EditOp(op='I', pos=1, old=None, new='b'), EditOp(op='I', pos=2, old=None, new='c')]
  >>> levenshtein_sequence('', 'abc')
  [EditOp(op='I', pos=0, old=None, new='abc')]
  >>> levenshtein_sequence('abc', '')
  [EditOp(op='D', pos=0, old='abc', new=None)]
  >>> levenshtein_sequence('abc', '', compress=False)
  [EditOp(op='D', pos=0, old='a', new=None), EditOp(op='D', pos=0, old='b', new=None), EditOp(op='D', pos=0, old='c', new=None)]
  >>> levenshtein_sequence('ba', 'ab')
  [EditOp(op='S', pos=0, old='ba', new='ab')]
  >>> levenshtein_sequence('ba', 'ab', compress=False)
  [EditOp(op='S', pos=0, old='b', new='a'), EditOp(op='S', pos=1, old='a', new='b')]
  >>> levenshtein_sequence('eat', 'seat')
  [EditOp(op='I', pos=0, old=None, new='s')]
  >>> levenshtein_sequence('abcdefghijk','cdefghijklm')
  [EditOp(op='D', pos=0, old='ab', new=None), EditOp(op='I', pos=9, old=None, new='lm')]
  >>> levenshtein_sequence('abcdefghijk','cdefghijklm', compress=False)
  [EditOp(op='D', pos=0, old='a', new=None), EditOp(op='D', pos=0, old='b', new=None), EditOp(op='I', pos=9, old=None, new='l'), EditOp(op='I', pos=10, old=None, new='m')]

  It works with arbitrary sequences too:
  >>> levenshtein_sequence('abcd', ['b', 'a', 'c', 'd', 'e'])
  [EditOp(op='S', pos=0, old='ab', new='ba'), EditOp(op='I', pos=4, old=None, new='e')]
  '''
  # Reduce problem size by removing any common prefix or suffix, noting the
  # prefix offset to adjust the edit sequence
  s1,s2,offset = reduce_match(s1,s2)

  # Fast-path: If both s1 and s2 are is empty, report a perfect match
  if not s1 and not s2:
    return []

  # Fast-path: If s1 is empty, then construct a series of insertions to create s2
  if not s1:
    return [ EditOp('I',offset,  None,s2) ] if compress \
      else [ EditOp('I',offset+j,None,c2) for j,c2 in enumerate(s2) ]

  # Fast-path: If s2 is empty, then construct a series of deletions to remove all of s1
  if not s2:
    return [ EditOp('D',offset,s1,None) ] if compress \
      else [ EditOp('D',offset,c1,None) for c1 in s1 ]

  # Otherwise, prepare storage for the edit costs and operations Only store
  # the current and previous cost rows, adding an element to the beginning
  # for the cost of inserting i characters
  n        = len(s1)
  m        = len(s2)

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  current  = range(m+1)

  # Initialize previous costs to zero, since the values are not used and are
  # going to be overwritten
  previous = [0]*(m+1)

  # Allocate an empty character matrix to track the best edits at each step
  # in order to reconstruct an optimal sequence at the end
  edits    = np.zeros((n,m), dtype='S1')

  # For each location in s1
  for i,c1 in enumerate(s1):
    # Swap current and previous cost row storage to recycle the old
    # 'previous' as the new 'current'
    previous,current = current,previous

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match        = previous[j]   + (c1 != c2)

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert       = current[j]    + 1

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete       = previous[j+1] + 1

      # Take minimum cost operation
      current[j+1] = mcost = min(match, insert, delete)

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if mcost==match:
        edits[i,j]='M'
      elif mcost==insert:
        edits[i,j]='I'
      else:
        edits[i,j]='D'

  # Build and return a minimal edit sequence using the saved operations
  seq = _edit_sequence(s1, s2, edits, offset)

  # Compress sequential substitution, deletion, and insertion operations
  if compress:
    seq = list(_compress_edit_sequence(seq))

  return seq


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

  It works with arbitrary sequences too:
  >>> damerau_levenshtein_distance('abcd', ['b', 'a', 'c', 'd', 'e'])
  2
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
  m         = len(s2)

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  current   = range(m+1)

  # Initialize previous two cost rows to zero, since the values are not used
  # and are going to be overwritten
  previous1 = [0]*(m+1)
  previous2 = [0]*(m+1)

  # For each location in s1
  for i,c1 in enumerate(s1):
    # Swap current and previous cost row storage to recycle the old
    # 'previous2' as the new 'current'
    previous2,previous1,current = previous1,current,previous2

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    # Update minimum match, insertion, and deletion costs for each
    # location in s2
    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match  = previous1[j]   + (c1 != c2)

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
      current[j+1] = min(match, insert, delete, trans)

  # Return the minimum edit cost for both complete sequences.  i.e. the cost
  # transforming s1[0:i+1]->s2[0:j+1], which is current[-1].
  return current[-1]


def damerau_levenshtein_sequence(s1, s2, compress=True):
  '''
  Calculate a minimum sequence of edit operations to transform s1 into s2
  based on the Damerau-Levenshtein distance.

  This sequence is comprised of the minimum number of insertions, deletions,
  substitutions, and transpositions needed to transform the first sequence
  into the second.  Transpositions are exchanges of two *consecutive*
  characters and may not overlap with other transpositions.

  See: http://en.wikipedia.org/wiki/Damerau%E2%80%93Levenshtein_distance

       В.И. Левенштейн (1965). "Двоичные коды с исправлением выпадений, вставок и
       замещений символов".  Доклады Академий Наук СCCP 163 (4): 845–8.  Appeared
       in English as: Levenshtein VI (1966).  "Binary codes capable of correcting
       deletions, insertions, and reversals".  Soviet Physics Doklady 10:
       707–10.  http://www.scribd.com/full/18654513?access_key=key-10o99fv9kcwswflz9uas

       Damerau F (1964). "A technique for computer detection and correction
       of spelling errors".  Communications of the ACM 7 (3):171-6.
       http://www.cis.uni-muenchen.de/~heller/SuchMasch/apcadg/literatur/data/damerau_distance.pdf

  This function is related to damerau_levenshtein_sequence since the length
  of the uncompressed edit sequence is equal to the minimum edit distance:

    len(damerau_levenshtein_sequence(s1, s2, False))
                               == damerau_levenshtein_distance(s1,s2)

  This implementation is based on a standard dynamic programming algorithm,
  requiring O(N*M) time and space, where N and M are the lengths of the two
  sequences.  See the following for more information on these concepts:

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
  cost[i,j] and the edit sequence is obtained by keeping note of which
  operation was selected at each step and backtracking from the end to the
  beginning.  This implementation saves space by only storing the last three
  cost rows at any given time (cost[i-2], cost[i-1], and cost[i]).

  >>> damerau_levenshtein_sequence('ba', 'ab')
  [EditOp(op='T', pos=0, old='ba', new='ab')]
  >>> damerau_levenshtein_sequence('ba', 'abc')
  [EditOp(op='I', pos=0, old=None, new='a'), EditOp(op='S', pos=2, old='a', new='c')]
  >>> damerau_levenshtein_sequence('fee', 'deed')
  [EditOp(op='S', pos=0, old='f', new='d'), EditOp(op='I', pos=3, old=None, new='d')]
  >>> damerau_levenshtein_sequence('eat', 'seat')
  [EditOp(op='I', pos=0, old=None, new='s')]

  It works with arbitrary sequences too:
  >>> damerau_levenshtein_sequence('abcd', ['b', 'a', 'c', 'd', 'e'])
  [EditOp(op='T', pos=0, old='ab', new=['b', 'a']), EditOp(op='I', pos=4, old=None, new='e')]
  '''
  # Reduce problem size by removing any common prefix or suffix, noting the
  # prefix offset to adjust the edit sequence
  s1,s2,offset = reduce_match(s1,s2)

  # Fast-path: If both s1 and s2 are is empty, report a perfect match
  if not s1 and not s2:
    return []

  # Fast-path: If s1 is empty, then construct a series of insertions to create s2
  if not s1:
    return [ EditOp('I',offset,  None,s2) ] if compress \
      else [ EditOp('I',offset+j,None,c2) for j,c2 in enumerate(s2) ]

  # Fast-path: If s2 is empty, then construct a series of deletions to remove all of s1
  if not s2:
    return [ EditOp('D',offset,s1,None) ] if compress \
      else [ EditOp('D',offset,c1,None) for c1 in s1 ]

  # Otherwise, prepare storage for the edit costs
  # Only store the current and previous cost rows, adding a single
  # element to the beginning for the cost of inserting i characters
  n         = len(s1)
  m         = len(s2)

  # Start by assigning the cost of transforming an empty s1 into s2[0:m] by
  # inserting 0..m elements
  current   = range(m+1)

  # Initialize previous two cost rows to zero, since the values are not used
  # and are going to be overwritten
  previous1 = [0]*(m+1)
  previous2 = [0]*(m+1)

  # Allocate an empty character matrix to track the best edits at each step
  # in order to reconstruct an optimal sequence at the end
  edits     = np.zeros((n,m), dtype='S1')

  # For each location in s1
  for i,c1 in enumerate(s1):
    # Swap current and previous cost row storage to recycle the old
    # 'previous2' as the new 'current'
    previous2,previous1,current = previous1,current,previous2

    # Initialize the cost of inserting characters up to and including i
    current[0] = i+1

    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match:        transform s1[0:i]->s2[0:j] + 0, if s1[i]==s2[j]
      # Substitution: transform s1[0:i]->s2[0:j] + 1, if s1[i]!=s2[j]
      match        = previous1[j]   + (c1 != c2)

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert       = current[j]     + 1

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete       = previous1[j+1] + 1

      # Transpose: transform s1[0:i-1]->s2[0:j-1] + 1,
      #            if s1[i]==s2[j-1] and s1[i-1]==s2[j]
      trans        = match
      if i and j and c1==s2[j-1] and s1[i-1]==c2:
        trans      = previous2[j-1]+1

      # Take minimum cost operation
      current[j+1] = mcost = min(match, insert, delete, trans)

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions over transpositions.  This ambiguity for
      # equal cost options implies that there may not be a unique optimum
      # edit sequence, but one or more sequences of equal length.
      if mcost==match:
        edits[i,j]='M'
      elif mcost==insert:
        edits[i,j]='I'
      elif mcost==delete:
        edits[i,j]='D'
      else:
        edits[i,j]='T'

  # Build and return a minimal edit sequence using the saved operations
  seq = _edit_sequence(s1, s2, edits, offset)

  # Compress sequential substitution, deletion, and insertion operations
  if compress:
    seq = list(_compress_edit_sequence(seq))

  return seq


def _roll_cigar(i, j, edits, anchor_left):
  '''
  Compute the sequence of edits required to transform sequence s1 to s2
  using the operations encoded in the supplied matrix of edit operations.
  Used internally for edit matrices generated by levenshtein_sequence and
  damerau_levenshtein_sequence
  '''
  n,m = edits.shape

  igar = []

  while i>=0 and j>=0:
    op = edits[i,j]

    if not op:
      break

    igar.append(op)
    if op in 'M=X':
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

  if anchor_left:
    if j>=0:
      igar.extend('I'*(j+1))
      j = -1

    if i>=0:
      igar.extend('D'*(i+1))
      i = -1

  cigar = [ CigarOp(code,len(list(run))) for code,run in groupby(reversed(igar)) ]

  return i+1,j+1,cigar


def cigar_to_string(cigar):
  return ''.join('%d%s' % (n,op) for op,n in cigar)


def cigar_alignment(s1,s2,cigar,hide_match=True):
  i  = j = 0
  a1 = []
  a2 = []

  for op,n in cigar:
    if op in 'M=X':
      a1.extend(s1[i:i+n])

      if not hide_match:
        a2.extend(s2[j:j+n])
      else:
        a2.extend( ('.' if c1==c2 else c2) for c1,c2 in izip(s1[i:i+n],s2[j:j+n]))

      i += n
      j += n

    elif op=='I':
      a1.extend('-'*n)
      a2.extend(s2[j:j+n])
      j += n

    elif op=='N':
      a1.extend(' '*n)
      a2.extend(s2[j:j+n])
      j += n

    elif op=='D':
      a1.extend(s1[i:i+n])
      a2.extend('-'*n)
      i += n

    elif op=='S':
      a1.extend(s1[i:i+n])
      a2.extend(' '*n)
      i += n

    else:
      raise ValueError('Unknown CIGAR operation')

  return ''.join(a1),''.join(a2)


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
  if local_align is not None and not local_align:
    anchor_left = anchor_right = True
  elif anchor_left or anchor_right:
    local_align = False

  # Prepare storage for the edit costs and operations
  # Allocate an empty character matrix to track the best edits at each step
  # in order to reconstruct an optimal sequence at the end
  n         = len(s1)
  m         = len(s2)

  # Transform scalar scores into position-specific scores
  if isinstance(match_score,int):
    match_score = [match_score]*n

  if isinstance(mismatch_score,int):
    mismatch_score = [mismatch_score]*n

  if isinstance(gap_score,int):
    gap_score = [gap_score]*n

  curr_cost = range(m+1) if anchor_left else [0]*(m+1)
  prev_cost = [0]*(m+1)

  edits     = np.zeros((n,m), dtype='S1')

  end_i     = 0
  end_j     = 0
  score     = 0

  for i,c1 in enumerate(s1):

    curr_cost,prev_cost = prev_cost,curr_cost

    curr_cost[0] = i+1 if anchor_left else 0

    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      match  = prev_cost[j] + (match_score[i] if c1==c2 else mismatch_score[i])

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert = curr_cost[j] + gap_score[i]

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete = prev_cost[j+1] + gap_score[i]

      # Take best costing operation
      mcost = max(match, insert, delete)

      if not anchor_left:
        mcost = max(0, mcost)

      curr_cost[j+1] = mcost

      if mcost>score:
        end_i = i
        end_j = j
        score = mcost

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if not anchor_left and mcost==0:
        pass
      elif mcost==match:
        edits[i,j]='=' if c1==c2 else 'X'
      elif mcost==insert:
        edits[i,j]='I'
      else:
        edits[i,j]='D'

  # Build and return a minimal edit sequence using the saved operations
  if anchor_right:
    end_i = n-1
    end_j = m-1
    score = mcost

  start_i,start_j,cigar = _roll_cigar(end_i,end_j,edits, anchor_left)

  return slice(start_i,end_i+1),slice(start_j,end_j+1),score,cigar


def smith_waterman_gotoh(s1, s2, match_score=10, mismatch_score=-9,
                                 gap_open_score=-15, gap_extend_score=-6,
                                 s2_homopolymer_penalty=None,
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
  if gap_open_score==gap_extend_score:
    return smith_waterman(s1, s2, match_score=match_score, mismatch_score=mismatch_score,
                                  gap_score=gap_open_score, anchor_left=anchor_left,
                                  anchor_right=anchor_right, local_align=local_align)

  if local_align is not None and not local_align:
    anchor_left = anchor_right = True

  # Prepare storage for the edit costs and operations
  # Allocate an empty character matrix to track the best edits at each step
  # in order to reconstruct an optimal sequence at the end
  n         = len(s1)
  m         = len(s2)

  # Transform scalar scores into position-specific scores
  if isinstance(match_score,int):
    match_score = [match_score]*n

  if isinstance(mismatch_score,int):
    mismatch_score = [mismatch_score]*n

  if isinstance(gap_open_score,int):
    gap_open_score = [gap_open_score]*n

  if isinstance(gap_extend_score,int):
    gap_extend_score = [gap_extend_score]*n

  curr_cost = [0]*(m+1)
  prev_cost = [0]*(m+1)
  curr_gap1 = [0]*(m+1)  # Not needed
  prev_gap1 = [0]*(m+1)  # Can use a scalar
  curr_gap2 = [0]*(m+1)
  prev_gap2 = [0]*(m+1)  # Not needed, since curr_gap2 does not overlap

  if anchor_left:
    curr_cost[0] = curr_gap1[0] = curr_gap2[0] = 0
    curr_cost[1] = curr_gap1[1] = curr_gap2[1] = gap_open_score[0]
    for i in range(1,m):
      curr_cost[i+1] = curr_gap1[i+1] = curr_gap2[i+1] = curr_cost[i] + gap_extend_score[0]

  edits     = np.zeros((n,m), dtype='S1')

  end_i     = 0
  end_j     = 0
  score     = 0

  for i,c1 in enumerate(s1):

    curr_cost,prev_cost = prev_cost,curr_cost
    curr_gap1,prev_gap1 = prev_gap1,curr_gap1
    curr_gap2,prev_gap2 = prev_gap2,curr_gap2

    curr_cost[0] = prev_cost[0]+gap_extend_score[i] if anchor_left else 0
    curr_gap1[0] = prev_gap1[0]+gap_extend_score[i] if anchor_left else 0
    curr_gap2[0] = prev_gap2[0]+gap_extend_score[i] if anchor_left else 0

    for j,c2 in enumerate(s2):
      # Compute cost of transforming s1[0:i+1]->s2[0:j+1] allowing the
      # following edit operations:

      # Match/Mismatch: transform s1[0:i]->s2[0:j] + match/mismatch cost
      match  = prev_cost[j] + (match_score[i] if c1==c2 else mismatch_score[i])

      # Insert: transform s1[0:i+1]->s2[0:j] and insert s2[j]
      insert = curr_gap1[j+1] = max(curr_gap1[j] + gap_extend_score[i],
                                    curr_cost[j] + gap_open_score[i])

      # Delete: transform s1[0:i]->s2[0:j+1] and delete s1[i]
      delete = curr_gap2[j+1] = max(prev_gap2[j+1] + gap_extend_score[i],
                                    prev_cost[j+1] + gap_open_score[i])

      # Take best costing operation
      mcost = max(match, insert, delete)

      if not anchor_left:
        mcost = max(0, mcost)

      curr_cost[j+1] = mcost

      if mcost>score:
        end_i = i
        end_j = j
        score = mcost

      # Record the operation chosen, with preference for (mis)matches over
      # insertions over deletions.  This ambiguity for equal cost options
      # implies that there may not be a unique optimum edit sequence, but
      # one or more sequences of equal length.
      if not anchor_left and mcost==0:
        pass
      elif mcost==match:
        edits[i,j]='=' if c1==c2 else 'X'
      elif mcost==insert:
        edits[i,j]='I'
      else:
        edits[i,j]='D'

  # Build and return a minimal edit sequence using the saved operations
  if anchor_right:
    end_i = n-1
    end_j = m-1
    score = mcost

  start_i,start_j,cigar = _roll_cigar(end_i,end_j,edits,anchor_left)

  return slice(start_i,end_i+1),slice(start_j,end_j+1),score,cigar


try:
  from glu.lib.seqlib._edits import *
except ImportError:
  pass


def _benchmark():
  import time

  s1='AGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAGAGACCAAGTCTCTGCTACCGTACATACTCGTACTGAGACTGCCAAGGCACACAGGGGATAG'
  s2='GCTGGTGCGACACATGCTGGTGCGACACGTGCGACACATGCTGGGGTGCGACACATGCTGGTGCGACACGTGCGACACATGCTGGTGCGACACAT'

  N = 20000
  t0 = time.time()
  for i in range(N):
    smith_waterman_gotoh(s1,s2,mismatch_score=-20)
  t1 = time.time()

  print 'Standard:  %8d aligns/s' % (N/(t1-t0))

  t0 = time.time()
  for i in range(N):
    smith_waterman_gotoh2_score(s1,s2,mismatch_score=-20)
  t1 = time.time()

  print 'New score: %8d aligns/s' % (N/(t1-t0))

  t0 = time.time()
  for i in range(N):
    smith_waterman_gotoh2_align(s1,s2,mismatch_score=-20)
  t1 = time.time()

  print 'New align: %8d aligns/s' % (N/(t1-t0))


def _test():
  import doctest
  doctest.testmod()
  try:
    import glu.lib.seqlib._edits as _edits
    doctest.testmod(_edits)
  except ImportError:
    pass


if __name__ == '__main__':
  _test()
