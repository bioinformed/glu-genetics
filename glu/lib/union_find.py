# -*- coding: utf-8 -*-
'''
File:          union_find.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


class union_find(object):
  '''
  A union_find instance X maintains a family of disjoint sets of
  hashable objects, supporting the following two primary methods:

  - X[item] returns a name for the set containing the given item.
    Each set is named by an arbitrarily-chosen one of its members; as
    long as the set remains unchanged it will keep the same name. If
    the item is not yet part of a set in X, a new singleton set is
    created for it.

  - X.union(item1, item2, ...) merges the sets containing each item
    into a single larger set.  If any item is not yet part of a set
    in X, it is added to X as one of the members of the merged set.

  Based on code from Josiah Carlson with additional suggestions from D.
  Eppstein.  See:

      http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912

  >>> uf = union_find()
  >>> uf.union(1,2,3)
  >>> uf[1] is uf[2]
  True
  >>> uf[2] is uf[3]
  True

  >>> uf.union(4,5)
  >>> (4,5) in uf
  True
  >>> (1,4) in uf
  False

  >>> uf.union(3,5)
  >>> (1,4) in uf
  True

  >>> uf.union(10,11)
  >>> sorted( (ex,sorted(s)) for ex,s in uf.setmap().iteritems())
  [(1, [1, 2, 3, 4, 5]), (10, [10, 11])]
  '''
  def __init__(self, labels=None):
    '''
    Construct a new Union_find object with given sets of elements if any.
    '''
    self.labels = dict( (l,l) for l in labels or [])

  def find(self, i, create=True):
    '''
    Return the label for the specified element.  If create==False then the
    item is not tracked as part of the universe of possible objects.
    '''
    labels = self.labels

    # Assign new label
    if i not in labels:
      if create:
        labels[i] = i
      return i

    # Find exemplar
    j = i
    while labels[j] != j:
      j = labels[j]

    # Apply path compression
    while labels[i] != i:
      labels[i],i = j,labels[i]

    return j

  def __getitem__(self, i):
    '''
    Return the label for the specified element.
    '''
    return self.find(i)

  def union(self, *items):
    '''
    Union a sequence of items (len>1)
    '''
    if len(items) < 2:
      raise ValueError('Union requires at least two items')

    exemplars = (self.find(i) for i in items)
    exemplar  = exemplars.next()
    for i in exemplars:
      self.labels[i] = exemplar

  def has_key(self, i):
    '''
    Return whether or not the element that was passed in exists in the union_find object
    '''
    return i in self.labels

  def sets(self):
    '''
    Return the separate, nonoverlapping sets
    '''
    return self.setmap().values()

  def setmap(self):
    '''
    Return a dictionary of exemplars to the contents of each disjoint set
    '''
    sets = {}
    for i in self:
      j = self.find(i)
      d = sets.get(j)
      if not d:
        d = sets[j] = set([j])
      d.add(i)
    return sets

  def __contains__(self,seq):
    '''Returns true if all items in seq belong to the same set'''
    if len(seq) < 2:
      return True
    seq = iter(seq)
    j = self.find(seq.next(),False)
    return all(j == self.find(i,False) for i in seq)

  def __iter__(self):
    '''Iterator over all items seen'''
    return iter(self.labels)


def build_dupsets(dups):
  '''
  Given a set of elements that were passed in, break them up into a number of separate, nonoverlapping sets.

  @param dups: set of elements
  @type  dups: sequence
  @return    : a number of separate, nonoverlapping sets
  @rtype     : list of sets

  >>> dups = [(1,1),(2,2),(3,3),(4,4),(1,2),(3,4)]
  >>> build_dupsets(dups)
  [set([1, 2]), set([3, 4])]
  >>> dups = [('a','b'),('a','c'),('a','a'),('d','e'),('d','f')]
  >>> build_dupsets(dups)
  [set(['a', 'c', 'b']), set(['e', 'd', 'f'])]
  '''
  uf = union_find()
  for i1,i2 in dups:
    uf.union(i1,i2)
  return uf.sets()


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
