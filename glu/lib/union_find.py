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
  A Class to implement union-find algorithm that performs two useful operations on a disjoint-set data structure:

    Find  : Determine which set a particular element is in and also if two elements are in the same set.
    Union : Combine or merge sets into a single set.
  '''
  def __init__(self, labels=None):
    '''
    Construct a new Union_find object with given sets of elements if any.
    '''
    if labels is not None:
      self.labels = dict( (l,l) for l in labels )
    else:
      self.labels = {}

  def find(self, i, create=True):
    '''
    Return the label for the specified element. If not existed, create a new label
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

  def has_key(self, i):
    '''
    Return whether or not the element that was passed in exists in the union_find object
    '''
    return i in self.labels

  def union(self, i1, i2):
    '''
    Merge sets into one set
    '''
    self.labels[self.find(i1)] = self.find(i2)

  def sets(self):
    '''
    Return the separate, nonoverlapping sets
    '''
    return self.setmap().values()

  def setmap(self):
    '''
    Return a dictionary which contains a number of separate, nonoverlapping sets
    '''
    sets = {}
    for i in self.labels:
      j = self.find(i)
      d = sets.get(j)
      if not d:
        d = sets[j] = set([j])
      d.add(i)
    return sets

  def __contains__(self,seq):
    '''Returns true if all items in seq belong to the same set'''
    seq = iter(seq)
    j = self.find(seq.next(),False)
    return all(j == self.find(i,False) for i in seq)


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
