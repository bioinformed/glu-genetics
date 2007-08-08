# -*- coding: utf-8 -*-
'''
File:          remap.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import unittest
from   itertools   import izip,chain,islice


def _remap_comp(a1s,a2s):
  return sum(1 for a1,a2 in izip(a1s,a2s) if a1 and a1==a2 )


def remap_category(allelemap):
  # FIXME: Use genoarray model once support is available
  complement = {'A':'T','T':'A','G':'C','C':'G'}
  a1s,a2s    = zip(*( (a1,a2) for a1,a2 in allelemap.iteritems() if a1 and a2 ))
  identity   = _remap_comp(a1s,a2s)
  complement = _remap_comp(a1s, (complement.get(a2) for a2 in a2s))
  swap       = _remap_comp(a1s, (allelemap.get(a2)  for a2 in a2s))

  n = len(allelemap)

  if identity == n:
    return 'identity'
  elif complement == n:
    return 'complement'
  elif not identity and swap == n or (identity,swap,complement,n) == (0,0,0,1):
    return 'swap'
  elif complement+identity == n:
    return 'partial complement'
  elif swap+identity > n:
    return 'partial swap'
  else:
    return 'other'


def eval_remap(gmap, genocounts):
  concord = 0
  for (g1,g2),n in genocounts.iteritems():
    a     = g1
    b1,b2 = g2
    b1,b2 = gmap.get(b1),gmap.get(b2)
    if a == (b1,b2) or a == (b2,b1):
      concord += n
  return concord


def remap_alleles(genocounts):
  '''
  Find the best(i.e. the highest concordance) allele mapping between two sets of genotypes

  @param:  genocounts: key is a tuple of left geno and right geno, value is the count of the pair
  @type:   genocounts: dictionary
  @return: the mapping from the right alleles to the left alleles
  @rtype:  dictionary
  '''
  genos = izip(*genocounts)
  left  = set(a for g in genos.next() for a in g)
  right = set(a for g in genos.next() for a in g)

  if len(left) > len(right):
    right.update( left - right )
  elif len(left) < len(right):
    left.update( right - left )

  left  = list(left)
  right = list(right)

  maps = (dict(izip(right,lperm)) for lperm in permutations(left))
  concord,bestmap = max( (eval_remap(map,genocounts),map) for map in maps )

  return concord,bestmap


def permutations(l):
  n = len(l)
  if not n:
    yield []
  elif n == 1:
    yield [l[0]]
  else:
    a,b = l[:1],l[1:]
    for p in permutations(b):
      for i in xrange(n):
        yield tuple(chain(islice(p,0,i),a,islice(p,i,None)))


def encode(s):
  for left,right,count in s:
    yield (tuple(left),tuple(right)),count


class TestRemap(unittest.TestCase):
  def test_remap_alleles(self):
    cases = [([('AA','TT',20),('AG','TC',10),('GG','CC',30),('AG','TT',50)], (60,{'T':'A','C':'G'})),
             ([('AA','AA',20),('AG','AG',10),('GG','GG',30),('AG','AA',50)], (60,{'A':'A','G':'G'})),
             ([('GA','GT',20),('AA','TT',10),('GG','GG',30)],                (60,{'T':'A','G':'G'})),
             ([('AT','GC',20),('AA','GG',10),('TT','CC',30)],                (60,{'C':'T','G':'A'})),
             ([('AT','AC',20),('AA','AA',10),('TT','CC',30)],                (60,{'A':'A','C':'T'})),
             ([('AA','TT',20),('GG','AA',10),('AG','AT',30)],                (60,{'A':'G','T':'A'})),
             ([('AA','AT',30),('AA','TT',40),('AA','AA',30)],                (40,{'A':'T','T':'A'})),
             ([('AA','AT',30),('AA','TT',20),('AA','AA',30)],                (30,{'A':'A','T':'T'})),
             ([('AT','AA',30),('TT','AA',40),('AA','AA',30)],                (40,{'A':'T','T':'A'})),
             ([('AT','AG',30),('TT','AT',40),('AA','GT',30),('AA','AA',10)], (40,{'A':'A','T':'G','G':'T'})),
             ([('AT','GC',30),('AA','AG',30),('TT','TG',20),('GT','AG',20)], (50,{'A':'G','C':'A','T':'C','G':'T'})),
             ([('AA','CC', 6),('AC','CC', 1)],                               ( 6,{'A':'C','C':'A'}))]
    for data,result in cases:
      self.assertEquals(remap_alleles(dict(encode(data))),result)

  def test_remap_category(self):
    cases = [({'A':'A','C':'C'},         'identity'),
             ({'A':'T','T':'A'},         'complement'),
             ({'A':'T','C':'C'},         'partial complement'),
             ({'A':'C','C':'A'},         'swap'),
             ({'A':'C'},                 'swap'),
             ({'G':'C'},                 'complement'),
             ({'A':'C','C':'A','G':'G'}, 'partial swap'),
             ({'A':'T','T':'A','G':'G'}, 'partial complement'),
             ({'A':'C','C':'T'},         'other'),
             ({'A':'G','G':'A','T':'T'}, 'partial swap'),
             ({'A':'G','G':'A','T':'T'}, 'partial swap'),
             ({'A':'A','G':'C'},         'partial complement'),
             ({'A':'G','G':'A','G':'C'}, 'other'),
             ({'A':'A','G':'G'},         'identity'),
             ({'A':'T','G':'C'},         'complement'),
             ({'A':'G','G':'A'},         'swap')]

    for data,result in cases:
      self.assertEquals(remap_category(data), result)


if __name__=='__main__':
  unittest.main()
