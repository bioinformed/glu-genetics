# -*- coding: utf-8 -*-
'''
File:          dataset.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from itertools import izip, repeat
from random    import shuffle


class Locus(object):
  def __init__(self, name, data, chromosome=None, location=None):
    self.name       = name
    self.chromosome = chromosome
    self.location   = location


class Sample(object):
  def __init__(self, name):
    self.name = name


class GenotypeData(object):
  def __init__(self, locusmajor=True):
    # Assume locus-major for now
    self._genos     = []
    self.locusmajor = locusmajor
    self.locusmap   = {}
    self.locusnames = []
    self.samplemap  = {}
    self.samplenames= []

  def locus_index(self, name):
    n = self.locusmap.get(name,None)
    if n is None:
      n = len(self.locusmap)
      self.locusmap[name] = n
      self.locusnames.append(name)
    return n

  def sample_index(self, name):
    n = self.samplemap.get(name,None)
    if n is None:
      n = len(self.samplemap)
      self.samplemap[name] = n
      self.samplenames.append(name)
    return n

  def has_locus(self, name):
    return name in self.locusnames

  def has_sample(self, name):
    return name in self.samplenames

  def loci(self):
    return (LocusData(name,self) for name in self.locusnames)

  def locus(self,name):
    self.locusmap[name]
    return LocusData(name, self)

  def samples(self):
    return (SampleData(name,self) for name in self.samplenames)

  def sample(self,name):
    self.samplemap[name]
    return SampleData(name, self)

  def transpose(self):
    newdata = GenotypeData(not self.locusmajor)
    newdata.locusmap = self.locusmap.copy()
    newdata.samplemap = self.samplemap.copy()
    if self.locusmajor:
      for sample in self.samples():
        newdata.sample(sample.name).set_genolist(sample.genolist())
    else:
      for locus in self.loci():
        newdata.locus(locus.name).set_genolist(locus.genolist())

    return newdata

  def genos(self):
    if self.locusmajor:
      for locus in self.loci():
        for g in locus.genos():
          yield g
    else:
      for sample in self.samples():
        for g in self.sample.genos():
          yield g

  def set_genos(self, genos):
    selfgenos = self._genos

    for locus,sample,geno in genos:
      i=self.locus_index(locus)
      j=self.sample_index(sample)

      if not self.locusmajor:
        i,j = j,i

      n = len(selfgenos)
      if i >= n:
        selfgenos.extend( [] for i in xrange(i-n+1) )

      g = selfgenos[i]

      m = len(g)
      if j >= m:
        g.extend( ['  ']*(j-m+1) )

      g[j] = geno

  def _set_genolist(self, index, rowmajor, genolist):
    n = len(self._genos)
    if rowmajor:
      if index == n:
        self._genos.append(genolist)
      elif index < n:
        self._genos[index] = genolist
      else:
        self._genos.extend( [] for i in xrange(index-n) )
        self._genos.append(genolist)
    else:
      m = len(genolist)
      if n <= m:
        self._genos.extend( [] for i in xrange(m-n) )

      assert len(self._genos) == len(genolist)

      for gs,g in izip(self._genos,genolist):
        n = len(gs)
        if index == n:
          gs.append(g)
        elif index < n:
          gs[index] = g
        else:
          gs.extend(['  ']*(index-n+1))
          gs[-1] = g


class LocusData(object):
  def __init__(self, name, data):
    self.name = name
    self.data = data

  def set_genolist(self, genos):
    data = self.data
    genos = list(genos)

    assert len(genos) == len(data.samplemap)

    i = data.locusmap[self.name]
    data._set_genolist(i, data.locusmajor, genos)

  def genolist(self):
    data = self.data
    i = data.locusmap[self.name]
    if data.locusmajor:
      return data._genos[i]
    else:
      return [ g[i] for g in data._genos ]

  def genos(self):
    return izip(repeat(self.name),self.data.samplenames,self.genolist())


class SampleData(object):
  def __init__(self, name, data):
    self.name = name
    self.data = data

  def set_genolist(self, genos):
    data = self.data
    genos = list(genos)

    assert len(genos) == len(data.locusmap)

    i = data.samplemap[self.name]
    data._set_genolist(i, not data.locusmajor, genos)

  def genolist(self):
    data = self.data
    i = data.samplemap[self.name]
    if not data.locusmajor:
      return data._genos[i]
    else:
      return [ g[i] for g in data._genos ]

  def genos(self):
    return izip(self.data.locusnames,repeat(self.name),self.genolist())



def completion(data):
  for item in data:
    m,n = genotype_completion(item.genolist())
    yield item.name(),m,n


def completion_by_locus(data):
  return completion(data.loci())


def completion_by_sample(data):
  return completion(data.samples())


def completion_all(data):
  samples = defaultdict(lambda: [0,0])
  loci    = defaultdict(lambda: [0,0])
  for locus,sample,geno in data.genos():
    s = samples[sample]
    l = loci[locus]
    s[1] += 1
    l[1] += 1
    if not s.missing():
      s[0] += 1
      l[0] += 1
  return dict(locus),dict(samples)


def concordance(data1,data2):
  samples1 = list(data1.samples())
  samples2 = list(data2.samples())
  names1 = set(sample.name for sample in samples1)
  names2 = set(sample.name for sample in samples2)
  common = names1&names2

  for name in common:
    m,n = genotype_concordance(samples1[name],samples2[name])
    yield name,m,n


def dupcheck(data):
  items = list(data)
  n = len(items)
  for i in range(0,n):
    for j in range(i,n):
      a = items[i]
      b = items[j]
      m,n = genotype_concordance(a.genotypes(),b.genotypes())
      if is_duplicate(m,n):
        yield a.name(),b.name(),m,n


def main():
  loci    = [ Locus('rs1', '1', 12345), Locus('rs2', '2', 23456), Locus('rs3', '3', 34567) ]
  samples = [ Sample('s1'), Sample('s2'), Sample('s3') ]
  genos   = [ ['AA', 'AB', 'BB'], ['  ','CD','CD'], ['EE','FF',' E'] ]
  data = GenotypeData()

  def init(data):
    for sample in samples:
      data.sample_index(sample.name)
    for locus in loci:
      data.locus_index(locus.name)

  init(data)
  lgs = zip(loci,genos)*2
  shuffle(lgs)
  for locus,gs in lgs:
    data.locus(locus.name).set_genolist(gs)

  assert data._genos == [['AA','AB','BB'],['  ','CD','CD'],['EE','FF',' E']]

  xloci = data.loci()
  locus = xloci.next()
  assert locus.name == 'rs1' and locus.genolist() == ['AA','AB','BB']
  locus = xloci.next()
  assert locus.name == 'rs2' and locus.genolist() == ['  ','CD','CD']
  locus = xloci.next()
  assert locus.name == 'rs3' and locus.genolist() == ['EE','FF',' E']

  xsamples = data.samples()
  sample = xsamples.next()
  assert sample.name == 's1' and sample.genolist() == ['AA', '  ', 'EE']
  sample = xsamples.next()
  assert sample.name == 's2' and sample.genolist() == ['AB', 'CD', 'FF']
  sample = xsamples.next()
  assert sample.name == 's3' and sample.genolist() ==  ['BB', 'CD', ' E']

  assert list(data.genos()) == [('rs1', 's1', 'AA'), ('rs1', 's2', 'AB'), ('rs1', 's3', 'BB'),
                                ('rs2', 's1', '  '), ('rs2', 's2', 'CD'), ('rs2', 's3', 'CD'),
                                ('rs3', 's1', 'EE'), ('rs3', 's2', 'FF'), ('rs3', 's3', ' E')]

  assert data.transpose()._genos == [['AA', '  ', 'EE'], ['AB', 'CD', 'FF'], ['BB', 'CD', ' E']]
  datat = data.transpose()

  d2 = GenotypeData()
  init(d2)
  g = list(data.genos())
  shuffle(g*2)
  d2.set_genos(g)

  assert data._genos == d2._genos

  d3 = GenotypeData(False)
  init(d3)
  g = list(data.genos())
  shuffle(g*2)
  d3.set_genos(g)
  d3 = d3.transpose()
  assert data._genos == d3._genos

  d4 = GenotypeData()
  init(d4)

  for sample in data.samples():
    d4.sample(sample.name).set_genolist(sample.genolist())

  assert d4._genos == data._genos

  d5 = GenotypeData(False)
  init(d5)

  for sample in data.samples():
    d5.sample(sample.name).set_genolist(sample.genolist())

  assert d5._genos == datat._genos

  print 'All tests passed'


if __name__ == '__main__':
  main()
