# -*- coding: utf-8 -*-
'''
File:          genoarray2.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-04-10

Abstract:      Efficient bit-packed genotype array representation

Requires:      Python 2.5

Revision:      $Id$
'''

from   __future__ import division

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


__all__ = ['Genotype','UnphasedGenotypeModel','GenotypeArray','GenotypeError']


class GenotypeError(ValueError): pass


try:
  from   _genoarray import GenotypeArray, Genotype, GenotypeArrayDescriptor,UnphasedMarkerModel

except ImportError:
  from   array     import array
  from   math      import log, ceil

  from   bitarray  import getbits,setbits

  MISSING,HEMIZYGOTE,HETEROZYGOTE,HOMOZYGOTE=range(4)

  class Genotype(object):
    __slots__ = ('model','allele1','allele2','index','gclass')

    def __init__(self, model, allele1, allele2, index):
      self.model   = model
      self.allele1 = allele1
      self.allele2 = allele2
      self.index   = index

      missing1 = allele1 is None
      missing2 = allele2 is None

      if missing1 and missing2:
        self.gclass = MISSING
      elif missing1 or missing2:
        self.gclass = HEMIZYGOTE
      elif allele1 is allele2:
        self.gclass = HOMOZYGOTE
      else:
        if allele1 == allele2:
          raise GenotypeError('Attempt to add non-singleton alleles')
        self.gclass = HETEROZYGOTE

    def alleles(self):
      return (self.allele1,self.allele2)

    def heterozygote(self):
      return self.gclass == HETEROZYGOTE

    def homozygote(self):
      return self.gclass == HOMOZYGOTE

    def hemizygote(self):
      return self.gclass == HEMIZYGOTE

    def missing(self):
      return self.gclass == MISSING

    def __nonzero__(self):
      return self.gclass != MISSING

    def __getitem__(self,i):
      return self.alleles()[i]

    def __len__(self):
      return len([a for a in self.alleles() if a is not None])

    def __repr__(self):
      return '<Genotype: %s/%s at 0x%X>' % (self.allele1,self.allele2,id(self))


  def genotype_bit_width(n,allow_hemizygote):
    if allow_hemizygote:
      m = (n+1)*(n+2)//2
    else:
      m = n*(n+1)//2 + 1

    return int(ceil(log(m)/log(2.0)))


  def byte_array_size(nbits):
    return int(ceil(nbits/8))


  class GenotypeArrayDescriptor(object):
    __slots__ = ('models','offsets','byte_size','bit_size')

    def __init__(self, models, initial_offset=0):
      n = len(models)
      offsets = [0]*(n+1)

      offsets[0] = initial_offset
      for i,m in enumerate(models):
        offsets[i+1] = offsets[i] + m.bit_width

      self.models    = models
      self.offsets   = offsets
      self.bit_size  = offsets[-1]
      self.byte_size = byte_array_size(self.bit_size)

    def __len__(self):
      return len(self.models)


  class GenotypeArray(object):
    __slots__ = ('descriptor','data')

    def __init__(self, descriptor, genos=None):
      if isinstance(descriptor, GenotypeArrayDescriptor):
        self.descriptor = descriptor
      elif isinstance(descriptor, GenotypeArray):
        self.descriptor = descriptor.descriptor

      self.data = array('B', [0]*self.descriptor.byte_size)

      if genos is not None:
        self[:] = genos

    def __len__(self):
      return len(self.descriptor.models)

    def __getitem__(self, i):
      descr = self.descriptor

      if isinstance(i,slice):
        x = xrange(*i.indices(len(descr.models)))
        return [ self[i] for i in x ]

      model    = descr.models[i]
      startbit = descr.offsets[i]
      width    = model.bit_width
      j        = getbits(self.data, startbit, width)

      return model.genotypes[j]

    def __setitem__(self, i, geno):
      descr = self.descriptor

      if isinstance(i,slice):
        x = xrange(*i.indices(len(descr.models)))
        try:
          n = len(geno)
        except TypeError:
          geno = list(geno)
          n = len(geno)
        if len(x) != n:
          raise IndexError('Invalid slice')
        for i,j in enumerate(x):
          self[j] = geno[i]
        return
      elif isinstance(geno,tuple) and len(geno) == 2:
        geno = descr.models[i][geno]
      elif not isinstance(geno,Genotype):
        raise GenotypeError('Invalid genotype: %s' % geno)

      model    = descr.models[i]
      startbit = descr.offsets[i]
      width    = model.bit_width

      assert geno.model is model
      setbits(self.data, startbit, geno.index, width)


  class UnphasedMarkerModel(object):
    __slots__ = ('alleles','genotypes','genomap','bit_width', 'allow_hemizygote')

    def __init__(self, allow_hemizygote=False, max_alleles=2):
      self.genomap          = {}
      self.genotypes        = []
      self.alleles          = []
      self.max_alleles      = max_alleles
      self.bit_width        = genotype_bit_width(max_alleles,allow_hemizygote)
      self.allow_hemizygote = allow_hemizygote
      self.add_genotype( (None,None) )

    def get_allele(self, allele):
      return self.alleles.index(allele)

    def add_allele(self, allele):
      if allele in self.alleles:
        return self.alleles.index(allele)

      n = len(self.alleles)
      new_width = genotype_bit_width(n,self.allow_hemizygote)
      if new_width > self.bit_width:
        raise GenotypeError('Allele cannot be added to model due to fixed bit width')
      self.alleles.append(a)

      return n

    def get_genotype(self, geno):
      return self.genomap[geno]

    __getitem__ = get_genotype

    def add_genotype(self, geno):
      g = self.genomap.get_genotype(geno)

      # If the genotype has not already been seen for this locus
      if g is not None:
        return g

      allele1,allele2 = sorted(geno)

      index1 = self.add_allele(allele1)
      index2 = self.add_allele(allele2)

      allele1 = self.alleles[index1]
      allele2 = self.alleles[index2]

      # Create and save new genotype
      g = Genotype(self, index1, index2, len(self.genotypes))

      if not self.allow_hemizygote and g.hemizygote():
        raise GenotypeError('Genotype model does not all hemizygous genotypes')

      self.genotypes.append(g)
      self.genomap[allele1,allele2] = g
      self.genomap[allele2,allele1] = g

      return g


def model_from_alleles(alleles, allow_hemizygote=False, max_alleles=None):
  alleles = sorted(set(a for a in alleles if a is not None))
  n = len(alleles)

  if max_alleles is None:
    max_alleles = n

  if allow_hemizygote:
    alleles = [None]+alleles

  n = len(alleles)
  genos = [ (alleles[i],alleles[j]) for i in range(n) for j in range(i,n) ]

  model = UnphasedMarkerModel(allow_hemizygote=allow_hemizygote,max_alleles=max_alleles)
  for g in genos:
    model.add_genotype(g)

  return model

def model_from_genotypes(genotypes, allow_hemizygote=None, max_alleles=None):
  alleles = sorted(set(a for g in genoset for a in g if a is not None))
  return model_from_alleles_and_genotypes(alleles, genotypes, allow_hemizygote, max_alleles)

def model_from_alleles_and_genotypes(alleles, genotypes, allow_hemizygote=False, max_alleles=None):
  genoset = set(genotypes)

  if not allow_hemizygote:
    def hemi(g):
      return (g[0] is None) ^ (g[1] is None)

    hemi = any(hemi(g) for g in genoset)

    if allow_hemizygote is not None and hemi:
      raise GenotypeError('Genotype model does not allow hemizygous genotypes')

    allow_hemizygote = hemi

  n = len(set(alleles))

  if max_alleles is None:
    max_alleles = n
  elif n > max_alleles:
    raise GenotypeError('Genotype model supports at most %d alleles, %d specified' % (max_alleles,n))

  model = UnphasedMarkerModel(allow_hemizygote=allow_hemizygote,max_alleles=max_alleles)

  for a in alleles:
    model.add_allele(a)

  for g in genotypes:
    model.add_genotype(g)

  if allow_hemizygote:
    alleles = [None]+alleles

  all_genos = ( (alleles[i],alleles[j]) for i in range(n) for j in range(i,n) )

  for g in all_genos:
    model.add_genotype(g)

  return model

def model_from_complete_alleles_and_genotypes(alleles, genotypes, allow_hemizygote=False, max_alleles=None):
  if max_alleles is None:
    max_alleles = len(set(alleles))

  model = UnphasedMarkerModel(allow_hemizygote=allow_hemizygote,max_alleles=max_alleles)

  for a in alleles:
    model.add_allele(a)

  for g in genotypes:
    model.add_genotype(g)

  return model


def main():
  import time
  import genoarray

  from   itertools import izip,repeat,imap
  from   operator  import getitem

  def parse_geno(g):
    g = g.strip()
    if not g:
      return None,None
    if len(g) == 1:
      return None,g
    else:
      return tuple(g)

  n = 500
  m = 5000

  genos = ['AA','AC','CC','AA','CA','',' A','C ']*m

  t1 = time.clock()

  if 1:
    model   = model_from_alleles('ACGT',allow_hemizygote=True)
    descr   = GenotypeArrayDescriptor( [model]*len(genos) )
    genomap = dict( (g,m.get_genotype(parse_geno(g))) for m,g in izip(descr.models,genos) )
    print descr.bit_size,descr.byte_size, float(descr.bit_size)/len(descr)

  if 1:
    for i in range(n):
      e = GenotypeArray(descr, imap(getitem, repeat(genomap), genos))
      f = GenotypeArray(e, e)
      d = map(Genotype.alleles, e)

  t2 = time.clock()

  if 1:
    for i in range(n):
      e = genoarray.snp_acgt2.pack_strs(genos)
      f = genoarray.snp_acgt2.pack_reps(e)
      d = genoarray.snp_acgt2.genos_from_reps(e)

  t3 = time.clock()

  if 1:
    for i in range(n):
      e = genoarray.snp_acgt.pack_strs(genos)
      f = genoarray.snp_acgt.pack_reps(e)
      d = genoarray.snp_acgt.genos_from_reps(e)

  t4 = time.clock()

  if 1:
    m = n//10
    s=float(n)/m
    for i in range(m):
      e = genoarray.snp_marker.pack_strs(genos)
      f = genoarray.snp_marker.pack_reps(e)
      d = genoarray.snp_marker.genos_from_reps(e)

  t5 = time.clock()

  if 1:
    print '_genoarray:',t2-t1
    print 'snp_acgt2 :',t3-t2
    print 'snp_acgt  :',t4-t3
    print 'snp_marker:',(t5-t4)*s
    print
    print '_genoarray/snp_acgt  :',(t2-t1)/(t4-t3),(t4-t3)/(t2-t1)
    print 'snp_acgt2 /snp_acgt  :',(t3-t2)/(t4-t3),(t4-t3)/(t3-t2)
    print
    print 'snp_acgt2 /snp_marker:',(t3-t2)/(t5-t4)/s,(t5-t4)/(t3-t2)*s
    print '_genoarray/snp_marker:',(t2-t1)/(t5-t4)/s,(t5-t4)/(t2-t1)*s

  if 0:
    snp_ab = model_from_alleles('AB',allow_hemizygote=True)
    print snp_ab.bit_width
    print snp_ab.alleles
    print snp_ab.genotypes
    print snp_ab.genomap

  if 0:
    print len(snp_ab.get_genotype( (None,None) ) )
    print len(snp_ab.get_genotype( ('A',None)  ) )
    print len(snp_ab.get_genotype( (None,'A')  ) )
    print len(snp_ab.get_genotype( ('A', 'A')  ) )
    print bool(snp_ab.get_genotype( (None,None) ) )
    print bool(snp_ab.get_genotype( ('A',None)  ) )
    print bool(snp_ab.get_genotype( (None,'A')  ) )
    print bool(snp_ab.get_genotype( ('A', 'A')  ) )
    print 'A'  in snp_ab.get_genotype( (None,None) )
    print 'A'  in snp_ab.get_genotype( ('A',None)  )
    print 'A'  in snp_ab.get_genotype( (None,'A')  )
    print 'A'  in snp_ab.get_genotype( ('A', 'A')  )
    print 'B'  in snp_ab.get_genotype( (None,None) )
    print 'B'  in snp_ab.get_genotype( ('A',None)  )
    print 'B'  in snp_ab.get_genotype( (None,'A')  )
    print 'B'  in snp_ab.get_genotype( ('A', 'A')  )
    print None in snp_ab.get_genotype( (None,None) )
    print None in snp_ab.get_genotype( ('A',None)  )
    print None in snp_ab.get_genotype( (None,'A')  )
    print None in snp_ab.get_genotype( ('A', 'A')  )
    print snp_ab.get_genotype( (None,None) )[0],snp_ab.get_genotype( (None,None) )[1]
    print snp_ab.get_genotype( ('A',None)  )[0],snp_ab.get_genotype( ('A',None)  )[1]
    print snp_ab.get_genotype( (None,'A')  )[0],snp_ab.get_genotype( (None,'A')  )[1]
    print snp_ab.get_genotype( ('A', 'A')  )[0],snp_ab.get_genotype( ('A', 'A')  )[1]

  if 0:
    genos = ['AA','AA','','BB','BA','AA']
    print snp_ab.genotypes
    descr = GenotypeArrayDescriptor( [snp_ab]*len(genos) )
    e = GenotypeArray(descr, (parse_geno(g) for g in genos))
    print descr.byte_size
    print e
    print e.data
    print buffer(e.data)


if __name__ == '__main__':
  if 0:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
  else:
    main()
