# -*- coding: utf-8 -*-

from   __future__ import division

__abstract__  = 'Efficient bit-packed genotype array representation'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


__all__ = ['Genotype','UnphasedMarkerModel','GenotypeArray','GenotypeArrayDescriptor',
           'GenotypeLookupError','GenotypeRepresentationError']


from   itertools     import izip,chain

import numpy as np

from   glu.lib.utils import izip_exact


MISSING,HEMIZYGOTE,HOMOZYGOTE,HETEROZYGOTE=range(4)

# Set to test pure-Python implementation
FORCE_PYTHON=False


try:
  if FORCE_PYTHON:
    raise ImportError

  from   glu.lib.genolib._genoarray import (GenotypeArray, Genotype, GenotypeArrayDescriptor,
                                            UnphasedMarkerModel, genotype_indices,
                                            count_genotypes, genotype_categories,
                                            locus_summary, sample_summary, genoarray_concordance,
                                            genoarray_ibs,
                                            GenotypeLookupError, GenotypeRepresentationError,
                                            pick, pick_columns, place, place_list)

  GENO_ARRAY_VERSION='C'

except ImportError:
  GENO_ARRAY_VERSION='Python'

  import sys
  from   math      import log, ceil

  from   glu.lib.genolib.bitarray  import getbits,setbits

  sys.stderr.write('[WARNING] Using slow Python genoarray\n')

  def _hemi(geno):
    return (geno[0] is None) ^ (geno[1] is None)


  class GenotypeLookupError(KeyError): pass
  class GenotypeRepresentationError(ValueError): pass


  class Genotype(object):
    '''
    Python implementation of a genotype representation object.  Genotype
    objects belong to a specific model, which is responsible for allocation
    genotypes, managing alleles, and assigning a unique index to each
    possible genotype within a model.
    '''

    __slots__ = ('model','allele1_index','allele2_index','index','category')

    def __init__(self, model, allele1, allele2, index):
      '''
      The Genotype constructor should not be called by users, only by model
      objects.

      @param    model: genotype representation
      @type     model: UnphasedMarkerRepresentation or similar object
      @param  allele1: the first allele
      @type   allele1: str
      @param  allele2: the second allele
      @type   allele2: str
      @param    index: index of genotype within the model
      @type     index: int
      '''
      self.model   = model
      self.index   = index

      self.allele1_index = model.alleles.index(allele1)
      self.allele2_index = model.alleles.index(allele2)

      missing1 = allele1 is None
      missing2 = allele2 is None

      if missing1 and missing2:
        self.category = MISSING
      elif missing1 or missing2:
        self.category = HEMIZYGOTE
      elif allele1 is allele2:
        self.category = HOMOZYGOTE
      else:
        if allele1 == allele2:
          raise GenotypeRepresentationError('Attempt to add non-singleton alleles')
        self.category = HETEROZYGOTE

    @property
    def allele1(self):
      return self.model.alleles[self.allele1_index]

    @property
    def allele2(self):
      return self.model.alleles[self.allele2_index]

    def alleles(self):
      '''
      Return a tuple of alleles
      '''
      alleles = self.model.alleles
      allele1 = alleles[self.allele1_index]
      allele2 = alleles[self.allele2_index]
      return (allele1,allele2)

    def heterozygote(self):
      '''
      Return whether this genotype is heterozygous
      '''
      return self.category == HETEROZYGOTE

    def homozygote(self):
      '''
      Return whether this genotype is homozygous
      '''
      return self.category == HOMOZYGOTE

    def hemizygote(self):
      '''
      Return whether this genotype is hemizygous
      '''
      return self.category == HEMIZYGOTE

    def missing(self):
      '''
      Return whether this genotype is the missing genotype
      '''
      return self.category == MISSING

    def __nonzero__(self):
      '''
      Return whether this genotype is not the missing genotype
      '''
      return self.category != MISSING

    def __getitem__(self,i):
      '''
      Return the i'th allele, for i in [0,1]
      '''
      return self.alleles()[i]

    def __len__(self):
      '''
      Return the number of non-missing alleles
      '''
      return len([a for a in self.alleles() if a is not None])

    def __repr__(self):
      '''
      Return a string representation of the alleles
      '''
      return repr(self.alleles())

    def __hash__(self):
      return hash(self.alleles())

    def __eq__(self,other):
      '''
      Compare this genotype with another genotype or tuple for equality.

      If other is a genotype, then the comparison is by object identity.
      If other is a tuple, then the comparison is by allele tuple equivalence.

      @param other: other genotype
      @type  other: Genotype or tuple
      '''
      geno_self  = isinstance(self,Genotype)
      geno_other = isinstance(other,Genotype)

      if geno_self:
        self = self.alleles()
      if geno_other:
        other = other.alleles()

      if not isinstance(self, tuple) or len(self) !=2 or \
         not isinstance(other,tuple) or len(other)!=2:
        return NotImplemented

      return sorted(self)==sorted(other)

    def __ne__(self,other):
      '''
      Compare this genotype with another genotype or tuple for inequality.

      If other is a genotype, then the comparison is by object identity.
      If other is a tuple, then the comparison is by allele tuple equivalence.

      @param other: other genotype
      @type  other: Genotype or tuple
      '''
      geno_self  = isinstance(self,Genotype)
      geno_other = isinstance(other,Genotype)

      if geno_self:
        self = self.alleles()
      if geno_other:
        other = other.alleles()

      if not isinstance(self, tuple) or len(self) !=2 or \
         not isinstance(other,tuple) or len(other)!=2:
        return NotImplemented

      return sorted(self)!=sorted(other)

    def __lt__(self,other):
      '''
      Test if this genotype is less than with another genotype or tuple.
      Comparison is performed on allele tuples.

      @param other: other genotype
      @type  other: Genotype or tuple
      '''
      if isinstance(self,Genotype):
        self = self.alleles()
      if isinstance(other,Genotype):
        other = other.alleles()

      if not isinstance(self, tuple) or len(self) !=2 or \
         not isinstance(other,tuple) or len(other)!=2:
        return NotImplemented

      return sorted(self)<sorted(other)

    def __le__(self,other):
      '''
      Test if this genotype is less than or equal to another genotype or
      tuple.  Comparison is performed on allele tuples.

      @param other: other genotype
      @type  other: Genotype or tuple
      '''

      if isinstance(self,Genotype):
        self = self.alleles()
      if isinstance(other,Genotype):
        other = other.alleles()

      if not isinstance(self, tuple) or len(self) !=2 or \
         not isinstance(other,tuple) or len(other)!=2:
        return NotImplemented

      return sorted(self)<=sorted(other)

  def genotype_bit_size(n,allow_hemizygote):
    '''
    Return the number of bits of storage required to represent a genotypes
    with n alleles, possibly supporting hemizygote genotypes.

    @param                 n: number of alleles
    @type                  n: int
    @param  allow_hemizygote: flag indicating if hemizygote genotypes are to be supported
    @type   allow_hemizygote: bool
    @return                 : number of bits required
    @rtype                  : int
    '''

    if allow_hemizygote:
      m = (n+1)*(n+2)//2
    else:
      m = n*(n+1)//2 + 1

    return int(ceil(log(m)/log(2.0)))

  def byte_array_size(n):
    '''
    Return the number of whole bytes required to represent n bits.

    @param  nbits: bit size
    @type   nbits: int
    @return      : byte size
    @rtype       : int
    '''
    return int(ceil(n/8))


  class GenotypeArrayDescriptor(object):
    __slots__ = ('_models','offsets','byte_size','bit_size','max_bit_size','homogeneous')

    def __init__(self, models, initial_offset=0):
      '''
      Construct a new GenotypeArrayDescriptor
      '''
      n = len(models)

      homogeneous = models[0].bit_size if models else 0

      offsets = [0]*(n+1)

      offsets[0] = initial_offset
      max_bit_size = 0
      for i,m in enumerate(models):
        offsets[i+1] = offsets[i] + m.bit_size
        max_bit_size = max(max_bit_size,m.bit_size)
        if m.bit_size != homogeneous:
          homogeneous = 0

      self._models      = models
      self.offsets      = offsets
      self.bit_size     = offsets[-1]
      self.byte_size    = byte_array_size(self.bit_size)
      self.max_bit_size = max_bit_size
      self.homogeneous  = homogeneous

    def __len__(self):
      return len(self._models)

    def __iter__(self):
      return iter(self._models)

    def __getitem__(self, i):
      return self._models[i]

    def __setitem__(self, i, new_model):
      '''
      Replace a genotype model with another provided that:

        1. bit size is identical between new and old models
        2. new model defines an equal or greater number of alleles as the old
        3. new model defines an equal or greater number of genotypes as the old

      Replacing a model will change the decoding of all genotypes and
      genotype arrays that refer to this descriptor

      @param         i: index of model to replace
      @type          i: index type
      @param new_model: new model
      @type  new_model: UnphasedMarkerModel
      '''
      if isinstance(i,slice):
        x = xrange(*i.indices(len(self)))
        try:
          n    = len(new_model)
        except TypeError:
          new_model = list(new_model)
          n         = len(new_model)

        if len(x) != n:
          raise IndexError('Invalid slice')

        for j,n in izip(x,new_model):
          self[j] = n

        return

      old_model = self._models[i]

      if old_model.bit_size != new_model.bit_size:
        raise GenotypeRepresentationError('Model may not be replaced due to incompatible size')

      if len(old_model.alleles) > len(new_model.alleles):
        raise GenotypeRepresentationError('Model may not be replaced since new model defines too few alleles')

      if len(old_model.genotypes) > len(new_model.genotypes):
        raise GenotypeRepresentationError('Model may not be replaced since new model defines too few genotypes')

      n = len(old_model.genotypes)
      if old_model.genotypes != new_model.genotypes[:n]:
        raise GenotypeRepresentationError('Model may not be replaced since new model defines incompatible genotypes')

      self._models[i] = new_model


  class GenotypeArray(object):
    __slots__ = ('descriptor','data')

    def __init__(self, descriptor, genos=None):
      '''
      Construct a new GenotypeArray with given GenotypeArrayDescriptor and
      optionally initialized with a sequence of genotypes.

      @param   descriptor: genotype array representation
      @type    descriptor: GenotypeArrayDescriptor or GenotypeArray
      @param        genos: initial genotypes
      @type         genos: sequence of genotypes objects or tuples
      '''
      if isinstance(descriptor, GenotypeArrayDescriptor):
        self.descriptor = descriptor
      elif isinstance(descriptor, GenotypeArray):
        self.descriptor = descriptor.descriptor

      self.data = np.zeros(self.descriptor.byte_size, dtype=np.ubyte)

      if genos is not None:
        self[:] = genos

    def __len__(self):
      '''
      Return the length of this array

      @return: array length
      @rtype : int
      '''
      return len(self.descriptor)

    def __getitem__(self, i):
      '''
      Return the i'th genotype or slice or list of genotypes by index

      @param i: index or slice or list of indexes
      @type  i: int or slice or list of int
      @return : i'th genotype if an index is specified, sequence of
                genotypes if a slice or list is specified
      @rtype  : Genotype object or sequence of Genotype objects
      '''
      descr = self.descriptor

      if isinstance(i,slice):
        x = xrange(*i.indices(len(descr)))
        return [ self[i] for i in x ]
      elif isinstance(i,(list,tuple)):
        return [ self[j] for j in i ]

      model    = descr[i]
      startbit = descr.offsets[i]
      width    = model.bit_size
      j        = getbits(self.data, startbit, width)

      return model.genotypes[j]

    def __setitem__(self, i, geno):
      '''
      Set the i'th genotype to geno or slice

      @param    i: index or slice
      @type     i: int or slice
      @param geno: genotype if an index is specified, sequence of genotypes
                   if a slice is specified
      @type  geno: slice, tuple, Genotype object
      '''
      descr = self.descriptor

      if isinstance(i,slice):
        x = xrange(*i.indices(len(descr)))
        try:
          n    = len(geno)
        except TypeError:
          geno = list(geno)
          n    = len(geno)

        if len(x) != n:
          raise IndexError('Invalid slice')

        for g,j in izip(geno,x):
          self[j] = g

        return

      elif isinstance(geno,tuple) and len(geno) == 2:
        geno = descr[i][geno]
      elif not isinstance(geno,Genotype):
        raise GenotypeRepresentationError('Invalid genotype: %s' % geno)

      model    = descr[i]
      startbit = descr.offsets[i]
      width    = model.bit_size

      if isinstance(geno,Genotype):
        if geno.model is not model:
          geno = model[geno.alleles()]
      else:
        geno = model[geno]

      setbits(self.data, startbit, geno.index, width)

    def __repr__(self):
      return repr(list(self))

    def indices(self):
      '''
      Return an array of integer genotype indices
      '''
      data    = self.data
      offsets = self.descriptor.offsets
      inds    = [ getbits(data, offsets[i], offsets[i+1]-offsets[i]) for i in xrange(len(self)) ]
      return np.asarray(inds,dtype=int)

    def counts(self,counts=None,model=None):
      '''
      Return an array of genotype counts by index
      '''
      return count_genotypes(self,counts,model)

    def categories(self,counts=None):
      '''
      Count the number of occurrences of each genotypes category
      '''
      return genotype_categories(self,counts)

    def tolist(self):
      '''
      Return a list of genotypes
      '''
      return self[:]


  class UnphasedMarkerModel(object):
    '''
    bit-packed unphased marker representation where the genotype representation
    and internal representation are the same.
    '''

    __slots__ = ('alleles','genotypes','genomap','bit_size','allow_hemizygote','max_alleles')

    def __init__(self, allow_hemizygote=False, max_alleles=None):
      '''
      Construct a new UnphasedMarkerModel

      This class represents bidirectional mappings of genotypes between
      strings and Python objects.  The object representation of a genotype
      is a list of two alleles or up to the max_alleles. Given this
      representation, the alleles names need not be known in advance, only
      the maximum number that may be stored within the model.

      @param  allow_hemizygote: flag indicating if hemizygote is allowed
      @type   allow_hemizygote: bool
      @param       max_alleles: maximum number of alleles allowed. Default is 2
      @type        max_alleles: int or None
      '''

      self.genomap          = {}
      self.genotypes        = []
      self.alleles          = []
      self.max_alleles      = max(2,max_alleles)
      self.bit_size         = genotype_bit_size(self.max_alleles,allow_hemizygote)
      self.allow_hemizygote = allow_hemizygote
      self.add_genotype( (None,None) )

    def get_allele(self, allele):
      '''
      Return the representation for a given allele.  If the specified allele
      does not belong to this model, a GenotypeLookupError is raised.

      @param allele: allele object
      @type  allele: object, usually str
      @return      : allele representation
      @rtype       : int
      '''
      try:
        return self.alleles.index(allele)
      except ValueError:
        raise GenotypeLookupError,allele

    def add_allele(self, allele):
      '''
      Return the representation for a given allele.  If the specified
      allele does not belong to this model it will be added. A GenotypeRepresentationError
      is raised if there is not sufficient space within the model.

      @param allele: allele object
      @type  allele: object, usually str
      @return      : allele representation
      @rtype       : int
      '''
      if allele in self.alleles:
        return self.alleles.index(allele)

      n = len(self.alleles)
      new_width = genotype_bit_size(n,self.allow_hemizygote)
      if new_width > self.bit_size:
        raise GenotypeRepresentationError('Allele cannot be added to model due to fixed bit width')
      self.alleles.append(allele)

      return n

    def get_genotype(self, geno):
      '''
      Return the representation for a given genotype tuple.  If the
      specified genotype tuple does not belong to this model, a GenotypeLookupError is
      raised.

      @param geno: genotype tuple
      @type  geno: 2-tuple of allele objects
      @return    : genotype object
      @rtype     : Genotype
      '''
      if isinstance(geno,Genotype):
        if geno.model is self:
          return geno
        geno = geno.alleles()

      try:
        return self.genomap[geno]
      except KeyError:
        if (geno[0] in self.alleles and geno[1] in self.alleles and
           (self.allow_hemizygote or not _hemi(geno))):
          return self.add_genotype(geno)
        raise GenotypeLookupError,geno

    __getitem__ = get_genotype

    def __contains__(self, geno):
      '''
      Return if a given genotype tuple or object is contained within this model.

      @param geno: genotype tuple
      @type  geno: 2-tuple of allele objects
      @return    : True if geno is contained within this model, otherwise False
      @rtype     : bool
      '''
      try:
        self.get_genotype(geno)
        return True
      except GenotypeLookupError:
        return False

    def add_genotype(self, geno):
      '''
      Return the representation for a given genotype tuple.  If the
      specified genotype tuple does not belong to this model, it will be added.

      @param geno: genotype tuple
      @type  geno: 2-tuple of allele objects
      @return    : genotype object
      @rtype     : Genotype
      '''
      if isinstance(geno,Genotype):
        if geno.model is self:
          return geno
        geno = geno.alleles()

      g = self.genomap.get(geno)

      # If the genotype has not already been seen for this locus
      if g is not None:
        return g

      allele1,allele2 = sorted(geno)

      index1 = self.add_allele(allele1)
      index2 = self.add_allele(allele2)

      allele1 = self.alleles[index1]
      allele2 = self.alleles[index2]

      # Create and save new genotype
      g = Genotype(self, allele1, allele2, len(self.genotypes))

      if not self.allow_hemizygote and g.hemizygote():
        raise GenotypeRepresentationError('Genotype model does not allow hemizygous genotypes')

      self.genotypes.append(g)
      self.genomap[allele1,allele2] = g
      self.genomap[allele2,allele1] = g

      return g

    def replaceable_by(self, other):
      '''
      Return model could replace the other without altering existing
      encodings.  To be compatible this model must:

        1. Have identical bit width, maximum alleles, and hemizygosity
        2. the other model must have a superset of the alleles
        3. this model's genotypes must be a prefix of genotypes of other

      >>> model1 = build_model('A')
      >>> model2 = build_model('AB')
      >>> model3 = build_model('AB',allow_hemizygote=True)
      >>> model4 = build_model('AB',genotypes=[('B','B')])

      >>> model1.replaceable_by(model1)
      True
      >>> model2.replaceable_by(model2)
      True
      >>> model1.replaceable_by(model2)
      True
      >>> model2.replaceable_by(model1)
      False
      >>> model2.replaceable_by(model3)
      False
      >>> model1.replaceable_by(model4)
      False
      >>> model2.replaceable_by(model4)
      False
      '''
      if self is other:
        return True

      if (self.bit_size != other.bit_size or self.max_alleles      != other.max_alleles
                                          or self.allow_hemizygote != other.allow_hemizygote):
        return False

      if not (set(other.alleles)>=set(self.alleles)):
        return False

      if len(self.genotypes) > len(other.genotypes):
        return False

      if any( g1!=g2 for g1,g2 in izip(self.genotypes,other.genotypes) ):
        return False

      return True


  def genoarray_concordance(genos1, genos2):
    '''
    Generate simple concordance statistics from two genotype arrays

    @param genos1: first input
    @type  genos1: sequence of genotypes
    @param genos2: second input
    @type  genos2: sequence of genotypes
    @return      : number of concordant genotypes and number of non-missing comparisons made
    @rtype       : (int,int)
    '''
    if len(genos1) != len(genos2):
      raise ValueError("genotype vector sizes do not match: %zd != %zd" % (len(genos1),len(genos2)))

    concordant = comparisons = 0
    for a,b in izip(genos1,genos2):
      if a and b:
        if a==b:
          concordant += 1
        comparisons += 1

    return concordant,comparisons


  def genoarray_ibs(genos1, genos2):
    '''
    Generate counts of alleles shared IBS between two genotype arrays

    @param genos1: first input
    @type  genos1: sequence of genotypes
    @param genos2: second input
    @type  genos2: sequence of genotypes
    @return      : Tuple (IBS0, IBS1, IBS2) where IBSn is the count of genotypes that share n alleles IBS
    @rtype       : (int,int,int)
    '''
    if len(genos1) != len(genos2):
      raise ValueError("genotype vector sizes do not match: %zd != %zd" % (len(genos1),len(genos2)))

    ibs0,ibs1,ibs2 = 0,0,0
    for a,b in izip(genos1,genos2):
      if a and b:
        if a==b:
          ibs2 += 1
        elif (a[0]==b[0] or a[0]==b[1] or
              a[1]==b[0] or a[1]==b[1]):
          ibs1 += 1
        else:
          ibs0 += 1

    return ibs0,ibs1,ibs2


  def genotype_indices(genos):
    '''
    Return an array of integer genotype indices
    '''
    if isinstance(genos, GenotypeArray):
      return genos.indices()
    return [ g.index for g in genos ]


  def count_genotypes(genos,counts=None,model=None):
    '''
    Count the number of occurrences of each genotypes belonging to a single
    model given a sequence of genotypes.  Counts are returned as a list of
    integers corresponding to the number of each genotype in model.genotype
    observed.

    @param model: model for all genotypes
    @type  model: UnphasedMarkerModel
    @param genos: sequence of genotype objects belonging to model
    @type  genos: sequence of Genotype instances
    @return     : count of each genotype
    @rtype      : list of integers

    >>> import random
    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*1400)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes*200+model.genotypes[1:]*200
    >>> random.shuffle(genos)

    >>> count_genotypes(genos)
    array([200, 400, 400, 400])
    >>> count_genotypes(GenotypeArray(descr,genos))
    array([200, 400, 400, 400])

    >>> c=count_genotypes(genos)
    >>> count_genotypes(genos,counts=c)
    array([400, 800, 800, 800])
    >>> count_genotypes(genos,counts=c)
    array([ 600, 1200, 1200, 1200])
    '''
    if not genos and not model:
      raise GenotypeRepresentationError('empty genotype sequence with unknown model')

    if not model:
      model = genos[0].model

    if counts is None:
      counts = [0]*len(model.genotypes)
    else:
      if len(counts) != len(model.genotypes):
        raise ValueError('invalid count array')

    for geno in genos:
      if geno.model is not model:
        raise GenotypeRepresentationError('genotype models do not match')
      counts[geno.index] += 1

    return np.asarray(counts,dtype=int)


  def genotype_categories(genos,counts=None):
    '''
    Count the number of occurrences of each genotypes category given a
    sequence of genotypes.  Counts are returned as a list of integers
    corresponding to the number of missing, hemizygote, homozygote, and
    heterozygote genotype observed.

    @param genos: sequence of genotype objects belonging to model
    @type  genos: sequence of Genotype instances
    @return     : missing, hemizygote, heterozygote, and homozygote counts
    @rtype      : list of integers

    >>> import random
    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*1400)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes*200+model.genotypes[1:]*200
    >>> random.shuffle(genos)

    >>> genotype_categories(genos)
    array([200,   0, 800, 400])
    >>> genotype_categories(GenotypeArray(descr,genos))
    array([200,   0, 800, 400])
    '''
    if counts is None:
      counts = [0]*4
    elif len(counts) != 4:
      raise ValueError('invalid count array')

    for geno in genos:
      counts[geno.category] += 1
    return np.array(counts,dtype=int)


  def locus_summary(genos,sample_counts=None,locus_counts=None,model=None):
    '''
    Count the number of occurrences of each n genotypes for a given locus and
    the category for each sample given a sequence of genotypes at a locus.
    Genotype counts are returned as a list of integers counts corresponding
    to genotype index.  Sample category counts are returned as an n x 4 array
    with the number of missing, hemizygote, homozygote, and heterozygote
    genotypes observed for each sample.

    @param         genos: genotypes belonging to a single model
    @type          genos: sequence
    @param sample_counts: optional array of category counts to accumulate
    @type  sample_counts: same type as np.zeros( (len(genos),4), dtype=int)
    @param  locus_counts: optional array of genotype counts to accumulate
    @type   locus_counts: same type as np.zeros(len(genos), dtype=int)
    @return             : 2-tuple of locus_counts and sample_counts
    @rtype              : tuple of np.array

    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*4)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes

    >>> l,s = locus_summary(genos)
    >>> l
    array([1, 1, 1, 1])
    >>> s
    array([[1, 0, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1],
           [0, 0, 1, 0]])
    >>> l,s=locus_summary(GenotypeArray(descr,reversed(genos)),s,l)
    >>> l
    array([2, 2, 2, 2])
    >>> s
    array([[1, 0, 1, 0],
           [0, 0, 1, 1],
           [0, 0, 1, 1],
           [1, 0, 1, 0]])
    '''
    if not genos and not model:
      raise GenotypeRepresentationError('empty genotype sequence with unknown model')

    if not model:
      model = genos[0].model

    n = len(genos)

    if locus_counts is None:
      locus_counts = [0]*(1<<model.bit_size)
    elif len(locus_counts) < len(model.genotypes):
      raise ValueError('invalid locus count array')

    if sample_counts is None:
      sample_counts = np.zeros((n,4), dtype=int)
    elif sample_counts.shape != (n,4):
      raise ValueError('invalid sample count array')

    for i,geno in enumerate(genos):
      if geno.model is not model:
        raise GenotypeRepresentationError('genotype models do not match')
      locus_counts[geno.index]       += 1
      sample_counts[i,geno.category] += 1

    locus_counts = np.asarray(locus_counts,dtype=int)

    return locus_counts,sample_counts


  def sample_summary(genos,locus_counts=None,sample_counts=None):
    '''
    Count the number of occurrences of each genotype category and the of
    occurrences of each of n genotypes for a given sample.  Genotype category
    counts are returned as a list of integers counts corresponding to
    genotype index.  Per-locus genotype counts are returned as an n x
    max(genos) array with the number of missing, hemizygote, homozygote, and
    heterozygote genotypes observed for each sample.

    @param         genos: genotypes belonging to a single model
    @type          genos: sequence
    @param  locus_counts: optional array of genotype counts to accumulate
    @type   locus_counts: same type as np.zeros((len(genos),len(maxgenos)), dtype=int)
    @param sample_counts: optional array of category counts to accumulate
    @type  sample_counts: same type as np.zeros(4, dtype=int)
    @return             : 2-tuple of  sample_counts and locus_counts
    @rtype              : tuple of np.array

    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*4)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes

    >>> s,l = sample_summary(genos)
    >>> s
    array([1, 0, 2, 1])
    >>> l
    array([[1, 0, 0, 0],
           [0, 1, 0, 0],
           [0, 0, 1, 0],
           [0, 0, 0, 1]])
    >>> s,l=sample_summary(GenotypeArray(descr,reversed(genos)),l,s)
    >>> s
    array([2, 0, 4, 2])
    >>> l
    array([[1, 0, 0, 1],
           [0, 1, 1, 0],
           [0, 1, 1, 0],
           [1, 0, 0, 1]])
    '''
    n = len(genos)

    if sample_counts is None:
      sample_counts = [0]*4
    elif len(sample_counts) != 4:
      raise ValueError('invalid sample count array')

    if locus_counts is None:
      if isinstance(genos,GenotypeArray):
        max_bits = genos.descriptor.max_bit_size
      else:
        max_bits = max(g.model.bit_size for g in genos) if n else 0

      locus_counts = np.zeros((n,1<<max_bits), dtype=int)

    elif locus_counts.shape[0] != n:
      raise ValueError('invalid locus count array')

    for i,geno in enumerate(genos):
      sample_counts[geno.category] += 1
      locus_counts[i,geno.index]   += 1

    sample_counts = np.asarray(sample_counts,dtype=int)

    return sample_counts,locus_counts


  def pick(row, indices):
    '''
    Pick a list of genotypes from a list of indices

    @param      row: genotype sequence
    @type       row: sequence
    @param  indices: sequence of indices
    @type   indices: sequence
    @return        : sequence with elements from the specified indices
    @rtype         : sequence

    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*4)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes

    >>> pick(model.genotypes, [3,2,0])
    [('B', 'B'), ('A', 'B'), (None, None)]

    >>> pick(GenotypeArray(descr,model.genotypes),[3,2,0])
    [('B', 'B'), ('A', 'B'), (None, None)]
    '''
    if isinstance(row,GenotypeArray):
      return row[indices]
    return [ row[i] for i in indices ]


  def pick_columns(rows, index=None):
    '''
    Pick a list of genotypes from a list of indices

    @param   rows: sequence of sequences
    @type    rows: sequence
    @param  index: index
    @type   index: int
    @return      : list of items from the index'th column
    @rtype       : list

    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],0)
    [0, 3, 6]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],1)
    [1, 4, 7]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],2)
    [2, 5, 8]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],[0,2])
    [[0, 3, 6], [2, 5, 8]]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]])
    [[0, 3, 6], [1, 4, 7], [2, 5, 8]]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(2))
    [[0, 3, 6], [1, 4, 7]]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(1,2))
    [[1, 4, 7]]
    >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(0,3,2))
    [[0, 3, 6], [2, 5, 8]]
    '''
    if not rows:
      return []
    elif isinstance(index, int):
      return [ row[index] for row in rows ]
    else:
      if index == None:
        index = xrange(len(rows[0]))
      elif isinstance(index, slice):
        index = xrange(*index.indices(len(rows[0])))
      return [ [ row[i] for row in rows ] for i in index ]


  def place(dest, src, indices):
    '''
    Place items from a sequence at indices into a destination sequence

    @param     dest: destination sequence
    @type      dest: sequence
    @param      src: source sequence
    @type       src: sequence
    @param  indices: sequence of indices
    @type   indices: sequence
    @return        : dest updated with elements of src
    @rtype         : sequence

    >>> model = build_model('AB')
    >>> descr = GenotypeArrayDescriptor([model]*4)
    >>> model.genotypes
    [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
    >>> genos = model.genotypes

    >>> place([None]*4, model.genotypes, [3,2,0,1])
    [('A', 'B'), ('B', 'B'), ('A', 'A'), (None, None)]

    >>> place([None]*4, GenotypeArray(descr,model.genotypes),[3,2,0,1])
    [('A', 'B'), ('B', 'B'), ('A', 'A'), (None, None)]
    '''
    if len(src) != len(indices):
      raise ValueError('source data and index length mismatch')

    for v,i in izip(src,indices):
      dest[i] = v

    return dest


  def place_list(dest, src, indices):
    '''
    Place items from a sequence at indices into a destination sequence
    concatenating into lists if items are not None

    @param     dest: destination sequence
    @type      dest: sequence
    @param      src: source sequence
    @type       src: sequence
    @param  indices: sequence of indices
    @type   indices: sequence
    @return        : dest updated with elements of src
    @rtype         : sequence

    >>> place_list([None]*4, ['A',['A','B'],None,['A']], [2,2,0,1])
    [None, 'A', ['A', 'B', 'A'], None]
    '''
    if len(src) != len(indices):
      raise ValueError('source data and index length mismatch')

    for s,i in izip(src,indices):
      d = dest[i]

      if s is None:
        continue

      s_list = isinstance(s,list)
      d_list = isinstance(d,list)

      if d is None:
        if s_list and len(s)==1:
          dest[i] = s[0]
        else:
          dest[i] = s
      elif s_list and d_list:
        dest[i] += s
      elif d_list:
        dest[i].append(s)
      elif s_list:
        s.append(d)
        dest[i] = s
      else:
        dest[i] = [d,s]

    return dest


###############################################################################


def count_genotypes2(genos1,genos2):
  '''
  Count the two-locus genotypes belonging to the two specified models and
  genotype sequences.  Counts are returned as a list of integers m*n
  integers corresponding to i*n+j, where i,n and j,m are the genotype index
  and total number of genotypes in model1 and model2, respectively.

  @param model1: model for genos1 genotypes
  @type  model1: UnphasedMarkerModel
  @param genos1: genotypes at first locus
  @type  genos1: sequence of Genotype instances
  @param model2: model for genos2 genotypes
  @type  model2: UnphasedMarkerModel
  @param genos2: genotypes at second locus
  @type  genos2: sequence of Genotype instances
  @return      : count of each genotype
  @rtype       : sequence of two-locus genotype and integer count

  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*1400)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos1 = model.genotypes*350
  >>> genos2 = model.genotypes*200+model.genotypes[1:]*200

  >>> for (g1,g2),n in count_genotypes2(genos1,genos2):
  ...   print g1,g2,n
  (None, None) (None, None) 200
  (None, None) ('A', 'A') 50
  (None, None) ('A', 'B') 50
  (None, None) ('B', 'B') 50
  ('A', 'A') (None, None) 0
  ('A', 'A') ('A', 'A') 250
  ('A', 'A') ('A', 'B') 50
  ('A', 'A') ('B', 'B') 50
  ('A', 'B') (None, None) 0
  ('A', 'B') ('A', 'A') 50
  ('A', 'B') ('A', 'B') 250
  ('A', 'B') ('B', 'B') 50
  ('B', 'B') (None, None) 0
  ('B', 'B') ('A', 'A') 50
  ('B', 'B') ('A', 'B') 50
  ('B', 'B') ('B', 'B') 250

  >>> genos1 = GenotypeArray(descr,genos1)
  >>> genos2 = GenotypeArray(descr,genos2)
  >>> for (g1,g2),n in count_genotypes2(genos1,genos2):
  ...   print g1,g2,n
  (None, None) (None, None) 200
  (None, None) ('A', 'A') 50
  (None, None) ('A', 'B') 50
  (None, None) ('B', 'B') 50
  ('A', 'A') (None, None) 0
  ('A', 'A') ('A', 'A') 250
  ('A', 'A') ('A', 'B') 50
  ('A', 'A') ('B', 'B') 50
  ('A', 'B') (None, None) 0
  ('A', 'B') ('A', 'A') 50
  ('A', 'B') ('A', 'B') 250
  ('A', 'B') ('B', 'B') 50
  ('B', 'B') (None, None) 0
  ('B', 'B') ('A', 'A') 50
  ('B', 'B') ('A', 'B') 50
  ('B', 'B') ('B', 'B') 250
  '''
  if not genos1 or not genos2:
    return []

  if len(genos1) != len(genos2):
    raise ValueError("genotype vector sizes do not match: %zd != %zd" % (len(genos1),len(genos2)))

  model1 = genos1[0].model
  model2 = genos2[0].model

  n,m = len(model1.genotypes),len(model2.genotypes)
  counts = [0]*n*m
  for geno1,geno2 in izip(genos1,genos2):
    counts[geno1.index*n+geno2.index] += 1
  genos = ( (g1,g2) for g1 in model1.genotypes
                    for g2 in model2.genotypes )
  return zip(genos,counts)


def locus_genotype_completion_rate(genocounts):
  '''
  Estimate the genotype completion rate, the proportion of non-missing
  genotypes among all genotypes observed, for a given set of genotype counts
  derived from a single locus model.

  # FIXME: Add parameters and doctests
  '''
  n = sum(genocounts)
  if not n:
    return 0.0
  return 1-genocounts[0]/n


def locus_genotype_missing_rate(genocounts):
  '''
  Estimate the genotype completion rate, the proportion of missing genotypes
  among all genotypes observed, for a given set of genotype counts derived
  from a single locus model.

  # FIXME: Add parameters and doctests
  '''
  n = sum(genocounts)
  if not n:
    return 1.0
  return genocounts[0]/n


def count_alleles_from_genocounts(model,genocounts):
  '''
  Count the number of occurrences of each allele belonging to the specified
  locus model given a set of genotype counts.  Counts are returned as a list
  of integers corresponding to the number of each allele in model.alleles
  observed.

  @param model: model for all genotypes
  @type  model: UnphasedMarkerModel
  @param genos: sequence of genotype objects belonging to model
  @type  genos: sequence of Genotype instances
  @return     : count of each genotype
  @rtype      : list of integers

  >>> import random
  >>> model = build_model('AB')
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes*200+model.genotypes[2:]*200
  >>> random.shuffle(genos)

  >>> count_alleles_from_genocounts(model,count_genotypes(genos))
  [400, 800, 1200]
  '''
  counts = [0]*len(model.alleles)
  for geno,n in izip_exact(model.genotypes,genocounts):
    counts[geno.allele1_index] += n
    counts[geno.allele2_index] += n
  return counts


def count_alleles_from_genos(model,genos):
  '''
  Count the number of occurrences of each allele belonging to the specified
  locus model given a set of genotype counts.  Counts are returned as a list
  of integers corresponding to the number of each allele in model.alleles
  observed.

  @param model: model for all genotypes
  @type  model: UnphasedMarkerModel
  @param genos: sequence of genotype objects belonging to model
  @type  genos: sequence of Genotype instances
  @return     : count of each genotype
  @rtype      : list of integers

  >>> import random
  >>> model = build_model('AB')
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes*200+model.genotypes[2:]*200
  >>> random.shuffle(genos)

  >>> count_alleles_from_genos(model,genos)
  [400, 800, 1200]
  '''
  return count_alleles_from_genocounts(model,count_genotypes(genos))


def minor_allele_from_allelecounts(model,allelecounts):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def alleles_from_genos(nn,aa,ab,bb):
  ...   return count_alleles_from_genocounts(model,count_genotypes([NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb))
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(0,1,2,1))
  ('A', 0.5)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(0,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(10000,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(10000,0,2000,2000))
  ('A', 0.25)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(10000,2000,2000,0))
  ('B', 0.25)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(0,1000,0,0))
  ('B', 0.0)
  >>> minor_allele_from_allelecounts(model, alleles_from_genos(1000,0,0,0))
  (None, 0.0)
  '''
  if len(model.alleles) != len(allelecounts):
    raise ValueError('allele counts to not match model alleles')

  allelecounts = [ (n,a) for a,n in izip(model.alleles[1:],allelecounts[1:]) if n ]
  n = len(allelecounts)

  if n > 3:
    raise ValueError('minor allele frequency is defined only for biallelic loci')
  elif not n:
    return None,0.0
  elif n==1:
    b = None
    # If the model has only two possible alleles, choose the other one
    if len(model.alleles)==3:
      a = allelecounts[0][1]
      if a==model.alleles[1]:
        b = model.alleles[2]
      else:
        b = model.alleles[1]
    return b,0.0

  m = sum(n for n,a in allelecounts)
  f,a = min(allelecounts)
  return a,f/m


def minor_allele_from_genocounts(model,counts):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def mkgenos(nn,aa,ab,bb):
  ...   g =[NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb
  ...   return count_genotypes(g)
  >>> minor_allele_from_genocounts(model,mkgenos(0,1,2,1))
  ('A', 0.5)
  >>> minor_allele_from_genocounts(model,mkgenos(0,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_genocounts(model,mkgenos(10000,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_genocounts(model,mkgenos(10000,0,2000,2000))
  ('A', 0.25)
  >>> minor_allele_from_genocounts(model,mkgenos(10000,2000,2000,0))
  ('B', 0.25)
  >>> minor_allele_from_genocounts(model,mkgenos(0,1000,0,0))
  ('B', 0.0)
  >>> minor_allele_from_genocounts(model,mkgenos(1000,0,0,0))
  (None, 0.0)
  '''
  counts = count_alleles_from_genocounts(model,counts)
  return minor_allele_from_allelecounts(model, counts)


def minor_allele_from_genos(genos):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def mkgenos(nn,aa,ab,bb):
  ...   return [NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb
  >>> minor_allele_from_genos(mkgenos(0,1,2,1))
  ('A', 0.5)
  >>> minor_allele_from_genos(mkgenos(0,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_genos(mkgenos(10000,1000,2000,1000))
  ('A', 0.5)
  >>> minor_allele_from_genos(mkgenos(10000,0,2000,2000))
  ('A', 0.25)
  >>> minor_allele_from_genos(mkgenos(10000,2000,2000,0))
  ('B', 0.25)
  >>> minor_allele_from_genos(mkgenos(0,1000,0,0))
  ('B', 0.0)
  >>> minor_allele_from_genos(mkgenos(1000,0,0,0))
  (None, 0.0)
  '''
  if not genos:
    return (None, 0.0)

  model = genos[0].model

  counts = count_alleles_from_genocounts(model,count_genotypes(genos))
  return minor_allele_from_allelecounts(model, counts)


def major_allele_from_allelecounts(model,allelecounts):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def alleles_from_genos(nn,aa,ab,bb):
  ...   return count_alleles_from_genocounts(model,count_genotypes([NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb))
  >>> major_allele_from_allelecounts(model, alleles_from_genos(0,1,2,1))
  ('B', 0.5)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(0,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(10000,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(10000,0,2000,2000))
  ('B', 0.75)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(10000,2000,2000,0))
  ('A', 0.75)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(0,1000,0,0))
  ('A', 1.0)
  >>> major_allele_from_allelecounts(model, alleles_from_genos(1000,0,0,0))
  (None, 0.0)
  '''
  if len(model.alleles) != len(allelecounts):
    raise ValueError('allele counts to not match model alleles')

  allelecounts = [ (n,a) for a,n in izip(model.alleles[1:],allelecounts[1:]) if n ]
  n = len(allelecounts)

  if n > 3:
    raise ValueError('minor allele frequency is defined only for biallelic loci')
  elif not n:
    return None,0.0

  m = sum(n for n,a in allelecounts)
  f,a = max(allelecounts)
  return a,f/m


def major_allele_from_genocounts(model,counts):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def mkgenos(nn,aa,ab,bb):
  ...   g =[NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb
  ...   return count_genotypes(g)
  >>> major_allele_from_genocounts(model,mkgenos(0,1,2,1))
  ('B', 0.5)
  >>> major_allele_from_genocounts(model,mkgenos(0,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_genocounts(model,mkgenos(10000,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_genocounts(model,mkgenos(10000,0,2000,2000))
  ('B', 0.75)
  >>> major_allele_from_genocounts(model,mkgenos(10000,2000,2000,0))
  ('A', 0.75)
  >>> major_allele_from_genocounts(model,mkgenos(0,1000,0,0))
  ('A', 1.0)
  >>> major_allele_from_genocounts(model,mkgenos(1000,0,0,0))
  (None, 0.0)
  '''
  counts = count_alleles_from_genocounts(model,counts)
  return major_allele_from_allelecounts(model, counts)


def major_allele_from_genos(genos):
  '''
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def mkgenos(nn,aa,ab,bb):
  ...   return [NN]*nn+[AA]*aa+[AB]*ab+[BB]*bb
  >>> major_allele_from_genos(mkgenos(0,1,2,1))
  ('B', 0.5)
  >>> major_allele_from_genos(mkgenos(0,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_genos(mkgenos(10000,1000,2000,1000))
  ('B', 0.5)
  >>> major_allele_from_genos(mkgenos(10000,0,2000,2000))
  ('B', 0.75)
  >>> major_allele_from_genos(mkgenos(10000,2000,2000,0))
  ('A', 0.75)
  >>> major_allele_from_genos(mkgenos(0,1000,0,0))
  ('A', 1.0)
  >>> major_allele_from_genos(mkgenos(1000,0,0,0))
  (None, 0.0)
  '''
  if not genos:
    return (None, 0.0)

  model = genos[0].model

  counts = count_alleles_from_genocounts(model,count_genotypes(genos))
  return major_allele_from_allelecounts(model, counts)

##################################################################################


def genotype_count_matrix(genos):
  '''
  Count the number of occurrences of each n genotypes for a series of loci.
  Genotype counts are returned as an m x n matrix of genotype counts where m
  is the number of loci and n-1 is the maximum number of genotypes at any
  locus, and each row contains integers counts corresponding to genotype
  index.  The first column corresponds to the count of missing genotypes.

  @param         genos: genotype stream
  @type          genos: GenotypeStream instance
  @return             : tuple of loci, samples, and m x n matrix of genotype counts
  @rtype              : list of str, set of str, same type as np.zeros( (len(genos.loci),max(genosize)), dtype=int)
  '''
  if genos.format not in ('sdat','ldat'):
    genos = genos.as_sdat()

  if genos.format == 'ldat':
    locus_counts = []
    loci = []
    for lname,geno in genos:
      loci.append(lname)

      count = count_genotypes(geno)

      # FIXME: Keep a running maxgeno instead of assuming allelic SNPs
      if len(count) < 4:
        count = count.tolist()+[0]*(4-len(count))

      locus_counts.append(count)

    samples = set(genos.samples)
    locus_counts = np.asarray(locus_counts,dtype=int)

  else:
    samples = set()
    loci = genos.loci
    locus_counts = None
    for sample,geno in genos:
      sample_count,locus_counts = sample_summary(geno,locus_counts)
      samples.add(sample)

  return loci,samples,locus_counts


##################################################################################

class ModelBuilder(object):
  '''
  Build and cache genotype model objects that may be requested by any:
    * allele list (unordered)
    * genotype list (ordered)
    * maximum number of alleles storable (determines bit width)
    * flag to allow hemizygote genotypes (boolean)

  Any combination of alleles, genotypes, and maximum number of alleles may
  be specified, provided they are consistent and at least one is non-missing
  or empty.  Alleles are treated as unordered sets, while genotypes are
  assumed to be ordered based on the desired encoding.

  Model objects are returned fully populated with all possible genotypes
  assigned an encoding index (see Genotype class).  Caching is extremely
  aggressive, since it is assumed that models are fixed upon creation and no
  new alleles may be added.

  >>> mcache = ModelBuilder()
  >>> m1 = mcache.get(alleles='AB')
  >>> m1.alleles
  [None, 'A', 'B']
  >>> int(m1.max_alleles)
  2
  >>> m2 = mcache.get(alleles='BA',max_alleles=2)
  >>> m2.alleles
  [None, 'A', 'B']
  >>> int(m2.max_alleles)
  2
  >>> m1 is m2
  True
  >>> m3 = mcache.get(alleles='AB',max_alleles=3)
  >>> m1 is not m3
  True
  >>> m4 = mcache.get(genotypes=map(tuple,['AA','AB','BB']))
  >>> m1 is m4
  True
  >>> m5 = mcache.get(genotypes=map(tuple,['AA','BA','BB']))
  >>> m1 is m5
  True
  >>> m6 = mcache.get(genotypes=map(tuple,['BB','BA','AA']))
  >>> m1 is m6
  False
  '''
  def __init__(self):
    self.cache = {}

  def get(self, alleles=None, genotypes=None, max_alleles=None, allow_hemizygote=None, base=None):
    if alleles is None and genotypes is None and not max_alleles:
      raise ValueError('Insufficient information to build genotype model')

    # FASTPATH: Determine the form of key based on mixture of genotype,
    #           allele, and size information

    # Extend base model, if specified, with new alleles and genotypes
    if base is not None:
      if max_alleles is not None and base.max_alleles != max_alleles:
        raise GenotypeRepresentationError('Incompatible model maximum number of alleles')
      if allow_hemizygote is not None and base.allow_hemizygote != allow_hemizygote:
        raise GenotypeRepresentationError('Incompatible model allow_hemizygote')

      if genotypes:
        if not all(g1==g2 for g1,g2 in izip(base.genotypes[1:],genotypes)):
          raise GenotypeRepresentationError('Incompatible model genotypes')
        if len(genotypes) >= len(base.genotypes):
          return base

      alleles          = base.alleles[1:] + list(alleles or [])
      genotypes        = [ g.alleles() for g in base.genotypes[1:] ] + list(genotypes or [])
      max_alleles      = base.max_alleles
      allow_hemizygote = base.allow_hemizygote

    if allow_hemizygote is None:
      allow_hemizygote = False

    if genotypes is not None:
      genotypes = tuple(tuple(sorted(g)) for g in genotypes)

      if alleles is None:
        alleles = (a for g in genotypes for a in g if a is not None)

      alleles = tuple(sorted(set(alleles)))

      if max_alleles is None:
        max_alleles = len(alleles)

      key = (max_alleles,allow_hemizygote)+alleles+genotypes

    else:
      alleles = tuple(sorted(set(alleles)))
      if max_alleles is None:
        max_alleles = len(alleles)
      key = (max_alleles,allow_hemizygote)+alleles

    # Query cache
    model = self.cache.get(key)

    # Handle miss
    if model is None:
      if genotypes and not alleles:
        alleles = tuple(sorted(set(a for g in genotypes for a in g if a is not None)))

      if alleles and not max_alleles:
        max_alleles = len(alleles)

      model = UnphasedMarkerModel(allow_hemizygote,max_alleles)

      if allow_hemizygote:
        alleles = (None,)+alleles

      n = len(alleles)
      genos = ((alleles[i],alleles[j]) for i in range(n) for j in range(i,n))

      if genotypes:
        genos = chain(genotypes,genos)

      for g in genos:
        model.add_genotype(g)

      self.cache[key] = model

      # Store model under alleles with unspecified genotype encoding if no such version exists
      if genotypes:
        key = (max_alleles,allow_hemizygote)+alleles
        if key not in self.cache:
          self.cache[key] = model
      else:
        key += tuple(g.alleles() for g in model.genotypes[1:])
        if key not in self.cache:
          self.cache[key] = model

    return model


class DescrBuilder(object):
  def __init__(self):
    self.cache = {}

  def get_homogeneous(self, model, n):
    key = model,n
    descr = self.cache.get(key)
    if descr is None:
      descr = GenotypeArrayDescriptor( [model]*n )
      self.cache[key] = descr
    return descr

  get_ldat = get_homogeneous

  def get_heterogeneous(self, models):
    return GenotypeArrayDescriptor(models)

  get_sdat = get_heterogeneous


model_builder = ModelBuilder()
build_model = model_builder.get

descr_builder = DescrBuilder()
build_descr = descr_builder.get_homogeneous

##################################################################################

def test_genotypes():
  '''
  >>> model = build_model('AB')
  >>> model.alleles
  [None, 'A', 'B']
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> [ int(g.index)    for g in model.genotypes ]
  [0, 1, 2, 3]
  >>> [ g.category for g in model.genotypes ]
  [0, 2, 3, 2]

  >>> model = build_model('AB',allow_hemizygote=True)
  >>> model.alleles
  [None, 'A', 'B']
  >>> model.genotypes
  [(None, None), (None, 'A'), (None, 'B'), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> [ int(g.index)    for g in model.genotypes ]
  [0, 1, 2, 3, 4, 5]
  >>> [ g.category for g in model.genotypes ]
  [0, 1, 1, 2, 3, 2]

  # Test EQ

  >>> model2 = build_model('AB')
  >>> all(model[geno]==geno for geno in model2.genotypes)
  True
  >>> model['A','B'] == ('B','A')
  True
  >>> model['A','B'] == ('A','B')
  True
  >>> model['B','B'] == ('A','A')
  False
  >>> model['B','B'] == ('A','B')
  False
  >>> model['A','B'] == model2['B','A']
  True
  >>> model['A','B'] == model2['A','B']
  True
  >>> model['B','B'] == model2['A','A']
  False
  >>> model['B','B'] == model2['A','B']
  False

  # Test NE
  >>> model['A','B'] != ('B','A')
  False
  >>> model['A','B'] != ('A','B')
  False
  >>> model['B','B'] != ('A','A')
  True
  >>> model['B','B'] != ('A','B')
  True
  >>> model['A','B'] != model2['B','A']
  False
  >>> model['A','B'] != model2['A','B']
  False
  >>> model['B','B'] != model2['A','A']
  True
  >>> model['B','B'] != model2['A','B']
  True

  # Test hash
  >>> all(hash(model[geno])==hash(geno) for geno in model2.genotypes)
  True
  >>> all( not ((g1==g2) ^ (hash(g1)==hash(g2))) for g1 in model.genotypes for g2 in model2.genotypes)
  True

  # Test contains
  >>> model['A','A'] in model
  True
  >>> ('A','A') in model
  True
  >>> ('A','G') in model
  False
  >>> model2['A','A'] in model
  True

  # Test model compatibility
  >>> model1 = build_model('A')
  >>> model2 = build_model('AB')
  >>> model3 = build_model('AB',allow_hemizygote=True)
  >>> model4 = build_model('AB',genotypes=[('B','B')])

  >>> model1.replaceable_by(model1)
  True
  >>> model2.replaceable_by(model2)
  True
  >>> model1.replaceable_by(model2)
  True
  >>> model2.replaceable_by(model1)
  False
  >>> model2.replaceable_by(model3)
  False
  >>> model1.replaceable_by(model4)
  False
  >>> model2.replaceable_by(model4)
  False
  '''


def test_tolist():
  '''
  >>> def g(genos):
  ...   descr = GenotypeArrayDescriptor([model]*len(genos))
  ...   return GenotypeArray(descr,genos)

  Test 2 bit:

  >>> model = build_model('AB')
  >>> int(model.bit_size)
  2
  >>> NN,AA,AB,BB = model.genotypes
  >>> genos = g([NN,AA,AB,AB,BB,NN])
  >>> g(genos).tolist()
  [(None, None), ('A', 'A'), ('A', 'B'), ('A', 'B'), ('B', 'B'), (None, None)]

  Test 4 bit:

  >>> model = build_model('AB',max_alleles=5)
  >>> int(model.bit_size)
  4
  >>> NN,AA,AB,BB = model.genotypes
  >>> genos = g([NN,AA,AB,AB,BB,NN])
  >>> g(genos).tolist()
  [(None, None), ('A', 'A'), ('A', 'B'), ('A', 'B'), ('B', 'B'), (None, None)]

  >>> model = build_model('AB',max_alleles=16)
  >>> int(model.bit_size)
  8
  >>> NN,AA,AB,BB = model.genotypes
  >>> genos = g([NN,AA,AB,AB,BB,NN])
  >>> g(genos).tolist()
  [(None, None), ('A', 'A'), ('A', 'B'), ('A', 'B'), ('B', 'B'), (None, None)]
  '''


def test_indices():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> genos = g([NN,AA,AB,AB,BB,NN])
  >>> len(genos)
  6
  >>> int(genos.descriptor.homogeneous)
  2
  >>> list(genos.descriptor.offsets)
  [0, 2, 4, 6, 8, 10, 12]
  >>> len(genos.data)
  2
  >>> list(g([NN,AA,AB,AB,BB,NN]).indices())
  [0, 1, 2, 2, 3, 0]

  >>> print [ int(geno.index) for geno in g([NN,AA,AB,AB,BB,NN]) ]
  [0, 1, 2, 2, 3, 0]
  '''


def test_count_genotypes():
  '''
  >>> import random
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*1400)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes*200+model.genotypes[1:]*200
  >>> random.shuffle(genos)
  >>> genos = GenotypeArray(descr,genos)

  Test 2 bit:

  >>> int(model.bit_size)
  2
  >>> genos.counts()
  array([200, 400, 400, 400])
  >>> count_genotypes(genos)
  array([200, 400, 400, 400])
  >>> count_genotypes(genos[:])
  array([200, 400, 400, 400])

  Test 2 bit with tail:

  >>> descr = GenotypeArrayDescriptor([model]*1403)
  >>> genos = genos[:]+model.genotypes[1:]
  >>> random.shuffle(genos)
  >>> genos = GenotypeArray(descr,genos)

  >>> genos.counts()
  array([200, 401, 401, 401])
  >>> count_genotypes(genos)
  array([200, 401, 401, 401])
  >>> count_genotypes(genos[:])
  array([200, 401, 401, 401])
  '''


def test_categories():
  '''
  >>> import random
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*1400)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes*200+model.genotypes[1:]*200
  >>> random.shuffle(genos)

  >>> genotype_categories(genos)
  array([200,   0, 800, 400])
  >>> genotype_categories(GenotypeArray(descr,genos))
  array([200,   0, 800, 400])
  >>> GenotypeArray(descr,genos).categories()
  array([200,   0, 800, 400])
  '''


def test_locus_summary():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes

  >>> l,s = locus_summary(genos)
  >>> l
  array([1, 1, 1, 1])
  >>> s
  array([[1, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 1],
         [0, 0, 1, 0]])
  >>> l,s=locus_summary(GenotypeArray(descr,reversed(genos)),s,l)
  >>> l
  array([2, 2, 2, 2])
  >>> s
  array([[1, 0, 1, 0],
         [0, 0, 1, 1],
         [0, 0, 1, 1],
         [1, 0, 1, 0]])
  '''


def test_sample_summary():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes

  >>> s,l = sample_summary(genos)
  >>> s
  array([1, 0, 2, 1])
  >>> l
  array([[1, 0, 0, 0],
         [0, 1, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 1]])
  >>> s,l=sample_summary(GenotypeArray(descr,reversed(genos)),l,s)
  >>> s
  array([2, 0, 4, 2])
  >>> l
  array([[1, 0, 0, 1],
         [0, 1, 1, 0],
         [0, 1, 1, 0],
         [1, 0, 0, 1]])
  '''


def test_count_genotypes_4bit():
  '''
  >>> model = build_model('AB',max_alleles=5)
  >>> int(model.bit_size)
  4

  Test even number of genotypes

  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> genos = g([NN,AA,AB,AB,BB,NN])
  >>> genos.counts()
  array([2, 1, 2, 1])
  >>> count_genotypes(genos)
  array([2, 1, 2, 1])
  >>> count_genotypes(genos[:])
  array([2, 1, 2, 1])

  Test odd number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*7)
  >>> genos = g([NN,AA,AB,AB,BB,NN,AB])
  >>> genos.counts()
  array([2, 1, 3, 1])
  >>> count_genotypes(genos)
  array([2, 1, 3, 1])
  >>> count_genotypes(genos[:])
  array([2, 1, 3, 1])
  '''


def test_concordance_generic():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)[:]

  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (4, 4)
  >>> genoarray_concordance(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0)
  >>> genoarray_concordance(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (1, 3)
  '''


def test_concordance_4bit():
  '''
  >>> model = build_model('AB',max_alleles=5)
  >>> int(model.bit_size)
  4

  Test even number of genotypes

  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (4, 4)
  >>> genoarray_concordance(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0)
  >>> genoarray_concordance(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (1, 3)

  Test odd number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*7)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN,AB]),g([NN,AA,AB,AB,BB,NN,AB]))
  (5, 5)
  >>> genoarray_concordance(g([NN,NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN,AA]),g([NN,NN,NN,NN,NN,NN,NN]))
  (0, 0)
  >>> genoarray_concordance(g([AA,AB,AB,BB,NN,BB,AB]),g([NN,AA,AB,AB,BB,NN,BB]))
  (1, 4)
  '''


def test_concordance_2bit():
  '''
  >>> model = build_model('AB')
  >>> int(model.bit_size)
  2

  Test even number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (4, 4)
  >>> genoarray_concordance(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0)
  >>> genoarray_concordance(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (1, 3)

  Test odd number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*7)

  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN,AB]),g([NN,AA,AB,AB,BB,NN,AB]))
  (5, 5)
  >>> genoarray_concordance(g([NN,NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0)
  >>> genoarray_concordance(g([NN,AA,AB,AB,BB,NN,AA]),g([NN,NN,NN,NN,NN,NN,NN]))
  (0, 0)
  >>> genoarray_concordance(g([AA,AB,AB,BB,NN,BB,AB]),g([NN,AA,AB,AB,BB,NN,BB]))
  (1, 4)

  Test fractions of a byte

  >>> descr = GenotypeArrayDescriptor([model]*1)
  >>> genoarray_concordance(g([AA]),g([AA]))
  (1, 1)
  >>> descr = GenotypeArrayDescriptor([model]*2)
  >>> genoarray_concordance(g([AA,AB]),g([AA,BB]))
  (1, 2)
  >>> descr = GenotypeArrayDescriptor([model]*3)
  >>> genoarray_concordance(g([AA,AB,BB]),g([AA,BB,AA]))
  (1, 3)
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> genoarray_concordance(g([AA,AB,BB,AA]),g([AA,BB,AA,BB]))
  (1, 4)
  >>> descr = GenotypeArrayDescriptor([model]*5)
  >>> genoarray_concordance(g([AA,AB,BB,AA,AA]),g([AA,BB,AA,BB,AA]))
  (2, 5)
  '''

def test_ibs_generic():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)[:]

  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 4)
  >>> genoarray_ibs(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (0, 2, 1)
  >>> genoarray_ibs(g([AB,BB,BB,BB,AB,AA]),g([AB,AA,AB,AB,BB,BB]))
  (2, 3, 1)
  '''


def test_ibs_4bit():
  '''
  >>> model = build_model('AB',max_alleles=5)
  >>> int(model.bit_size)
  4

  Test even number of genotypes

  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 4)
  >>> genoarray_ibs(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (0, 2, 1)
  >>> genoarray_ibs(g([AB,BB,BB,BB,AB,AA]),g([AB,AA,AB,AB,BB,BB]))
  (2, 3, 1)

  Test odd number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*7)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN,AB]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0, 5)
  >>> genoarray_ibs(g([NN,NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0, 0)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN,AA]),g([NN,NN,NN,NN,NN,NN,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([AA,AB,AB,BB,NN,BB,AA]),g([NN,AA,AB,AB,BB,NN,BB]))
  (1, 2, 1)
  '''


def test_ibs_2bit():
  '''
  >>> model = build_model('AB')
  >>> int(model.bit_size)
  2

  Test even number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*6)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> NN,AA,AB,BB = model.genotypes

  >>> def g(genos):
  ...   return GenotypeArray(descr,genos)

  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 4)
  >>> genoarray_ibs(g([NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN]),g([NN,NN,NN,NN,NN,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([AA,AB,AB,BB,NN,BB]),g([NN,AA,AB,AB,BB,NN]))
  (0, 2, 1)
  >>> genoarray_ibs(g([AB,BB,BB,BB,AB,AA]),g([AB,AA,AB,AB,BB,BB]))
  (2, 3, 1)

  Test odd number of genotypes

  >>> descr = GenotypeArrayDescriptor([model]*7)

  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN,AB]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0, 5)
  >>> genoarray_ibs(g([NN,NN,NN,NN,NN,NN,NN]),g([NN,AA,AB,AB,BB,NN,AB]))
  (0, 0, 0)
  >>> genoarray_ibs(g([NN,AA,AB,AB,BB,NN,AA]),g([NN,NN,NN,NN,NN,NN,NN]))
  (0, 0, 0)
  >>> genoarray_ibs(g([AA,AB,AB,BB,NN,BB,AA]),g([NN,AA,AB,AB,BB,NN,BB]))
  (1, 2, 1)

  Test fractions of a byte

  >>> descr = GenotypeArrayDescriptor([model]*1)
  >>> genoarray_ibs(g([AA]),g([AA]))
  (0, 0, 1)
  >>> descr = GenotypeArrayDescriptor([model]*2)
  >>> genoarray_ibs(g([AA,AB]),g([AA,BB]))
  (0, 1, 1)
  >>> descr = GenotypeArrayDescriptor([model]*3)
  >>> genoarray_ibs(g([AA,AB,BB]),g([AA,BB,AA]))
  (1, 1, 1)
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> genoarray_ibs(g([AA,AB,BB,AA]),g([AA,BB,AA,BB]))
  (2, 1, 1)
  >>> descr = GenotypeArrayDescriptor([model]*5)
  >>> genoarray_ibs(g([AA,AB,BB,AA,AA]),g([AA,BB,AA,BB,AA]))
  (2, 1, 2)
  '''

def test_pick():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes

  >>> pick(model.genotypes, [3,2,0])
  [('B', 'B'), ('A', 'B'), (None, None)]

  >>> pick(GenotypeArray(descr,model.genotypes),[3,2,0])
  [('B', 'B'), ('A', 'B'), (None, None)]

  >>> GenotypeArray(descr,model.genotypes)[ [3,2,0] ]
  [('B', 'B'), ('A', 'B'), (None, None)]
  '''


def test_pick_columns():
  '''
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],0)
  [0, 3, 6]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],1)
  [1, 4, 7]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],2)
  [2, 5, 8]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],[0,2])
  [[0, 3, 6], [2, 5, 8]]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]])
  [[0, 3, 6], [1, 4, 7], [2, 5, 8]]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(2))
  [[0, 3, 6], [1, 4, 7]]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(1,2))
  [[1, 4, 7]]
  >>> pick_columns([[0,1,2],[3,4,5],[6,7,8]],slice(0,3,2))
  [[0, 3, 6], [2, 5, 8]]
  '''


def test_place():
  '''
  >>> model = build_model('AB')
  >>> descr = GenotypeArrayDescriptor([model]*4)
  >>> model.genotypes
  [(None, None), ('A', 'A'), ('A', 'B'), ('B', 'B')]
  >>> genos = model.genotypes

  >>> place([None]*4, model.genotypes, [3,2,0,1])
  [('A', 'B'), ('B', 'B'), ('A', 'A'), (None, None)]

  >>> place([None]*4, GenotypeArray(descr,model.genotypes),[3,2,0,1])
  [('A', 'B'), ('B', 'B'), ('A', 'A'), (None, None)]
  '''


def test_place_list():
  '''
  >>> place_list([None]*4, ['A',['A','B'],None,['A']], [2,2,0,1])
  [None, 'A', ['A', 'B', 'A'], None]
  '''


def bench_locus_summary():
  import time

  m = 1000000
  n = 100 # 5000

  def bench(genos):

    genoslist = genos.tolist()

    sample_counts1 = np.zeros((len(genos),4), dtype=int)
    print sample_counts1.shape
    sample_counts2 = np.zeros((len(genos),4), dtype=int)

    t1 = time.clock()

    for i in xrange(n):
      locus_summary(genos,sample_counts1)

    t2 = time.clock()

    for i in xrange(n):
      locus_summary(genoslist,sample_counts2)

    t3 = time.clock()

    assert (sample_counts1==sample_counts2).all()

    print 'b=%d locus summary(g) = %f' % (descr.homogeneous,t2-t1)
    print 'b=%d locus summary(l) = %f' % (descr.homogeneous,t3-t2)
    print

  genos = [('A','A'),('A','C'),('C','C'),('A','A'),(None,None),(None,'A'),('C',None)]*(m//7)
  model = build_model('ACGT',allow_hemizygote=True)
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)

  model = build_model('AB')
  genos = model.genotypes*(m//len(model.genotypes))
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)


def bench_sample_summary():
  import time

  m = 1000000
  n = 100 # 5000

  def bench(genos):

    genoslist = genos.tolist()

    locus_counts1 = np.zeros((len(genos),1<<genos.descriptor.max_bit_size), dtype=int)
    print locus_counts1.shape
    locus_counts2 = np.zeros((len(genos),1<<genos.descriptor.max_bit_size), dtype=int)

    t1 = time.clock()

    for i in xrange(n):
      sample_summary(genos,locus_counts1)

    t2 = time.clock()

    for i in xrange(n):
      sample_summary(genoslist,locus_counts2)

    t3 = time.clock()

    assert (locus_counts1==locus_counts2).all()

    print 'b=%d sample summary(g) = %f' % (descr.homogeneous,t2-t1)
    print 'b=%d sample summary(l) = %f' % (descr.homogeneous,t3-t2)
    print

  genos = [('A','A'),('A','C'),('C','C'),('A','A'),(None,None),(None,'A'),('C',None)]*(m//7)
  model = build_model('ACGT',allow_hemizygote=True)
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)

  model = build_model('AB')
  genos = model.genotypes*(m//len(model.genotypes))
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)


def bench_categories():
  import time

  m = 10000
  n = 1000

  def bench(genos):

    genoslist = genos.tolist()

    t1 = time.clock()

    for i in xrange(n):
      genos.categories()

    t2 = time.clock()

    for i in xrange(n):
      genotype_categories(genos)

    t3 = time.clock()

    for i in xrange(n):
      genotype_categories(genoslist)

    t4 = time.clock()

    #for i in xrange(n):
    #  genotype_categories2(genos)

    #t5 = time.clock()

    #for i in xrange(n):
    #  genotype_categories2(genoslist)

    #t6 = time.clock()

    print 'b=%d .counts    = %f' % (descr.homogeneous,t2-t1)
    print 'b=%d counts(g)  = %f' % (descr.homogeneous,t3-t2)
    print 'b=%d counts(l)  = %f' % (descr.homogeneous,t4-t3)
    #print 'b=%d counts2(g) = %f' % (descr.homogeneous,t5-t4)
    #print 'b=%d counts2(l) = %f' % (descr.homogeneous,t6-t5)
    print

  genos = [('A','A'),('A','C'),('C','C'),('A','A'),(None,None),(None,'A'),('C',None)]*m
  model = build_model('ACGT',allow_hemizygote=True)
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)

  model = build_model('AB')
  genos = model.genotypes*m*2
  descr = GenotypeArrayDescriptor( [model]*len(genos) )
  genos = GenotypeArray(descr, genos)
  bench(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
  #bench_categories()
  #bench_locus_summary()
  #bench_sample_summary()
