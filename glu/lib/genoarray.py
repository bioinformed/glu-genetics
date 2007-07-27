# -*- coding: utf-8 -*-
'''
File:          genoarray.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import unittest

from   array      import array
from   operator   import getitem
from   itertools  import imap, repeat

from   genoarray2 import Genotype, GenotypeArray, GenotypeArrayDescriptor, model_from_alleles


class Nothing(object): pass

# Glossary of terms:
#   genotype tuple  = tuple of two alleles or None if missing.
#   internal rep    = opaque object representing a genotype,
#                     evaluates to False iff missing
#   genotype string = string representation
#   packed          = random access sequence of internal genotype reps


class UnphasedMarkerRepresentation(object):
  '''
  Generic unphased marker representation where the genotype representation
  and internal representation are the same.
  '''

  def __init__(self, missing_geno_str=Nothing, missing_geno_strs=(),
                     missing_allele_str='',    missing_allele_strs=(),
                     delimiter='/'):
    '''
    Construct a new UnphasedMarkerRepresentation

    This class represents bidirectional mappings of genotypes between
    strings and Python objects.  The object representation of a genotype is
    a tuple of two alleles, with the missing genotype being None or any
    other object that evaluates to False.  Given this representation,
    alleles need not be known in advance.  The "internal representation" is
    the same as the Python object representation (for an example of a non-
    trivial internal representation see UnphasedMarkerRepresentation
    below).

    @param    missing_geno_str: canonical missing genotype output string representation, optional
                                and defautls to Nothing
    @type     missing_geno_str: str
    @param   missing_geno_strs: missing genotype input strings
    @type    missing_geno_strs: sequence
    @param  missing_allele_str: canonical missing allele output string representation, default is ''
                                and defaults to an empty string
    @type   missing_allele_str: str
    @param missing_allele_strs: missing allele input strings
    @type  missing_allele_strs: sequence
    @param           delimiter: allele delimiter for genotypes, defaults is '/'
    @type            delimiter: str
    '''
    self.missing_geno_str    = missing_geno_str
    self.missing_geno_strs   = missing_geno_strs
    self.missing_allele_str  = missing_allele_str
    self.missing_allele_strs = missing_allele_strs
    self.delimiter           = delimiter

    self.allelerep  = {None:self.missing_allele_str}
    self.allelestr  = dict( (m,None) for m in missing_allele_strs )
    self.missingrep = (None,None)
    self.missing    = None
    self.strcache   = {}

  def homozygote_geno(self,g):
    '''
    Test if the genotype is homozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a homozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedMarkerRepresentation('AG')
    >>> model.homozygote_geno( ('A','A') )
    True
    >>> model.homozygote_geno( ('A','G') )
    False
    >>> model.homozygote_geno( ('A',None) )
    False
    >>> model.homozygote_geno(None)
    False
    '''
    return bool(g is not None and g[0]==g[1] and g[0])

  def heterozygote_geno(self,g):
    '''
    Test if the genotype is heterozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a heterozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedMarkerRepresentation('AG')
    >>> model.heterozygote_geno( ('A','A') )
    False
    >>> model.heterozygote_geno( ('A','G') )
    True
    >>> model.heterozygote_geno( ('A',None) )
    False
    >>> model.heterozygote_geno(None)
    False
    '''
    return bool(g is not None and g[0]!=g[1] and g[0] and g[1])

  def hemizygote_geno(self,g):
    '''
    Check if the genotype is hemizygote.
    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a hemizygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedMarkerRepresentation('AG')
    >>> model.hemizygote_geno( ('A','A') )
    False
    >>> model.hemizygote_geno( ('A','G') )
    False
    >>> model.hemizygote_geno( ('A',None) )
    True
    >>> model.hemizygote_geno(None)
    False
    '''
    return bool(g is not None and (bool(g[0]) ^ bool(g[1])))

  # Internal rep is a genotype tuple, so these functions are the same
  homozygote   = homozygote_geno
  heterozygote = heterozygote_geno
  hemizygote   = hemizygote_geno

  def geno_from_rep(self, rep):
    '''
    Transform a genotype in internal format to a genotype tuple

    @param rep: genotype in internal format
    @type  rep: object
    @return   : genotypes in tuple format
    @rtype    : tuple

    >>> snp_marker.geno_from_rep( ('A','A') )
    ('A', 'A')
    >>> snp_marker.geno_from_rep( ('A',None) )
    ('A', None)
    >>> snp_marker.geno_from_rep( ('A','C') )
    ('A', 'C')
    >>> snp_marker.geno_from_rep( (None,None) )
    >>> snp_marker.geno_from_rep( None )
    '''
    if not rep or rep == self.missingrep:
      return self.missing
    return rep

  def genos_from_reps(self, reps):
    '''
    Transform a sequence genotypes in internal format to a sequence of genotype tuples

    @param reps: genotypes in internal format
    @type  reps: sequence
    @return    : genotypes in tuple format
    @rtype     : sequence

    >>> g = [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C'),None,(None,None)]
    >>> snp_marker.genos_from_reps(g)
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C'), None, None]
    '''
    return map(self.geno_from_rep,reps)

  def str_from_rep(self, rep):
    '''
    Transform a genotype in internal format to a genotype string

    @param rep: genotype in internal format
    @type  rep: object
    @return   : genotype in string format
    @rtype    : str

    >>> snp_marker.str_from_rep( ('A','A') )
    'AA'
    >>> snp_marker.str_from_rep( (None,'A') )
    ' A'
    >>> snp_marker.str_from_rep( ('A',None) )
    'A '
    >>> snp_marker.str_from_rep(None)
    '  '
    >>> generic_marker.str_from_rep( ('A','A') )
    'A/A'
    >>> generic_marker.str_from_rep( (None,'A') )
    '/A'
    >>> generic_marker.str_from_rep(None)
    ''
    '''
    if not rep or rep == self.missingrep:
      if self.missing_geno_str is not Nothing:
        return self.missing_geno_str
      rep = self.missing
    return self.delimiter.join(self.allelerep.get(a,a) for a in rep)

  def strs_from_reps(self,reps):
    '''
    Transform a sequence of genotypes in internal format to a genotype strings

    @param rep: genotypes in internal format
    @type  rep: sequence
    @return   : genotypes in string format
    @rtype    : sequence

    >>> snp_marker.strs_from_reps( [('A','A'),(None,'A'),None] )
    ['AA', ' A', '  ']
    >>> generic_marker.strs_from_reps( [('A','A'),(None,'A'),None] )
    ['A/A', '/A', '']
    '''
    return map(self.str_from_rep,reps)

  def rep_from_geno(self,geno):
    '''
    Transform a genotype tuple to internal format

    @param geno: genotype tuple
    @type  geno: tuple
    @return    : genotype in internal format
    @rtype     : object

    >>> snp_marker.rep_from_geno( ('A','A') )
    ('A', 'A')
    >>> snp_marker.rep_from_geno( (None,'A') )
    (None, 'A')
    >>> snp_marker.rep_from_geno(None)
    >>> snp_marker.rep_from_geno( (None,None) )
    '''
    if not geno or geno == (None,None):
      return self.missing
    return geno

  def reps_from_genos(self,genos):
    '''
    Transform a sequence of genotype tuples to internal format

    @param geno: sequence of genotype tuples
    @type  geno: sequence
    @return    : genotypes in internal format
    @rtype     : sequence

    >>> snp_marker.reps_from_genos( [('A','A'),(None,'A'),(None,None),None] )
    [('A', 'A'), (None, 'A'), None, None]
    '''
    return map(self.rep_from_geno,genos)

  def rep_from_str(self, geno):
    '''
    Transform a genotype string to internal format

    @param geno: genotype string
    @type  geno: str
    @return    : genotype in internal format
    @rtype     : object

    >>> snp_marker.rep_from_str('AA')
    ('A', 'A')
    >>> snp_marker.rep_from_str('AB')
    ('A', 'B')
    >>> snp_marker.rep_from_str('BA')
    ('A', 'B')
    >>> snp_marker.rep_from_str(' A')
    (None, 'A')
    >>> snp_marker.rep_from_str('A ')
    (None, 'A')
    >>> snp_marker.rep_from_str('  ')
    >>> snp_marker.rep_from_str('')
    >>> generic_marker.rep_from_str('A/A')
    ('A', 'A')
    >>> generic_marker.rep_from_str('A/B')
    ('A', 'B')
    >>> generic_marker.rep_from_str('B/A')
    ('A', 'B')
    >>> generic_marker.rep_from_str(' /A')
    (None, 'A')
    >>> generic_marker.rep_from_str('/A')
    (None, 'A')
    >>> generic_marker.rep_from_str('A/ ')
    (None, 'A')
    >>> generic_marker.rep_from_str('A/')
    (None, 'A')
    >>> generic_marker.rep_from_str(' / ')
    >>> generic_marker.rep_from_str('/')
    >>> generic_marker.rep_from_str('/ ')
    >>> generic_marker.rep_from_str(' /')
    >>> generic_marker.rep_from_str('')
    '''
    strcache = self.strcache
    if geno in strcache:
      return strcache[geno]

    if geno in self.missing_geno_strs:
      self.strcache[geno] = None
      return None

    if self.delimiter:
      alleles = geno.split(self.delimiter)
    else:
      alleles = geno

    assert len(alleles) == 2

    rep = tuple(sorted(self.allelestr.get(a,a) for a in alleles))

    if not rep or rep == self.missingrep:
      self.strcache[geno] = None
      return None

    self.strcache[geno] = rep
    return rep

  def reps_from_strs(self,genos):
    '''
    Transform a sequence of genotype strings to internal format

    @param geno: sequence of genotype strings
    @type  geno: sequence of str
    @return    : sequence of genotype in internal format
    @rtype     : sequence of objects

    >>> snp_marker.reps_from_strs(['AA','AB','BA',' A','A ','  ',''])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), None, None]
    >>> generic_marker.reps_from_strs(['A/A','A/B','B/A',' /A','A/ ',' / ',''])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), None, None]
    '''
    return map(self.rep_from_str, genos)

  def pack_genos(self, genos):
    '''
    Transform a sequence of genotype tuples into a packed internal format

    @param geno: sequence of genotype strings
    @type  geno: sequence of str
    @return    : packed genotype in internal format
    @rtype     : packed sequence of objects

    >>> snp_marker.pack_genos( [('A','A'),('A','B'),('A','B'),(None,'A'),(None,'A'),None,None] )
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), None, None]
    '''
    return list(genos)

  def pack_reps(self, reps):
    '''
    Transform a sequence of genotype tuples into a packed internal format

    @param geno: sequence of genotype tuples
    @type  geno: sequence of str
    @return    : packed sequence of genotype in internal format
    @rtype     : packed sequence of objects

    >>> snp_marker.pack_reps([('A','A'),('A','B'),('A','B'),(None,'A'),(None,'A'),None,None])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), None, None]
    '''
    return list(reps)

  def pack_strs(self, genos):
    '''
    Transform a sequence of genotype strings into a packed internal format

    @param geno: sequence of genotype string
    @type  geno: sequence of str
    @return    : packed sequence of genotype in internal format
    @rtype     : packed sequence of objects

    >>> snp_marker.pack_strs(['AA','AB','BA',' A','A ','  ',''])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), None, None]
    '''
    return self.reps_from_strs(genos)

  # Backward compatibility with old API
  geno      = geno_from_rep
  genos     = genos_from_reps
  geno_str  = str_from_rep
  genos_str = strs_from_reps
  byte      = rep_from_geno
  bytes     = reps_from_genos
  byte_strs = reps_from_strs
  byte_str  = rep_from_str
  array     = pack_genos
  pack      = pack_reps
  array_str = pack_strs


class UnphasedPackedMarkerRepresentation(object):
  def __init__(self, alleles=(), missing_geno_str=Nothing, missing_geno_strs=(),
                                 missing_allele_str='',    missing_allele_strs=(),
                                 delimiter='/', max_alleles=None):
    '''
    Construct a new UnphasedPackedMarkerRepresentation

    This class represents bidirectional mappings of genotypes between
    strings, Python objects, and integers given a list of alleles, as well
    as string and object missing value representations.  This class can also
    build packed representations of genotype vectors, which utilize
    low-level arrays that are able to greatly reduce the amount of storage
    required.  These representations rely on a fixed single- or multi-byte
    encoding of each allele.  Thus, the genotype size is computed as twice
    the number bits required to store each allele.  Missing alleles are
    encoded as 0, all other alleles are assigned sequentially increasing
    integer values.

    Given n alleles, the width of the genotype is either 1, 2, or 4 bytes
    and a shift value is assigned that is the genotype byte width times 4 (8
    bits/byte divided by 2 alelles per genotype).

    When encoding a genotype, the representation is calulated as follows:
    the integer values of the two alleles are sorted and the smaller is
    shifted left by the shift constant and then combined with the integer
    value of second allele.  Thus, the high order bits store the lesser
    allele and the lower order bits the higher valued allele:

    ||=================== Genotype {1,2,4 bytes} ==================||
    ||--- allele 1 {4,8,16 bits} ---|--- allele 2 {4,8,16 bits} ---||
    ||       min(a1,a2)<<shift      |           max(a1,a2)         ||
    =================================================================

    @param             alleles: possible alleles
    @type              alleles: sequence
    @param    missing_geno_str: canonical missing genotype output string representation, optional
                                and defautls to Nothing
    @type     missing_geno_str: str
    @param   missing_geno_strs: missing genotype input strings
    @type    missing_geno_strs: sequence
    @param  missing_allele_str: canonical missing allele output string representation, default is ''
                                and defaults to an empty string
    @type   missing_allele_str: str
    @param missing_allele_strs: missing allele input strings
    @type  missing_allele_strs: sequence
    @param           delimiter: allele delimiter for genotypes, defaults is '/'
    @type            delimiter: str
    '''
    n = self.max_alleles = max(len(alleles),max_alleles)

    if n < 2:
      n = 255

    if n < 16:
      self.typecode = 'B'
    elif n < 256:
      self.typecode = 'H'
    elif n < 65536:
      self.typecode = 'L'
    else:
      raise ValueError,'Too many alleles for UnphasedPackedMarkerRepresentation: %d' % n

    self.delimiter           = delimiter
    self.alleles             = set()
    self.missing_geno_str    = missing_geno_str
    self.missing_geno_strs   = set()
    self.missing_allele_str  = missing_allele_str
    self.missing_allele_strs = set()
    self.allele_index        = {}
    self.geno_to_rep        = {None:0,(None,None):0}
    self.str_to_rep         = {}
    self.rep_to_geno        = [None]
    self.rep_to_str         = ['']
    self.rep_class          = [0]
    self.missing            = 0

    if missing_geno_str is not Nothing:
      self.add_missing_geno_str(missing_geno_str)

    for g in missing_geno_strs:
      self.add_missing_geno_str(g)

    self.add_missing_allele(missing_allele_str)
    for a in missing_allele_strs:
      self.add_missing_allele(a)

    for a in alleles:
      self.add_allele(a)

  def add_missing_geno_str(self,g):
    if g in self.missing_geno_strs:
      return

    if g in self.alleles or g in self.missing_allele_strs:
      raise ValueError, 'Missing genotype code may not be an allele'

    if g in self.str_to_rep or g in self.geno_to_rep:
      raise ValueError, 'Ambiguous missing genotype code specified'

    self.str_to_rep[g] = 0

  def _add_allele(self,a,n):
    self.allele_index[a] = n

    def _mis(a):
      if a in self.missing_allele_strs:
        return None
      return a

    if n:
      b = self._rep_rep(n,0)
      self.geno_to_rep[(a,None)] = b
      self.geno_to_rep[(None,a)] = b
      self.rep_class[b] = self._rep_class((a,None))

    for a2,b2 in self.allele_index.iteritems():
      g  = tuple(sorted([a,a2]))
      g2 = (g[1],g[0])
      g3 = (_mis(g[0]),_mis(g[1]))
      g4 = (g3[1],g3[0])
      b  = self._rep_rep(n,b2)

      if b:
        self.rep_to_geno[b] = g3

      self.geno_to_rep[g3] = b
      self.geno_to_rep[g4] = b

      self.rep_class[b] = self._rep_class(g3)

      g  = self.delimiter.join(g)
      g2 = self.delimiter.join(g2)
      self.str_to_rep[g]  = b
      self.str_to_rep[g2] = b

      if b or self.missing_geno_str is Nothing:
        self.rep_to_str[b] = g

  def _expand(self):
    # Expand the indexed arrays to accommodate a new allele
    # Let the amount we must expand be:
    # m =   new size   - old size
    #   = (n<<shift)+n - [ ((n-1)<<shift) + (n-1) ]
    #   = ((n-n+1)<<shift) + n-n+1
    #   = (1<<shift) + 1
    shift = self._shift[self.typecode]
    m = (1<<shift)+1
    self.rep_to_geno.extend([None]*m)
    self.rep_class.extend([0]*m)
    self.rep_to_str.extend([None]*m)

  def add_missing_allele(self,a):
    if a in self.missing_allele_strs:
      return

    if a in self.alleles:
      raise NotImplementedError, 'Setting an existing allele to missing not yet supported'

    self.missing_allele_strs.add(a)
    self._add_allele(a,0)

  def add_allele(self,a):
    if a in self.alleles:
      return

    if a in self.missing_allele_strs:
      raise NotImplementedError, 'Setting an existing missing allele to non-missing not yet supported'

    self.alleles.add(a)
    n = len(self.alleles)
    self._expand()
    self._add_allele(a,n)

  def _rep_class(self,g):
    if self.homozygote_geno(g):
      return 1
    elif self.heterozygote_geno(g):
      return 2
    elif self.hemizygote_geno(g):
      return 3
    else:
      return 0

  _shift = {'B':4,'H':8,'L':16}

  def _rep_rep(self,b1,b2):
    b1,b2 = min(b1,b2),max(b1,b2)
    shift = self._shift[self.typecode]
    return b1<<shift|b2

  def homozygote_geno(self,g):
    '''
    Test if the genotype is homozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a homozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.homozygote_geno( ('A','A') )
    True
    >>> model.homozygote_geno( ('A','G') )
    False
    >>> model.homozygote_geno( ('A',None) )
    False
    >>> model.homozygote_geno(None)
    False
    '''
    return bool(g is not None and g[0]==g[1] and g[0] is not None)

  def heterozygote_geno(self,g):
    '''
    Test if the genotype is heterozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a heterozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.heterozygote_geno( ('A','A') )
    False
    >>> model.heterozygote_geno( ('A','G') )
    True
    >>> model.heterozygote_geno( ('A',None) )
    False
    >>> model.heterozygote_geno(None)
    False
    '''
    return g is not None and None!=g[0]!=g[1]!=None

  def hemizygote_geno(self,g):
    '''
    Check if the genotype is hemizygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a hemizygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.hemizygote_geno( ('A','A') )
    False
    >>> model.hemizygote_geno( ('A','G') )
    False
    >>> model.hemizygote_geno( ('A',None) )
    True
    >>> model.hemizygote_geno(None)
    False
    '''
    return bool(g is not None and (bool(g[0]) ^ bool(g[1])))

  def homozygote(self,rep):
    '''
    Check if the encoded genotype is homozygote.

    @param rep: encoded genotype
    @type  rep: int
    @return   : True if homozygote, False otherwise
    @rtype    : bool

    >>> g = snp_acgt.rep_from_geno(('A','T'))
    >>> snp_acgt.homozygote(g)
    False
    >>> g = snp_acgt.rep_from_geno(('A','A'))
    >>> snp_acgt.homozygote(g)
    True
    >>> g = snp_acgt.rep_from_geno((None,'A'))
    >>> snp_acgt.homozygote(g)
    False
    '''
    return self.rep_class[rep] == 1

  def heterozygote(self,rep):
    '''
    Check if the encoded genotype is heterozygote.

    @param rep: the encoded genotype
    @type  rep: int
    @return    : True if heterozygote, False otherwise
    @rtype     : bool

    >>> g = snp_acgt.rep_from_geno(('A','T'))
    >>> snp_acgt.heterozygote(g)
    True
    >>> g = snp_acgt.rep_from_geno(('A','A'))
    >>> snp_acgt.heterozygote(g)
    False
    >>> g = snp_acgt.rep_from_geno((None,'A'))
    >>> snp_acgt.heterozygote(g)
    False
    '''
    return self.rep_class[rep] == 2

  def hemizygote(self,rep):
    '''
    Check if the encoded genotype is hemizygote.

    @param rep: the encoded genotype
    @type  rep: int
    @return    : True if hemizygote, False otherwise
    @rtype     : bool

    >>> g = snp_acgt.rep_from_geno((None,'A'))
    >>> snp_acgt.hemizygote(g)
    True
    >>> g = snp_acgt.rep_from_geno(('A','A'))
    >>> snp_acgt.hemizygote(g)
    False
    >>> g = snp_acgt.rep_from_geno(('A','G'))
    >>> snp_acgt.hemizygote(g)
    False
    '''
    return self.rep_class[rep] == 3

  def rep_from_geno(self,geno):
    '''
    Calculate encoded representation for the genotype.

    @param geno: the genotype
    @type  geno: tuple
    @return    : encoded genotype
    @rtype     : int

    >>> snp_acgt.rep_from_geno( (None, 'A') )
    1
    >>> snp_acgt.rep_from_geno( ('A', None) )
    1
    >>> snp_acgt.rep_from_geno(None)
    0
    >>> snp_acgt.rep_from_geno( (None,None) )
    0
    '''
    return self.geno_to_rep[geno]

  def reps_from_genos(self,genos):
    '''
    Calculate the encoded representation for the genotypes.

    @param genos: genotypes
    @type  genos: sequence
    @return     : encoded genotypes
    @rtype      : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> for rep in snp_acgt.reps_from_genos(g):
    ...   print rep
    17
    18
    18
    34
    '''
    return list(imap(getitem, repeat(self.geno_to_rep), genos))

  def pack_genos(self,genos):
    '''
    Build a packed array (compact representation) for the genotypes.

    @param genos: genotypes
    @type  genos: iterable
    @return     : genotypes represented in a encoded array
    @rtype      : array

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt.pack_genos(g)
    array('B', [17, 18, 18, 34])
    '''
    return array(self.typecode, imap(getitem, repeat(self.geno_to_rep), genos) )

  def rep_from_str(self,geno):
    '''
    Calculate encoded representation for the string genotype.

    @param geno: the genotype
    @type  geno: str
    @return    : encoded genotype
    @rtype     : int

    >>> snp_acgt.rep_from_str(' A')
    1
    >>> snp_acgt.rep_from_str('  ')
    0
    '''
    return self.str_to_rep[geno]

  def reps_from_strs(self,genos):
    '''
    Calculate the encoded representation for the string genotypes.

    @param genos: genotypes
    @type  genos: sequence of str
    @return     : encoded genotypes
    @rtype      : iterator

    >>> g = ['AA','AC','CA','CC']
    >>> for rep in snp_acgt.reps_from_strs(g):
    ...   print rep
    17
    18
    18
    34
    '''
    return list(imap(getitem, repeat(self.str_to_rep), genos))

  def pack_strs(self,genos):
    '''
    Build a packed array (compact representation) for the string genotypes.

    @param genos: sequence of str
    @type  genos: iterable
    @return     : encoded genotypes
    @rtype      : array

    >>> g = ['AA','AC','CA','CC']
    >>> snp_acgt.pack_strs(g)
    array('B', [17, 18, 18, 34])
    '''
    return array(self.typecode, imap(getitem, repeat(self.str_to_rep), genos) )

  def geno_from_rep(self,rep):
    '''
    Build the genotype tuple from an encoded genotype

    @param rep: encoded genotype
    @type  rep: object
    @return   : genotypes in tuple format
    @rtype    : tuple

    >>> snp_acgt.geno_from_rep(17)
    ('A', 'A')
    >>> snp_acgt.geno_from_rep(18)
    ('A', 'C')
    >>> snp_marker.geno_from_rep(0)
    '''
    return self.rep_to_geno[rep]

  def str_from_rep(self,rep):
    '''
    Build the genotype string out of encoded genotype

    @param rep: genotypes in encoded format
    @type  rep: object
    @return   : genotypes in string format
    @rtype    : str

    >>> snp_acgt.str_from_rep(17)
    'AA'
    >>> snp_acgt.str_from_rep(18)
    'AC'
    '''
    return self.rep_to_str[rep]

  def genos_from_reps(self,reps):
    '''
    Build the genotype tuples from encoded genotypes

    @param reps: genotypes in encoded format
    @type  reps: sequence
    @retur     : genotypes in tuple format
    @rtype     : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt.genos_from_reps(snp_acgt.pack_genos(g))
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    '''
    return list(imap(getitem, repeat(self.rep_to_geno), reps))

  def strs_from_reps(self,reps):
    '''
    Build the genotype strings out of encoded genotypes

    @param reps: genotypes in internal format
    @type  reps: sequence
    @return    : genotypes in string format
    @rtype     : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt.strs_from_reps(snp_acgt.pack_genos(g))
    ['AA', 'AC', 'AC', 'CC']
    '''
    return list(imap(getitem, repeat(self.rep_to_str), reps))

  def pack_reps(self,reps):
    '''
    Pack the reps into an array.

    @param reps: a sequence of encoded genotypes
    @type  reps: sequence
    @return    : a rep array (i.e. reps in compact form)
    @rtype     : array

    >>> snp_acgt.pack_reps( [17,18,18,34] )
    array('B', [17, 18, 18, 34])
    '''
    return array(self.typecode,reps)

  # Backward compatibility with old API
  geno      = geno_from_rep
  genos     = genos_from_reps
  geno_str  = str_from_rep
  genos_str = strs_from_reps
  byte      = rep_from_geno
  bytes     = reps_from_genos
  byte_strs = reps_from_strs
  byte_str  = rep_from_str
  array     = pack_genos
  pack      = pack_reps
  array_str = pack_strs


class UnphasedGenotypeArrayRepresentation(object):
  def __init__(self, alleles=(), missing_geno_str=Nothing, missing_geno_strs=(),
                                 missing_allele_str='',    missing_allele_strs=(),
                                 delimiter='/', max_alleles=None):
    '''
    Construct a new UnphasedGenotypeArrayRepresentation

    This class represents bidirectional mappings of genotypes between
    strings, and genotype tuples, as well as string and object missing value
    representations.  This class can also build packed representations of
    genotype vectors, which utilizes low-level bit-packed arrays that are able to
    greatly reduce the amount of storage required.

    @param             alleles: possible alleles
    @type              alleles: sequence
    @param    missing_geno_str: canonical missing genotype output string representation, optional
                                and defautls to Nothing
    @type     missing_geno_str: str
    @param   missing_geno_strs: missing genotype input strings
    @type    missing_geno_strs: sequence
    @param  missing_allele_str: canonical missing allele output string representation, default is ''
                                and defaults to an empty string
    @type   missing_allele_str: str
    @param missing_allele_strs: missing allele input strings
    @type  missing_allele_strs: sequence
    @param           delimiter: allele delimiter for genotypes, defaults is '/'
    @type            delimiter: str
    '''
    self.model = model_from_alleles(alleles, allow_hemizygote=True)
    self.missing = self.model.get_genotype( (None,None) )

    self.delimiter           = delimiter
    self.alleles             = set()
    self.missing_geno_str    = missing_geno_str
    self.missing_geno_strs   = set()
    self.missing_allele_str  = missing_allele_str
    self.missing_allele_strs = set()
    self.allele_index        = {}
    self.geno_to_rep         = {None:self.missing,(None,None):self.missing}
    self.str_to_rep          = {}
    self.rep_to_str          = {}
    self.rep_to_geno         = {self.missing:None}
    self._descr_cache        = {}

    if missing_geno_str is not Nothing:
      self.add_missing_geno_str(missing_geno_str)

    for g in missing_geno_strs:
      self.add_missing_geno_str(g)

    self.add_missing_allele(missing_allele_str)
    for a in missing_allele_strs:
      self.add_missing_allele(a)

    for a in alleles:
      self.add_allele(a)

  def add_missing_geno_str(self,g):
    if g in self.missing_geno_strs:
      return

    if g in self.alleles or g in self.missing_allele_strs:
      raise ValueError, 'Missing genotype code may not be an allele'

    if g in self.str_to_rep or g in self.geno_to_rep:
      raise ValueError, 'Ambiguous missing genotype code specified'

    self.str_to_rep[g] = self.model.get_genotype( (None,None) )

  def _add_allele(self,a,n):
    self.allele_index[a] = n

    def _mis(a):
      if a in self.missing_allele_strs:
        return None
      return a

    for a2,b2 in self.allele_index.iteritems():
      g  = tuple(sorted([a,a2]))
      g2 = (g[1],g[0])
      g3 = (_mis(g[0]),_mis(g[1]))
      g4 = (g3[1],g3[0])
      b  = self.model.add_genotype(g3)

      self.geno_to_rep[g3] = b
      self.geno_to_rep[g4] = b

      if b:
        self.rep_to_geno[b]  = g

      g  = self.delimiter.join(g)
      g2 = self.delimiter.join(g2)
      self.str_to_rep[g]  = b
      self.str_to_rep[g2] = b

      if b or self.missing_geno_str is Nothing:
        self.rep_to_str[b] = g

  def add_missing_allele(self,a):
    if a in self.missing_allele_strs:
      return

    if a in self.alleles:
      raise NotImplementedError, 'Setting an existing allele to missing not yet supported'

    self.missing_allele_strs.add(a)
    self._add_allele(a,0)

  def add_allele(self,a):
    if a in self.alleles:
      return

    if a in self.missing_allele_strs:
      raise NotImplementedError, 'Setting an existing missing allele to non-missing not yet supported'

    self.alleles.add(a)
    n = len(self.alleles)
    self._add_allele(a,n)

  def homozygote_geno(self,g):
    '''
    Test if the genotype is homozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a homozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.homozygote_geno( ('A','A') )
    True
    >>> model.homozygote_geno( ('A','G') )
    False
    >>> model.homozygote_geno( ('A',None) )
    False
    >>> model.homozygote_geno(None)
    False
    '''
    return bool(g is not None and g[0]==g[1] and g[0] is not None)

  def heterozygote_geno(self,g):
    '''
    Test if the genotype is heterozygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a heterozygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.heterozygote_geno( ('A','A') )
    False
    >>> model.heterozygote_geno( ('A','G') )
    True
    >>> model.heterozygote_geno( ('A',None) )
    False
    >>> model.heterozygote_geno(None)
    False
    '''
    return g is not None and None!=g[0]!=g[1]!=None

  def hemizygote_geno(self,g):
    '''
    Check if the genotype is hemizygote.

    @param g: genotype tuple
    @param g: tuple or None
    @return : True if g is a hemizygote, False otherwise
    @rtype  : bool

    >>> model = UnphasedPackedMarkerRepresentation('AG')
    >>> model.hemizygote_geno( ('A','A') )
    False
    >>> model.hemizygote_geno( ('A','G') )
    False
    >>> model.hemizygote_geno( ('A',None) )
    True
    >>> model.hemizygote_geno(None)
    False
    '''
    return bool(g is not None and (bool(g[0]) ^ bool(g[1])))

  def homozygote(self,rep):
    '''
    Check if the encoded genotype is homozygote.

    @param rep: encoded genotype
    @type  rep: int
    @return   : True if homozygote, False otherwise
    @rtype    : bool

    >>> g = snp_acgt2.rep_from_geno(('A','T'))
    >>> snp_acgt2.homozygote(g)
    False
    >>> g = snp_acgt2.rep_from_geno(('A','A'))
    >>> snp_acgt2.homozygote(g)
    True
    >>> g = snp_acgt2.rep_from_geno((None,'A'))
    >>> snp_acgt2.homozygote(g)
    False
    '''
    return rep.homozygote()

  def heterozygote(self,rep):
    '''
    Check if the encoded genotype is heterozygote.

    @param rep: the encoded genotype
    @type  rep: int
    @return    : True if heterozygote, False otherwise
    @rtype     : bool

    >>> g = snp_acgt2.rep_from_geno(('A','T'))
    >>> snp_acgt2.heterozygote(g)
    True
    >>> g = snp_acgt2.rep_from_geno(('A','A'))
    >>> snp_acgt2.heterozygote(g)
    False
    >>> g = snp_acgt2.rep_from_geno((None,'A'))
    >>> snp_acgt2.heterozygote(g)
    False
    '''
    return rep.heterozygote()

  def hemizygote(self,rep):
    '''
    Check if the encoded genotype is hemizygote.

    @param rep: the encoded genotype
    @type  rep: int
    @return    : True if hemizygote, False otherwise
    @rtype     : bool

    >>> g = snp_acgt2.rep_from_geno((None,'A'))
    >>> snp_acgt2.hemizygote(g)
    True
    >>> g = snp_acgt2.rep_from_geno(('A','A'))
    >>> snp_acgt2.hemizygote(g)
    False
    >>> g = snp_acgt2.rep_from_geno(('A','G'))
    >>> snp_acgt2.hemizygote(g)
    False
    '''
    return rep.hemizygote()

  def rep_from_geno(self,geno):
    '''
    Calculate encoded representation for the genotype.

    @param geno: the genotype
    @type  geno: tuple
    @return    : encoded genotype
    @rtype     : int

    >>> snp_acgt2.rep_from_geno( (None, 'A') )
    (None, 'A')
    >>> snp_acgt2.rep_from_geno( ('A', None) )
    (None, 'A')
    >>> snp_acgt2.rep_from_geno(None)
    (None, None)
    >>> snp_acgt2.rep_from_geno( (None,None) )
    (None, None)
    '''
    return self.geno_to_rep[geno]

  def reps_from_genos(self,genos):
    '''
    Calculate the encoded representation for the genotypes.

    @param genos: genotypes
    @type  genos: sequence
    @return     : encoded genotypes
    @rtype      : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt2.genos_from_reps( snp_acgt2.reps_from_genos(g) )
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    '''
    return list(imap(getitem, repeat(self.geno_to_rep), genos ))

  def _get_descr(self,n):
    if n in self._descr_cache:
      return self._descr_cache[n]
    descr = GenotypeArrayDescriptor([self.model]*n)
    self._descr_cache[n] = descr
    return descr

  def pack_genos(self,genos):
    '''
    Build a packed array (compact representation) for the genotypes.

    @param genos: genotypes
    @type  genos: iterable
    @return     : genotypes represented in a encoded array
    @rtype      : array

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> p = snp_acgt2.pack_genos(g)
    >>> p
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    >>> p.data
    array([ 86, 105], dtype=uint8)
    >>> snp_acgt2.strs_from_reps(p)
    ['AA', 'AC', 'AC', 'CC']
    '''
    try:
      n = len(genos)
    except ValueError:
      genos = list(genos)
      n = len(geons)

    descr = self._get_descr(n)
    return GenotypeArray(descr,imap(getitem, repeat(self.geno_to_rep), genos) )

  def rep_from_str(self,geno):
    '''
    Calculate encoded representation for the string genotype.

    @param geno: the genotype
    @type  geno: str
    @return    : encoded genotype
    @rtype     : int

    >>> snp_acgt2.rep_from_str(' A')
    (None, 'A')
    >>> snp_acgt2.rep_from_str('  ')
    (None, None)
    '''
    return self.str_to_rep[geno]

  def reps_from_strs(self,genos):
    '''
    Calculate the encoded representation for the string genotypes.

    @param genos: genotypes
    @type  genos: sequence of str
    @return     : encoded genotypes
    @rtype      : iterator

    >>> g = ['AA','AC','CA','CC']
    >>> snp_acgt2.strs_from_reps( snp_acgt2.reps_from_strs(g) )
    ['AA', 'AC', 'AC', 'CC']
    '''
    return list(imap(getitem, repeat(self.str_to_rep), genos))

  def pack_strs(self,genos):
    '''
    Build a packed array (compact representation) for the string genotypes.

    @param genos: sequence of str
    @type  genos: iterable
    @return     : encoded genotypes
    @rtype      : array

    >>> g = ['AA','AC','CA','CC']
    >>> p = snp_acgt2.pack_strs(g)
    >>> p
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    >>> p.data
    array([ 86, 105], dtype=uint8)
    '''
    try:
      n = len(genos)
    except ValueError:
      genos = list(genos)
      n = len(geons)

    descr = self._get_descr(n)
    return GenotypeArray(descr,imap(getitem, repeat(self.str_to_rep), genos) )

  def geno_from_rep(self,rep):
    '''
    Build the genotype tuple from an encoded genotype

    @param rep: encoded genotype
    @type  rep: object
    @return   : genotypes in tuple format
    @rtype    : tuple

    >>> snp_acgt2.geno_from_rep( snp_acgt2.rep_from_geno( ('A','A') ))
    ('A', 'A')
    >>> snp_acgt2.geno_from_rep( snp_acgt2.rep_from_geno( ('A','C') ))
    ('A', 'C')
    >>> snp_acgt2.geno_from_rep( snp_acgt2.rep_from_geno( (None,None) ))
    >>> snp_acgt2.geno_from_rep( snp_acgt2.rep_from_geno( None ))
    '''
    return self.rep_to_geno[rep]

  def str_from_rep(self,rep):
    '''
    Build the genotype string out of encoded genotype

    @param rep: genotypes in encoded format
    @type  rep: object
    @return   : genotypes in string format
    @rtype    : str

    >>> for g in ['AA','','AC']:
    ...   snp_acgt2.str_from_rep( snp_acgt2.rep_from_str(g) )
    'AA'
    '  '
    'AC'
    '''
    return self.rep_to_str[rep]

  def genos_from_reps(self,reps):
    '''
    Build the genotype tuples from encoded genotypes

    @param reps: genotypes in encoded format
    @type  reps: sequence
    @retur     : genotypes in tuple format
    @rtype     : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt2.genos_from_reps(snp_acgt2.pack_genos(g))
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    '''
    return list(imap(getitem, repeat(self.rep_to_geno), reps))

  def strs_from_reps(self,reps):
    '''
    Build the genotype strings out of encoded genotypes

    @param reps: genotypes in internal format
    @type  reps: sequence
    @return    : genotypes in string format
    @rtype     : iterator

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> snp_acgt2.strs_from_reps(snp_acgt2.pack_genos(g))
    ['AA', 'AC', 'AC', 'CC']
    '''
    return list(imap(getitem, repeat(self.rep_to_str), reps))

  def pack_reps(self,reps):
    '''
    Pack the reps into an array.

    @param reps: a sequence of encoded genotypes
    @type  reps: sequence
    @return    : a rep array (i.e. reps in compact form)
    @rtype     : array

    >>> g = map(tuple,['AA','AC','CA','CC'])
    >>> p = snp_acgt2.pack_reps(snp_acgt2.reps_from_genos(g))
    >>> p
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    >>> p.data
    array([ 86, 105], dtype=uint8)
    >>> list(p)
    [('A', 'A'), ('A', 'C'), ('A', 'C'), ('C', 'C')]
    '''
    try:
      n = len(reps)
    except ValueError:
      reps = list(reps)
      n = len(reps)

    descr = self._get_descr(n)
    return GenotypeArray(descr,reps)


snp_acgt       = UnphasedPackedMarkerRepresentation('ACGT', delimiter='', missing_geno_strs=['','  '],
                                                            missing_allele_str='',missing_allele_strs=['',' '])

snp_ab         = UnphasedPackedMarkerRepresentation('AB',   delimiter='', missing_geno_strs=['','  '],
                                                            missing_allele_str='',missing_allele_strs=['',' '])

snp_marker     = UnphasedMarkerRepresentation(delimiter='', missing_geno_str='  ',missing_geno_strs=['','  '],
                                                            missing_allele_str=' ',missing_allele_strs=['',' '])

generic_marker = UnphasedMarkerRepresentation(delimiter='/',missing_geno_str='',missing_geno_strs=['','/'],
                                                           missing_allele_str='',missing_allele_strs=['',' '])

snp_acgt2      = UnphasedGenotypeArrayRepresentation('ACGT', delimiter='', missing_geno_strs=['','  '],
                                                            missing_allele_str='',missing_allele_strs=['',' '])

snp_ab2        = UnphasedGenotypeArrayRepresentation('AB',   delimiter='', missing_geno_strs=['','  '],
                                                            missing_allele_str='',missing_allele_strs=['',' '])


def get_genorepr(reprname):
  '''
  Retrieve the supported genotype representation. Otherwise raises an ValueError exception

  @param reprname: internal genotype representation
  @type  reprname: str
  @return        : supported genotype representation
  @rtype         : genotype representation object
  '''

  reprs = { 'snp_marker' : snp_marker,
            'snp_acgt'   : snp_acgt,
            'snp_ab'     : snp_ab,
            'snp_acgt2'  : snp_acgt2,
            'snp_ab2'    : snp_ab2,
            'generic'    : generic_marker }

  if not reprname:
    reprname = 'snp_marker'

  try:
    return reprs[reprname]
  except KeyError:
    raise ValueError, 'Unknown genotype representation: %s' % reprname


class genoarrayTest(unittest.TestCase):
  def test_byte(self):
    self.assertRaises(KeyError, snp_acgt.rep_from_geno, ())
    self.assertEquals(snp_acgt.rep_from_geno(None), 0)

  def test_byte_str(self):
    self.assertEquals(list(snp_acgt.reps_from_strs(['AA','CT','GC','','  '])),[17, 36, 35, 0, 0])


def main():
  if 0:
    models = [snp_acgt,snp_ab,UnphasedPackedMarkerRepresentation(map(str,range(300,317)))]

    for snp in models:
      print 'alleles:',snp.alleles
      print 'genostr-to-rep dict:',snp.str_to_rep
      print 'geno-to-rep dict:',snp.geno_to_rep
      print 'None mapped:',None in snp.geno_to_rep
      print 'Empty 2 tuple mapped:',('','') in snp.geno_to_rep
      print 'Missing geno:',snp.rep_to_geno[0]
      print 'rep-to-geno list:',[ (i,g) for i,g in enumerate(snp.rep_to_geno) if g ]
      print 'rep-to-geno, rep_class:',[ (b,c) for b,c in zip(snp.rep_to_geno,snp.rep_class) if b ]

  import doctest
  doctest.testmod()
  unittest.main()


if __name__ == '__main__':
  main()
