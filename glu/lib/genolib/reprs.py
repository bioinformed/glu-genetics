# -*- coding: utf-8 -*-
'''
File:          reprs.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Genotype representation parsers for text file input/output

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from   operator   import getitem
from   itertools  import imap, repeat

from   genoarray  import model_from_alleles


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

  def to_string(self, rep):
    '''
    Transform a genotype in internal format to a genotype string

    @param rep: genotype in internal format
    @type  rep: object
    @return   : genotype in string format
    @rtype    : str

    >>> model=model_from_alleles('AG',allow_hemizygote=True)
    >>> snp.to_string( model['A','A'] )
    'AA'
    >>> snp.to_string( model[None, 'A'] )
    ' A'
    >>> snp.to_string( model['A', None] )
    ' A'
    >>> snp.to_string( model[None,None])
    '  '
    >>> hapmap.to_string( model['A','A'] )
    'AA'
    >>> hapmap.to_string( model[None, 'A'] )
    'NA'
    >>> hapmap.to_string( model['A', None] )
    'NA'
    >>> hapmap.to_string( model[None,None])
    'NN'
    >>> marker.to_string( model['A','A'] )
    'A/A'
    >>> marker.to_string( model[None, 'A'] )
    '/A'
    >>> marker.to_string( model['A', None] )
    '/A'
    >>> marker.to_string( model[None,None])
    ''
    '''
    rep = rep.alleles()

    strcache = self.strcache
    if rep in strcache:
      return strcache[rep]

    if not rep or rep == self.missingrep:
      if self.missing_geno_str is not Nothing:
        strcache[rep] = self.missing_geno_str
        return self.missing_geno_str
      trep = self.missingrep
    else:
      trep = rep

    trep = self.delimiter.join(self.allelerep.get(a,a) for a in trep)
    strcache[rep] = trep
    return trep

  def to_strings(self,reps):
    '''
    Transform a sequence of genotypes in internal format to a genotype strings

    @param rep: genotypes in internal format
    @type  rep: sequence
    @return   : genotypes in string format
    @rtype    : sequence

    >>> model=model_from_alleles('AG',allow_hemizygote=True)
    >>> genos = [ model[g] for g in [('A','A'),(None,'A'),(None,None)] ]
    >>> snp.to_strings(genos)
    ['AA', ' A', '  ']
    >>> hapmap.to_strings(genos)
    ['AA', 'NA', 'NN']
    >>> marker.to_strings(genos)
    ['A/A', '/A', '']
    '''
    try:
      strcache = self.strcache
      return [ strcache[r.alleles()] for r in reps ]
    except KeyError:
      return map(self.to_string, reps)

  def from_string(self, geno):
    '''
    Transform a genotype string to internal format

    @param geno: genotype string
    @type  geno: str
    @return    : genotype in internal format
    @rtype     : object

    >>> snp.from_string('AA')
    ('A', 'A')
    >>> snp.from_string('AB')
    ('A', 'B')
    >>> snp.from_string('BA')
    ('A', 'B')
    >>> snp.from_string(' A')
    (None, 'A')
    >>> snp.from_string('A ')
    (None, 'A')
    >>> snp.from_string('  ')
    (None, None)
    >>> snp.from_string('')
    (None, None)

    >>> hapmap.from_string('AA')
    ('A', 'A')
    >>> hapmap.from_string('AB')
    ('A', 'B')
    >>> hapmap.from_string('BA')
    ('A', 'B')
    >>> hapmap.from_string('NA')
    (None, 'A')
    >>> hapmap.from_string('AN')
    (None, 'A')
    >>> hapmap.from_string('NN')
    (None, None)

    >>> marker.from_string('A/A')
    ('A', 'A')
    >>> marker.from_string('A/B')
    ('A', 'B')
    >>> marker.from_string('B/A')
    ('A', 'B')
    >>> marker.from_string(' /A')
    (None, 'A')
    >>> marker.from_string('/A')
    (None, 'A')
    >>> marker.from_string('A/ ')
    (None, 'A')
    >>> marker.from_string('A/')
    (None, 'A')
    >>> marker.from_string(' / ')
    (None, None)
    >>> marker.from_string('/')
    (None, None)
    >>> marker.from_string('/ ')
    (None, None)
    >>> marker.from_string(' /')
    (None, None)
    >>> marker.from_string('')
    (None, None)
    '''
    strcache = self.strcache
    if geno in strcache:
      return strcache[geno]

    if geno in self.missing_geno_strs:
      missingrep = self.missingrep
      self.strcache[geno] = missingrep
      return missingrep

    if self.delimiter:
      alleles = geno.split(self.delimiter)
    else:
      alleles = geno

    if len(alleles) != 2:
      raise ValueError('Invalid genotype: %s' % geno)

    rep = tuple(sorted(self.allelestr.setdefault(a,a) for a in alleles))

    if not rep or rep == self.missingrep:
      missingrep = self.missingrep
      self.strcache[geno] = missingrep
      return missingrep

    self.strcache[intern(geno)] = rep
    return rep

  def from_strings(self,genos):
    '''
    Transform a sequence of genotype strings to internal format

    @param geno: sequence of genotype strings
    @type  geno: sequence of str
    @return    : sequence of genotype in internal format
    @rtype     : sequence of objects

    >>> snp.from_strings(['AA','AB','BA',' A','A ','  ',''])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), (None, None), (None, None)]
    >>> hapmap.from_strings(['AA','AB','BA','NA','AN','NN'])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), (None, None)]
    >>> marker.from_strings(['A/A','A/B','B/A',' /A','A/ ',' / ',''])
    [('A', 'A'), ('A', 'B'), ('A', 'B'), (None, 'A'), (None, 'A'), (None, None), (None, None)]
    '''
    try:
      return list(imap(getitem, repeat(self.strcache), genos))
    except KeyError:
      return map(self.from_string, genos)


snp    = UnphasedMarkerRepresentation(delimiter='', missing_geno_str='  ',missing_geno_strs=['','  '],
                                      missing_allele_str=' ',missing_allele_strs=['',' '])
hapmap = UnphasedMarkerRepresentation(delimiter='',missing_geno_str='NN',missing_geno_strs=['','NN'],
                                      missing_allele_str='N',missing_allele_strs=['N'],)
marker = UnphasedMarkerRepresentation(delimiter='/',missing_geno_str='',missing_geno_strs=['','/'],
                                      missing_allele_str='',missing_allele_strs=['',' '])

def get_genorepr(reprname):
  '''
  Retrieve the supported genotype representation. Otherwise raises an ValueError exception

  @param reprname: internal genotype representation
  @type  reprname: str
  @return        : supported genotype representation
  @rtype         : genotype representation object
  '''

  reprs = { 'snp'    : snp,
            'hapmap' : hapmap,
            'marker' : marker }

  if not reprname:
    reprname = 'snp'

  try:
    return reprs[reprname]
  except KeyError:
    raise ValueError, 'Unknown genotype representation: %s' % reprname


def main():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  main()
