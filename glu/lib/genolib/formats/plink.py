# -*- coding: utf-8 -*-
'''
File:          plink.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU Plink genotype format input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

from __future__ import with_statement

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import re
import csv

from   itertools                 import islice

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN, \
                                        PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


__all__ = ['PlinkPedWriter',  'PlinkTPedWriter',
           'save_plink_ped',  'load_plink_ped',
           'save_plink_tped', 'load_plink_tped' ]


ALLELE_MAP  = {'0':None}
ALLELE_RMAP = {None:'0'}

SEX_MAP    = {'1':SEX_MALE,'2':SEX_FEMALE}
SEX_RMAP   = {SEX_UNKNOWN:'0', SEX_MALE:'1', SEX_FEMALE:'2'}

PHENO_MAP  = {'1':PHENO_UNAFFECTED, '2':PHENO_AFFECTED}
PHENO_RMAP = {PHENO_UNKNOWN:'0',PHENO_UNAFFECTED:'1',PHENO_AFFECTED:'2'}

PARENT_MAP = {'0':None}


def load_plink_map(filename,genome):
  mfile = autofile(filename)

  for i,line in enumerate(mfile):
    line = line.rstrip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    if len(fields) != 4:
      raise ValueError('Invalid PLINK MAP record %d' % (i+1))

    chr   = fields[0] or None
    lname = fields[1]
    gdist = int(fields[2])      if fields[2] else None
    pdist = abs(int(fields[3])) if fields[3] else None

    if not lname:
      raise ValueError('Invalid PLINK MAP record %d' % (i+1))

    genome.merge_locus(lname, chromosome=chr, location=pdist)

    yield lname


def load_plink_ped(filename,genome=None,phenome=None,unique=True,extra_args=None,**kwargs):
  '''
  Load a Plink PED format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  loci = get_arg(args, ['loci'])
  lmap = get_arg(args, ['map'])

  if loci is None and lmap is None:
    raise ValueError('Map file or locus list must be specified when loading PLINK PED files')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if genome is None:
    genome = Genome()

  if phenome is None:
    phenome = Phenome()

  if loci is not None:
    if isinstance(loci,basestring):
      loci = list(load_locus_records(loci)[2])
      # Merge map data into genome
      populate_genome(genome,loci)
      loci  = [ intern(l[0]) for l in loci ]

  if lmap is not None:
    map_loci = list(load_plink_map(lmap,genome))
    if loci is not None:
      if loci != map_loci:
        raise ValueError('Locus list and PLINK MAP file are not identical')
    else:
      loci = map_loci

  gfile = autofile(filename)
  n     = 6 + 2*len(loci)

  def _load_plink():
    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled:
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        family,name,father,mother,sex,pheno = [ s.strip() for s in fields[:6] ]

        if name == '0':
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        father  = PARENT_MAP.get(father,father)
        mother  = PARENT_MAP.get(mother,mother)

        ename   = '%s:%s' % (family,name)
        efather = '%s:%s' % (family,father) if father else None
        emother = '%s:%s' % (family,mother) if mother else None
        sex     = SEX_MAP.get(sex,SEX_UNKNOWN)
        pheno   = PHENO_MAP.get(pheno,PHENO_UNKNOWN)

        if father:
          phenome.merge_phenos(efather, family, father, sex=SEX_MALE)
        if mother:
          phenome.merge_phenos(emother, family, mother, sex=SEX_FEMALE)

        phenome.merge_phenos(ename, family, name, efather, emother, sex, pheno)

        fields = [ ALLELE_MAP.get(a,a) for a in islice(fields,6,None) ]
        genos  = zip(islice(fields,0,None,2),islice(fields,1,None,2))

      yield ename,genos

  return GenomatrixStream.from_tuples(_load_plink(),'sdat',loci=loci,genome=genome,phenome=phenome,unique=unique)


class PlinkPedWriter(object):
  '''
  Object to write the genotype matrix data to a text file

  Genotype matrix files are delimited ASCII files with the following format:

  format	heading1	heading2	heading2	...
  rowkey1	G11		G12		G13		...
  rowkey2	G21		G22		G23		...
  rowkey3	G31		G32		G33		...
  ...

  All rows must have the same number of columns, as determined by the header
  supplied, with each subsequent data row conforming.  Headers and row keys
  are user-specified, although they are typically either sample/subject
  identifiers and locus identifiers, although which are mapped to rows and
  columns is arbitrary.  When loci are given in rows, the format is
  typically 'ldat' and when samples are given in rows the format is
  typically 'sdat'.  However, these formatting issues are handled at a
  higher level by callers of TextGenomatrixWriter.

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> m = StringIO()
  >>> with PlinkPedWriter(o,genos.loci,genos.genome,genos.phenome,m) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 ? 0 A A 0 0 C T
  s2 s2 0 0 ? 0 A G C G C C
  s3 s3 0 0 ? 0 G G 0 0 C T
  >>> print m.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0
  0 l2 0 0
  0 l3 0 0
  '''
  def __init__(self,filename,loci,genome,phenome,mapfile=None):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param       header: column headings
    @type        header: list or str
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    '''
    self.out       = autofile(filename,'wb')
    self.loci      = loci
    self.genome    = genome
    self.phenome   = phenome
    self.mapfile   = mapfile

  def writerow(self, sample, genos, phenome=None):
    '''
    Write a row of genotypes given the row key and list of genotypes

    @param rowkey: row identifier
    @type  rowkey: str
    @param  genos: sequence of genotypes in an internal representation
    @type   genos: sequence
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    if len(genos) != len(self.loci):
      raise ValueError('[ERROR] Internal error: Genotypes do not match header')

    if phenome is None:
      phenome = self.phenome
    if phenome is None:
      phenome = Phenome()

    phenos     = phenome.get_phenos(sample)
    family     = phenos.family
    individual = phenos.individual or sample
    parent1    = phenos.parent1
    parent2    = phenos.parent2

    if parent1 and parent2:
      p1 = phenome.get_phenos(parent1)
      p2 = phenome.get_phenos(parent2)
      if p1.sex is SEX_FEMALE or p2.sex is SEX_MALE:
        parent1,parent2 = parent2,parent1

    sex   = SEX_RMAP[phenos.sex]
    pheno = PHENO_RMAP[phenos.phenoclass]

    row = [family or individual,individual,parent1 or '0',parent2 or '0',sex,pheno]
    for g in genos:
      row += [ ALLELE_RMAP.get(a,a) for a in g ]
    out.write(' '.join(row))
    out.write('\r\n')

  def writerows(self, rows, phenome=None):
    '''
    Write rows of genotypes given pairs of row key and list of genotypes

    @param rows: sequence of pairs of row key and sequence of genotypes in
                 an internal representation
    @type  rows: sequence of (str,sequence)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    n = len(self.loci)

    if phenome is None:
      phenome = self.phenome
    if phenome is None:
      phenome = Phenome()

    for sample,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')

      phenos     = phenome.get_phenos(sample)
      family     = phenos.family
      individual = phenos.individual or sample
      parent1    = phenos.parent1
      parent2    = phenos.parent2

      if parent1 and parent2:
        p1 = phenome.get_phenos(parent1)
        p2 = phenome.get_phenos(parent2)
        if p1.sex is SEX_FEMALE or p2.sex is SEX_MALE:
          parent1,parent2 = parent2,parent1

      sex   = SEX_RMAP[phenos.sex]
      pheno = PHENO_RMAP[phenos.phenoclass]

      row = [family or individual,individual,parent1 or '0',parent2 or '0',sex,pheno]

      for g in genos:
        row += [ ALLELE_RMAP.get(a,a) for a in g ]

      out.write(' '.join(row))
      out.write('\r\n')

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')

    # FIXME: Closing out causes problems with StringIO objects used for
    #        testing
    #self.out.close()
    self.out = None

    if self.mapfile:
      out = autofile(self.mapfile,'wb')
      for locus in self.loci:
        loc = self.genome.get_locus(locus)
        out.write( ' '.join( map(str,[loc.chromosome or 0,locus,0,loc.location or 0] )) )
        out.write('\r\n')

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()


def save_plink_ped(filename,genos,mapfile=None):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...           ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...           ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> save_plink_ped(o,genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 ? 0 A A 0 0 C T
  s2 s2 0 0 ? 0 A G C G C C
  s3 s3 0 0 ? 0 G G 0 0 C T
  '''
  genos = genos.as_sdat()
  with PlinkPedWriter(filename, genos.loci, genos.genome, genos.phenome, mapfile) as writer:
    writer.writerows(genos)


###############################################################################################


def load_plink_tfam(filename,phenome):
  mfile = autofile(filename)

  for i,line in enumerate(mfile):
    line = line.rstrip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    if len(fields) != 6:
      raise ValueError('Invalid PLINK TFAM record %d' % (i+1))

    family,name,father,mother,sex,pheno = [ s.strip() for s in fields ]

    if name == '0':
      raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

    father  = PARENT_MAP.get(father,father)
    mother  = PARENT_MAP.get(mother,mother)

    ename   = '%s:%s' % (family,name)
    efather = '%s:%s' % (family,father) if father else None
    emother = '%s:%s' % (family,mother) if mother else None
    sex     = SEX_MAP.get(sex,SEX_UNKNOWN)
    pheno   = PHENO_MAP.get(pheno,PHENO_UNKNOWN)

    if father:
      phenome.merge_phenos(efather, family, father, sex=SEX_MALE)
    if mother:
      phenome.merge_phenos(emother, family, mother, sex=SEX_FEMALE)

    phenome.merge_phenos(ename, family, name, efather, emother, sex, pheno)

    yield ename


def load_plink_tped(filename,genome=None,phenome=None,unique=True,extra_args=None,**kwargs):
  '''
  Load a Plink TPED format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  loci = get_arg(args, ['loci'])
  tfam = get_arg(args, ['tfam'])

  if tfam is None:
    raise ValueError('A TFAM file must be specified when loading PLINK TPED data')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if genome is None:
    genome = Genome()

  if phenome is None:
    phenome = Phenome()

  if loci is not None:
    if isinstance(loci,basestring):
      loci = list(load_locus_records(loci)[2])
      # Merge map data into genome
      populate_genome(genome,loci)
      loci  = [ intern(l[0]) for l in loci ]

  samples = list(load_plink_tfam(tfam,phenome))

  gfile = autofile(filename)
  n     = 4 + 2*len(samples)

  def _load_plink():
    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled:
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        chr   = fields[0] or None
        lname = fields[1]
        gdist = int(fields[2])      if fields[2] else None
        pdist = abs(int(fields[3])) if fields[3] else None

        if not lname:
          raise ValueError('Invalid PLINK TPED record %d' % (i+1))

        genome.merge_locus(lname, chromosome=chr, location=pdist)

        fields = [ ALLELE_MAP.get(a,a) for a in islice(fields,4,None) ]
        genos  = zip(islice(fields,0,None,2),islice(fields,1,None,2))

      yield lname,genos

  return GenomatrixStream.from_tuples(_load_plink(),'ldat',samples=samples,genome=genome,phenome=phenome,unique=unique)


class PlinkTPedWriter(object):
  '''
  Object to write a PLINK TPED file

  See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_ldat()
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> m = StringIO()
  >>> with PlinkTPedWriter(o,genos.samples,genos.genome,genos.phenome,m) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0 A A A G G G
  0 l2 0 0 0 0 C G 0 0
  0 l3 0 0 C T C C C T
  >>> print m.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 ? 0
  s2 s2 0 0 ? 0
  s3 s3 0 0 ? 0
  '''
  def __init__(self,filename,samples,genome,phenome,tfamfile=None):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param       header: column headings
    @type        header: list or str
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    '''
    self.out       = autofile(filename,'wb')
    self.samples   = samples
    self.genome    = genome
    self.phenome   = phenome
    self.tfamfile  = tfamfile

  def writerow(self, locus, genos):
    '''
    Write a row of genotypes given the row key and list of genotypes

    @param rowkey: row identifier
    @type  rowkey: str
    @param  genos: sequence of genotypes in an internal representation
    @type   genos: sequence
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    if len(genos) != len(self.samples):
      raise ValueError('[ERROR] Internal error: Genotypes do not match header')

    loc = self.genome.get_locus(locus)

    row = map(str,[loc.chromosome or 0,locus,0,loc.location or 0] )

    for g in genos:
      row += [ ALLELE_RMAP.get(a,a) for a in g ]

    out.write(' '.join(row))
    out.write('\r\n')

  def writerows(self, rows):
    '''
    Write rows of genotypes given pairs of row key and list of genotypes

    @param rows: sequence of pairs of row key and sequence of genotypes in
                 an internal representation
    @type  rows: sequence of (str,sequence)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    n = len(self.samples)

    for locus,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')

      loc = self.genome.get_locus(locus)

      row = map(str,[loc.chromosome or 0,locus,0,loc.location or 0] )

      for g in genos:
        row += [ ALLELE_RMAP.get(a,a) for a in g ]

      out.write(' '.join(row))
      out.write('\r\n')

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')

    # FIXME: Closing out causes problems with StringIO objects used for
    #        testing
    #self.out.close()
    self.out = None

    if self.tfamfile:
      out = autofile(self.tfamfile,'wb')
      for sample in self.samples:
        phenos     = self.phenome.get_phenos(sample)
        family     = phenos.family
        individual = phenos.individual or sample
        parent1    = phenos.parent1
        parent2    = phenos.parent2

        if parent1 and parent2:
          p1 = phenome.get_phenos(parent1)
          p2 = phenome.get_phenos(parent2)
          if p1.sex is SEX_FEMALE or p2.sex is SEX_MALE:
            parent1,parent2 = parent2,parent1

        sex   = SEX_RMAP[phenos.sex]
        pheno = PHENO_RMAP[phenos.phenoclass]

        row = [family or individual,individual,parent1 or '0',parent2 or '0',sex,pheno]
        out.write( ' '.join(row))
        out.write('\r\n')

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()


def save_plink_tped(filename,genos,tfamfile=None):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...           ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...           ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> save_plink_tped(o,genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0 A A A G G G
  0 l2 0 0 0 0 C G 0 0
  0 l3 0 0 C T C C C T
  '''
  genos = genos.as_ldat()
  with PlinkTPedWriter(filename, genos.samples, genos.genome, genos.phenome, tfamfile) as writer:
    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
