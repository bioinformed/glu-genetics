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

from   itertools                 import islice,izip,tee

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import model_from_alleles
from   glu.lib.genolib.reprs     import snp
from   glu.lib.genolib.locus     import Genome


__all__ = ['PlinkWriter', 'PlinkWriter',
           'save_plink', 'load_plink_ped' ]


ALLELE_MAP = {None:'0'}


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


def load_plink_ped(filename,genome=None,unique=True,extra_args=None,**kwargs):
  '''
  Load a Plink format genotype data file.

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
    amap  = {'0':None,'A':'A','C':'C','G':'G','T':'T','1':'1','2':'2','3':'3','4':'4'}

    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled:
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        sample = '%s_%s' % (fields[0],fields[1])

        fields = [ amap.get(a,a) for a in islice(fields,6,None) ]
        genos  = zip(islice(fields,0,None,2),islice(fields,1,None,2))

      yield sample,genos

  return GenomatrixStream.from_tuples(_load_plink(),'sdat',loci=loci,genome=genome,unique=unique)


class PlinkWriter(object):
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
  >>> with PlinkWriter(o,genos.loci) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  sdat  l1      l2      l3
  s1    AA              CT
  s2    AG      CG      CC
  s3    GG              CT
  '''
  def __init__(self,filename,loci,mapfile,genome,dialect='tsv'):
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
    @param      dialect: csv module dialect name ('csv' or 'tsv', default is 'tsv')
    @type       dialect: str or csv dialect object
    '''
    if genorepr is None:
      raise ValueError('genotype representation must be specified when reading a text genotype format')

    self.out       = csv.writer(autofile(filename,'w'),dialect=dialect)
    self.loci      = loci

    self.write_map(mapfile,loci,genome)


    self.out.writerow( [format]+[h.strip() for h in self.header] )

  def writerow(self, sample, genos):
    '''
    Write a row of genotypes given the row key and list of genotypes

    @param rowkey: row identifier
    @type  rowkey: str
    @param  genos: sequence of genotypes in an internal representation, to
                   be converted to the appropiate string representation by
                   the supplied genorepr class.
    @type   genos: sequence
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    if len(genos) != len(self.loci):
      raise ValueError('[ERROR] Internal error: Genotypes do not match header')

    row = [sample,sample,'0','0','?','0']
    for g in samplegenos:
      row += [ ALLELE_MAP[a] for a in g ]
    out.write('%s\n' % ' '.join(row))

  def writerows(self, rows):
    '''
    Write rows of genotypes given pairs of row key and list of genotypes

    @param rows: sequence of pairs of row key and sequence of genotypes in
                 an internal representation, to be converted to the
                 appropiate string representation by the supplied genorepr
                 class.
    @type  rows: sequence of (str,sequence)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    n = len(self.loci)

    for rowkey,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')
      row = [sample,sample,'0','0','?','0']
      for g in genos:
        row += [ ALLELE_MAP[a] for a in g ]
      out.write('%s\n' % ' '.join(row))

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')
    self.out = None

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


def save_plink(filename,genos,genome):
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
  >>> save_plink(o,genos,snp)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  sdat	l1	l2	l3
  s1	AA	  	CT
  s2	AG	CG	CC
  s3	GG	  	CT
  '''
  genos = genos.as_sdat()
  with PlinkWriter(filename, genos.loci) as writer:
    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
