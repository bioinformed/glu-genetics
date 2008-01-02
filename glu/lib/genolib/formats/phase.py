# -*- coding: utf-8 -*-
'''
File:          phase.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU Phase genotype format output

Requires:      Python 2.5

Revision:      $Id$
'''

from __future__ import with_statement

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from   itertools                 import islice

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,compressed_filename,  \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.locus     import Genome


__all__ = ['PhaseWriter', 'save_phase']


class PhaseWriter(object):
  '''
  Object to write Phase data

  See http://stephenslab.uchicago.edu/instruct2.1.pdf

  Caveats:
    1) Uses back-patching to set the number of samples written.  This means
       that it can only write to real files.  Should buffer if this is not
       possible.
    2) Assumes all loci are SNPs, but should probably buffer if this is in
       doubt.  Since the missing genotype encoding is different for non-SNP
       data, this cannot be easily corrected via backpatching.
    3) Does not sort by location, nor warn if trying to write more than one
       chromosome.
    4) Does not check for commas in identifiers, since they are declared
       illegal.

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with PhaseWriter(o,genos.loci,genos.genome) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  00000003
  00000003
  SSS
  s1
  A ? C
  A ? T
  s2
  A C C
  G G C
  s3
  G ? C
  G ? T
  '''
  def __init__(self,filename,loci,genome,extra_args=None,**kwargs):
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
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if compressed_filename(filename):
      raise ValueError('PHASE genotype files must not have a compressed extension')

    self.out         = autofile(filename,'w')
    self.loci        = loci
    self.genome      = genome
    self.samplecount = 0

    locations = [ genome.get_locus(locus).location for locus in loci ]

    self.out.write('00000000\n')
    self.out.write('%08d\n' % len(loci))

    if None not in locations:
      self.out.write('P ')
      self.out.write(' '.join(map(str,locations)))
      self.out.write('\n')

    self.out.write( 'S'*len(loci) )
    self.out.write('\n')

  def writerow(self, sample, genos):
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

    row1 = [ g.allele1 or '?' for g in genos ]
    row2 = [ g.allele2 or '?' for g in genos ]

    out.write('%s\n' % sample)
    out.write(' '.join(row1))
    out.write('\n')
    out.write(' '.join(row2))
    out.write('\n')
    self.samplecount += 1

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

    n = len(self.loci)

    for sample,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')

      row1 = [ g.allele1 or '?' for g in genos ]
      row2 = [ g.allele2 or '?' for g in genos ]

      out.write('%s\n' % sample)
      out.write(' '.join(row1))
      out.write('\n')
      out.write(' '.join(row2))
      out.write('\n')
      self.samplecount += 1

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
    out,self.out = self.out,None

    out.seek(0)
    out.write('%08d' % self.samplecount)

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


def save_phase(filename,genos,extra_args=None,**kwargs):
  '''
  Write the genotype matrix data to PHASE format.

  See http://stephenslab.uchicago.edu/instruct2.1.pdf

  Caveats:
    1) Uses back-patching to set the number of samples written.  This means
       that it can only write to real files.  Should buffer if this is not
       possible.
    2) Assumes all loci are SNPs, but should probably buffer if this is in
       doubt.  Since the missing genotype encoding is different for non-SNP
       data, this cannot be easily corrected via backpatching.
    3) Does not sort by location, nor warn if trying to write more than one
       chromosome.
    4) Does not check for commas in identifiers, since they are declared
       illegal.

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
  >>> save_phase(o,genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  00000003
  00000003
  SSS
  s1
  A ? C
  A ? T
  s2
  A C C
  G G C
  s3
  G ? C
  G ? T
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_sdat(mergefunc)

  with PhaseWriter(filename, genos.loci, genos.genome, extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
