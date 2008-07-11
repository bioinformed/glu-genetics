# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'WTCCC genotype format output object'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   glu.lib.fileutils         import autofile,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream


__all__ = ['WTCCCWriter', 'save_wtccc']


GENOS = [ ['0','0','0'],
          ['1','0','0'],
          ['0','1','0'],
          ['0','0','1'] ]


def _encode_wtccc(model):
  assert model is not None

  if len(model.alleles) > 3:
    raise ValueError('WTCCC files support only biallic models')

  genovalues = [ GENOS[0] ]*4

  allele1,allele2 = (model.alleles[1:]+[None,None])[:2]

  if allele1:
    genovalues[ model[allele1,allele1].index ] = GENOS[1]
  if allele2:
    genovalues[ model[allele2,allele2].index ] = GENOS[3]
  if allele1 and allele2:
    genovalues[ model[allele1,allele2].index ] = GENOS[2]

  return allele1,allele2,genovalues


class WTCCCWriter(object):
  '''
  Object to write WTCCC data

  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_ldat()
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with WTCCCWriter(o,genos.samples,genos.genome) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 A G 1 0 0 0 1 0 0 0 1
  0 l2 0 C G 0 0 0 0 1 0 0 0 0
  0 l3 0 C T 0 1 0 1 0 0 0 1 0
  '''
  def __init__(self,filename,samples,genome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param      samples: column headings
    @type       samples: list of str
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    self.out       = autofile(filename,'wb')
    self.samples   = samples
    self.genome    = genome

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
    allele1,allele2,genovalues = _encode_wtccc(loc.model)

    row = [ loc.chromosome or '0', locus, str(loc.location or '0'), allele1 or '?', allele2 or '?' ]

    for g in genos:
      row += genovalues[g.index]

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
    genome = self.genome

    for locus,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')

      loc = genome.get_locus(locus)
      allele1,allele2,genovalues = _encode_wtccc(loc.model)

      row = [ loc.chromosome or '0', locus, str(loc.location or '0'), allele1 or '?', allele2 or '?' ]

      for g in genos:
        row += genovalues[g.index]

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


def save_wtccc(filename,genos,extra_args=None,**kwargs):
  '''
  Write genotype data to a WTCCC file

  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html

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
  >>> save_wtccc(o,genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 A G 1 0 0 0 1 0 0 0 1
  0 l2 0 C G 0 0 0 0 1 0 0 0 0
  0 l3 0 C T 0 1 0 1 0 0 0 1 0
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_ldat(mergefunc)

  with WTCCCWriter(filename, genos.samples, genos.genome, extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
