# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'FUD genotype format output object'
__copyright__ = 'Copyright (c) 2007-2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['FUDWriter', 'save_fud']

__genoformats__ = [
  #LOADER        SAVER          WRITER       PFORMAT  ALIAS    EXTS
  ('load_fud',   'save_fud',    'FUDWriter', 'ldat',  None,    'fud') ]


import csv

from   glu.lib.fileutils         import autofile,parse_augmented_filename,get_arg,trybool,namefile

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.phenos    import Phenome
from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.genoarray import GenotypeArray,build_model,build_descr


HEADER_START = ['SNP','Chromosome','Position','AlleleA','AlleleB']


def _encode_fud(model):
  assert model is not None

  if len(model.alleles) > 3:
    raise ValueError('WTCCC files support only biallelic models')

  genovalues = [-1]*4

  allele1,allele2 = (model.alleles[1:]+[None,None])[:2]

  if allele1:
    genovalues[ model[allele1,allele1].index ] = 0
  if allele2:
    genovalues[ model[allele2,allele2].index ] = 2
  if allele1 and allele2:
    genovalues[ model[allele1,allele2].index ] = 1

  return allele1,allele2,genovalues


def _decode_fud(model,a,b,fields):
  remap = [model[None,None]]*4

  remap[0]  = model[a,a]
  remap[1]  = model[a,b]
  remap[2]  = model[b,b]

  return ( remap[int(round(float(g)))] for g in fields )


def load_fud(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  FUD format

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO(
  ... 'SNP\\tChromosome\\tPosition\\tAlleleA\\tAlleleB\\ts1\\ts2\\ts3\\n'
  ... 'l1\\t0\\tl1\\tA\\tG\\t0\\t1\\t2\\n'
  ... 'l2\\t0\\tl2\\tC\\tG\\t-1\\t1\\t-1\\n'
  ... 'l3\\t0\\tl3\\tC\\tT\\t1\\t0\\t1\\n')
  >>> genos = load_fud(data,'fud')
  >>> genos.format
  'ldat'
  >>> genos.columns
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [(None, None), ('C', 'G'), (None, None)])
  ('l3', [('C', 'T'), ('C', 'C'), ('C', 'T')])
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  unique    = trybool(get_arg(args, ['unique'], True))

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  gfile   = csv.reader(autofile(filename),dialect='excel-tab')
  header  = gfile.next()

  assert header[:5] == HEADER_START
  samples = header[5:]

  if phenome is None:
    phenome = Phenome()

  loci   = []
  models = []
  file_genome = Genome()

  def _load_fud():
    m           = len(samples)
    n           = 5+m
    modelcache  = {}

    for line_num,row in enumerate(gfile):
      if len(row) != n:
        raise ValueError('Invalid FUD row on line %d of %s' % (line_num+1,namefile(filename)))

      lname   = row[0]
      chrom   = row[1]
      loc     = row[2]
      a,b     = row[3:5]
      alleles = a,b

      if chrom.startswith('chr'):
        chrom = chrom[3:]

      if chrom.upper()=='MT':
        chrom = 'M'

      try:
        loc = int(loc) or None
      except ValueError:
        loc = None

      model = modelcache.get(alleles)
      if model is None:
        model = modelcache[alleles] = build_model(genotypes=[(a,a),(a,b),(b,b)],max_alleles=2)

      descr = build_descr(model,m)

      loci.append(lname)
      models.append(model)
      file_genome.set_locus(lname,model,chrom,loc)

      genos = GenotypeArray(descr, _decode_fud(model,a,b,row[5:]))

      yield lname,genos

  genos = GenomatrixStream(_load_fud(),'ldat',samples=samples,loci=loci,models=models,
                                              genome=file_genome,phenome=phenome,
                                              unique=unique,packed=True)

  if unique:
    genos = genos.unique_checked()

  if genome:
    genos = genos.transformed(recode_models=genome)

  return genos


class FUDWriter(object):
  '''
  Object to write FUD data

  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_ldat()
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with FUDWriter(o,'fud',genos.samples,genos.genome,genos.phenome) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  SNP Chromosome Position AlleleA AlleleB s1   s2   s3
  l1                         A       G     0    1    2
  l2                         C       G    -1    1   -1
  l3                         C       T     1    0    1
  '''
  def __init__(self,filename,format,samples,genome,phenome,extra_args=None,**kwargs):
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

    self.out       = csv.writer(autofile(filename,'wb'),dialect='excel-tab')
    self.samples   = samples
    self.genome    = genome

    self.out.writerow( HEADER_START+[s.strip() for s in samples] )


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
    allele1,allele2,genovalues = _encode_fud(loc.model)

    row  = [ locus, loc.chromosome or '', str(loc.location or ''), allele1 or '', allele2 or '' ]
    row += (genovalues[g.index] for g in genos)

    out.writerow(row)

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
      allele1,allele2,genovalues = _encode_fud(loc.model)

      row  = [ locus, loc.chromosome or '', str(loc.location or ''), allele1 or '', allele2 or '' ]
      row += (genovalues[g.index] for g in genos)

      out.writerow(row)

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


def save_fud(filename,genos,format,extra_args=None,**kwargs):
  '''
  Write genotype data to a FUD file

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
  >>> save_fud(o,genos,'fud')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  SNP Chromosome Position AlleleA AlleleB s1   s2   s3
  l1                         A       G     0    1    2
  l2                         C       G    -1    1   -1
  l3                         C       T     1    0    1
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_ldat(mergefunc)

  with FUDWriter(filename, format, genos.samples, genos.genome, genos.phenome, extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
