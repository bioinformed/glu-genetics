# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'WTCCC genotype format output object'
__copyright__ = 'Copyright (c) 2007-2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['WTCCCWriter', 'save_wtccc']

__genoformats__ = [
  #LOADER        SAVER          WRITER       PFORMAT  ALIAS    EXTS
  ('load_wtccc', 'save_wtccc', 'WTCCCWriter', 'ldat', 'wtccc', None) ]


from   itertools                 import izip, imap

from   glu.lib.fileutils         import autofile,parse_augmented_filename,get_arg,guess_related_file,related_file,trybool,namefile

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.phenos    import Phenome
from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.genoarray import GenotypeArray,build_model,build_descr


GENOS = [ ['0','0','0'],
          ['1','0','0'],
          ['0','1','0'],
          ['0','0','1'] ]


def _encode_wtccc(model):
  assert model is not None

  if len(model.alleles) > 3:
    raise ValueError('WTCCC files support only biallelic models')

  genovalues = [ GENOS[0] ]*4

  allele1,allele2 = (model.alleles[1:]+[None,None])[:2]

  if allele1:
    genovalues[ model[allele1,allele1].index ] = GENOS[1]
  if allele2:
    genovalues[ model[allele2,allele2].index ] = GENOS[3]
  if allele1 and allele2:
    genovalues[ model[allele1,allele2].index ] = GENOS[2]

  return allele1,allele2,genovalues


def _decode_wtccc(model,a,b,threshold,fields):
  if threshold<0.50:
    raise ValueError
  aa = model[a,a]
  ab = model[a,b]
  bb = model[b,b]
  nn = model[None,None]

  for p_aa,p_ab,p_bb in izip(*[imap(float,fields)]*3):
    if p_aa>threshold:
      geno = aa
    elif p_ab>threshold:
      geno = ab
    elif p_bb>threshold:
      geno = bb
    else:
      geno = nn

    yield geno


def load_wtccc_samples(filename):
  sfile = autofile(filename)
  header = sfile.next().split()

  if header[:3] != ['ID_1','ID_2','missing']:
    raise ValueError('Invalid header found for WTCCC sample file %s' % namefile(filename))

  header2 = sfile.next().split()

  if header2[:3] != ['0','0','0']:
    raise ValueError('Invalid header2 found for WTCCC sample file %s' % namefile(filename))

  for line_num,line in enumerate(sfile):
    fields = line.split()

    if fields[0:1] and fields[1:2]:
      sample = '%s:%s' % (fields[0],fields[1])
    else:
      sample = fields[0] or (len(fields)>1 and fields[1])

    if not sample:
      raise ValueError('Invalid WTCCC sample data on line %d of %s' % (line_num+1,namefile(filename)))

    yield sample


def load_wtccc(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html

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
  ... '0 l1 0 A G 1 0 0 0 1 0 0 0 1\\n'
  ... '0 l2 0 C G 0 0 0 0 1 0 0 0 0\\n'
  ... '0 l3 0 C T 0 1 0 1 0 0 0 1 0\\n')
  >>> samples = StringIO('ID_1 ID_2 missing\\n0 0 0 0 0 0\\nf1 s1 0\\nf1 s2 0\\nf1 s3 0\\n')
  >>> genos = load_wtccc(data,'wtccc',samples=samples)
  >>> genos.format
  'ldat'
  >>> genos.columns
  ('f1:s1', 'f1:s2', 'f1:s3')
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
  samples   = get_arg(args, ['samples','s']) or guess_related_file(filename,['lst'])
  threshold = float(get_arg(args, ['threshold','t'], 0.95))

  if samples is None:
    raise ValueError('Sample file must be specified when loading WTCCC files')

  if threshold<0.50:
    raise ValueError('WTCCC genotype confidence threshold must be greater than 0.5')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  samples = list(load_wtccc_samples(samples))
  gfile   = autofile(filename)

  if phenome is None:
    phenome = Phenome()

  loci   = []
  models = []
  file_genome = Genome()

  def _load_wtccc():
    m           = len(samples)
    n           = 5+m*3
    modelcache  = {}

    for line_num,row in enumerate(gfile):
      fields = row.split()

      if len(fields) != n:
        raise ValueError('Invalid WTCCC row on line %d of %s' % (line_num+1,namefile(filename)))

      chrom   = fields[0]
      lname   = fields[1]
      loc     = fields[2]
      a,b     = fields[3:5]
      alleles = a,b

      if chrom.startswith('chr'):
        chrom = chrom[3:]

      if chrom.upper()=='MT':
        chrom = 'M'
      elif chrom not in ('X','Y','XY','M'):
        try:
          chrom = int(chrom)
          if chrom < 1 or chrom > 22:
            chrom = None
        except ValueError:
          chrom = None

      try:
        loc = int(loc) or None
      except ValueError:
        loc = None

      model = modelcache.get(alleles)
      if model is None:
        model = modelcache[alleles] = build_model(alleles=alleles,max_alleles=2)

      descr = build_descr(model,m)

      loci.append(lname)
      models.append(model)
      file_genome.set_locus(lname,model,chr,loc)

      genos = GenotypeArray(descr, _decode_wtccc(model,a,b,threshold,fields[5:]))

      yield lname,genos

  genos = GenomatrixStream(_load_wtccc(),'ldat',samples=samples,loci=loci,models=models,
                                                genome=file_genome,phenome=phenome,
                                                unique=unique,packed=True)

  if unique:
    genos = genos.unique_checked()

  if genome:
    genos = genos.transformed(recode_models=genome)

  return genos


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
  >>> with WTCCCWriter(o,'wtccc',genos.samples,genos.genome,genos.phenome) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 A G 1 0 0 0 1 0 0 0 1
  0 l2 0 C G 0 0 0 0 1 0 0 0 0
  0 l3 0 C T 0 1 0 1 0 0 0 1 0
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

    samplefile  = get_arg(args, ['samplefile','samples','sample'])

    # Careful: samplefile=<blank> is intended to suppress output
    if samplefile is None:
      samplefile = related_file(filename,'sample')

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    self.out        = autofile(filename,'wb')
    self.samples    = samples
    self.genome     = genome
    self.samplefile = samplefile

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

    if self.samplefile:
      out = autofile(self.samplefile,'wb')
      out.write('ID_1 ID_2 missing\r\n')
      out.write('0 0 0\r\n')
      for sample in self.samples:
        out.write('%s %s 0\r\n' % (sample,sample))

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


def save_wtccc(filename,genos,format,extra_args=None,**kwargs):
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
  >>> save_wtccc(o,genos,'wtccc')
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

  with WTCCCWriter(filename, format, genos.samples, genos.genome, genos.phenome, extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
