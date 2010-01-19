# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'PLINK genotype format input/output objects'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['PlinkPedWriter',  'save_plink_ped',  'load_plink_ped',
                 'PlinkTPedWriter', 'save_plink_tped', 'load_plink_tped',
                 'PlinkBedWriter',  'save_plink_bed',  'load_plink_bed']

__genoformats__ = [
  #      LOADER           SAVER             WRITER        PFORMAT   ALIAS                      EXTS
  ('load_plink_ped', 'save_plink_ped', 'PlinkPedWriter',  'sdat',  'plink_ped',               'ped'),
  ('load_plink_tped','save_plink_tped','PlinkTPedWriter', 'ldat',  'plink_tped',             'tped'),
  ('load_plink_bed', 'save_plink_bed', 'PlinkBedWriter',  'ldat', ['plink_bed','lbed',
                                                                   'plink_lbed'],             'bed'),
  ('load_plink_bed', 'save_plink_bed', 'PlinkBedWriter',  'sdat', ['plink_sbed','sbed'],      None ) ]


import string

from   operator                  import getitem
from   itertools                 import islice,izip

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,               \
                                        guess_related_file,related_file, \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import build_model, build_descr, GenotypeArray, GenotypeArrayDescriptor
from   glu.lib.genolib.locus     import Genome,load_locus_records,populate_genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN, \
                                        PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


ALLELE_MAP  = {'0':None,'1':'1','2':'2','A':'A','C':'C','G':'G','T':'T','B':'B','+':'+','-':'-'}
ALLELE_RMAP = {None:'0'}

SEX_MAP    = {'1':SEX_MALE,'2':SEX_FEMALE}
SEX_RMAP   = {SEX_UNKNOWN:'0', SEX_MALE:'1', SEX_FEMALE:'2'}

PHENO_MAP  = {'1':PHENO_UNAFFECTED, '2':PHENO_AFFECTED}
PHENO_RMAP = {PHENO_UNKNOWN:'0',PHENO_UNAFFECTED:'1',PHENO_AFFECTED:'2'}

PARENT_MAP = {'0':None}

CHR_MAP    = {'0':None,'23':'X','24':'Y','25':'XY','26':'M'}
CHR_RMAP   = {None:'0','':'0','X':'23','Y':'24','XY':'25','M':'26','MT':'26'}


# Only ASCII is currently supported
_ident_ws_trans = string.maketrans('\t\n\x0b\x0c\r ','______')


def escape_ident_ws(ident):
  '''
  Strip leading and trailing whitespace and replace all remaining whitespace
  with _ in identifiers.  Only ASCII/Latin1 is currently supported.  This
  will not work if Unicode identifiers are ever allowed.

  >>> escape_ident_ws(' BW 012 213 ')
  'BW_012_213'
  >>> escape_ident_ws('foo\\tbar\\r')
  'foo_bar'
  '''
  return string.translate(ident.strip(),_ident_ws_trans)


def load_plink_map(filename,genome):
  mfile = autofile(filename)

  for i,line in enumerate(mfile):
    line = line.rstrip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    if len(fields) != 4:
      raise ValueError('Invalid PLINK MAP record %d' % (i+1))

    chr   = fields[0]
    lname = fields[1]
    gdist = int(fields[2])      if fields[2] else None
    pdist = abs(int(fields[3])) if fields[3] else None

    if not lname:
      raise ValueError('Invalid PLINK MAP record %d' % (i+1))

    if chr == '0':
      chr = None

    genome.merge_locus(lname, chromosome=chr, location=pdist)

    yield lname


def load_plink_ped(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a PLINK PED format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param      phenome: phenome descriptor
  @type       phenome: Phenome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
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

  unique = get_arg(args, ['unique'], True)
  loci   = get_arg(args, ['loci'])
  lmap   = get_arg(args, ['map']) or guess_related_file(filename,['map'])

  if loci is None and lmap is None:
    raise ValueError('Map file or locus list must be specified when loading PLINK PED files')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if genome is None:
    genome = Genome()

  if phenome is None:
    phenome = Phenome()

  if loci and isinstance(loci,basestring):
    loci = list(load_locus_records(loci)[2])
    # Merge map data into genome
    populate_genome(genome,loci)
    loci = [ intern(l[0]) for l in loci ]

  if lmap:
    map_loci = list(load_plink_map(lmap,genome))
    if not loci:
      loci = map_loci
    elif loci != map_loci:
      raise ValueError('Locus list and PLINK MAP file are not identical')

  loci = loci or []

  gfile = autofile(filename)
  n     = 6 + 2*len(loci)

  def _load_plink():
    aget = ALLELE_MAP.get

    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled():
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        family,name,father,mother,sex,pheno = fields[:6]

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

        a1 = [ aget(a,a) for a in islice(fields,6,None,2) ]
        a2 = [ aget(a,a) for a in islice(fields,7,None,2) ]
        genos = zip(a1,a2)

      yield ename,genos

  genos = GenomatrixStream.from_tuples(_load_plink(),'sdat',loci=loci,genome=genome,phenome=phenome,unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos


class PlinkPedWriter(object):
  '''
  Object to write a matrix data to a PLINK PED file

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> m = StringIO()
  >>> with PlinkPedWriter(o,'ped',genos.loci,genos.genome,genos.phenome,mapfile=m) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 0 0 A A 0 0 C T
  s2 s2 0 0 0 0 A G C G C C
  s3 s3 0 0 0 0 G G 0 0 C T
  >>> print m.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0
  0 l2 0 0
  0 l3 0 0
  '''
  def __init__(self,filename,format,loci,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param         loci: locus names
    @type          loci: list of str
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    mapfile  = get_arg(args, ['mapfile','map'])

    # Careful: mapfile=<blank> is intended to suppress output
    if mapfile is None:
      mapfile = related_file(filename,'map')

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if phenome is None:
      phenome = Phenome()

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

    family     = escape_ident_ws(family or individual)
    individual = escape_ident_ws(individual)
    parent1    = escape_ident_ws(parent1 or '0')
    parent2    = escape_ident_ws(parent2 or '0')

    sex   = SEX_RMAP[phenos.sex]
    pheno = PHENO_RMAP[phenos.phenoclass]

    row = [family,individual,parent1,parent2,sex,pheno]
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

      family     = escape_ident_ws(family or individual)
      individual = escape_ident_ws(individual)
      parent1    = escape_ident_ws(parent1 or '0')
      parent2    = escape_ident_ws(parent2 or '0')

      sex   = SEX_RMAP[phenos.sex]
      pheno = PHENO_RMAP[phenos.phenoclass]

      row = [family,individual,parent1,parent2,sex,pheno]

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

    # FIXME: map file writer should be refactored
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


def save_plink_ped(filename,genos,format,extra_args=None,**kwargs):
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
  >>> save_plink_ped(o,genos,'ped')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 0 0 A A 0 0 C T
  s2 s2 0 0 0 0 A G C G C C
  s3 s3 0 0 0 0 G G 0 0 C T
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_sdat(mergefunc)

  with PlinkPedWriter(filename, format, genos.loci, genos.genome, genos.phenome,
                                extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

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
      raise ValueError('Invalid record on line %d of %s' % (i+1,namefile(filename)))

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


def load_plink_tped(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a PLINK TPED format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param      phenome: phenome descriptor
  @type       phenome: Phenome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
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

  unique = get_arg(args, ['unique'], True)
  loci   = get_arg(args, ['loci'])
  tfam   = get_arg(args, ['tfam','fam']) or guess_related_file(filename,['tfam','fam'])

  if tfam is None:
    raise ValueError('A TFAM file must be specified when loading PLINK TPED data')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if genome is None:
    genome = Genome()

  if loci and isinstance(loci,basestring):
    loci = list(load_locus_records(loci)[2])
    # Merge map data into genome
    populate_genome(genome,loci)
    loci = [ intern(l[0]) for l in loci ]

  loci = loci or []

  samples = list(load_plink_tfam(tfam,phenome))

  gfile = autofile(filename)
  n     = 4 + 2*len(samples)

  def _load_plink():
    aget = ALLELE_MAP.get

    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled():
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        chr   = fields[0] or None
        lname = fields[1]
        gdist = int(fields[2])      if fields[2] else None
        pdist = abs(int(fields[3])) if fields[3] else None

        if not lname:
          raise ValueError('Invalid PLINK TPED record %d' % (line_num+1))

        if chr == '0':
          chr = None

        genome.merge_locus(lname, chromosome=chr, location=pdist)

        a1 = [ aget(a,a) for a in islice(fields,4,None,2) ]
        a2 = [ aget(a,a) for a in islice(fields,5,None,2) ]
        genos = zip(a1,a2)

      yield lname,genos

  genos = GenomatrixStream.from_tuples(_load_plink(),'ldat',samples=samples,genome=genome,phenome=phenome,unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos


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
  >>> with PlinkTPedWriter(o,'tped',genos.samples,genos.genome,genos.phenome,tfamfile=m) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0 A A A G G G
  0 l2 0 0 0 0 C G 0 0
  0 l3 0 0 C T C C C T
  >>> print m.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 0 0
  s2 s2 0 0 0 0
  s3 s3 0 0 0 0
  '''
  def __init__(self,filename,format,samples,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param      samples: sample names
    @type       samples: list of str
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    tfamfile = get_arg(args, ['tfamfile','tfam'])

    # Careful: mapfile=<blank> is intended to suppress output
    if tfamfile is None:
      tfamfile = related_file(filename,'tfam')

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if phenome is None:
      phenome = Phenome()

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

    # FIXME: tfam file writer should be refactored
    if self.tfamfile:
      phenome = self.phenome
      out = autofile(self.tfamfile,'wb')
      for sample in self.samples:
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


def save_plink_tped(filename,genos,format,extra_args=None,**kwargs):
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
  >>> save_plink_tped(o,genos,'tped')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  0 l1 0 0 A A A G G G
  0 l2 0 0 0 0 C G 0 0
  0 l3 0 0 C T C C C T
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_ldat(mergefunc)

  with PlinkTPedWriter(filename, format, genos.samples, genos.genome, genos.phenome,
                                 extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


##############################################################################################


def load_plink_bim(filename,genome):
  mfile = autofile(filename)
  modelcache = {}

  for i,line in enumerate(mfile):
    line = line.rstrip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    if len(fields) != 6:
      raise ValueError('Invalid PLINK BIM record %d' % (i+1))

    chr     = CHR_MAP.get(fields[0].upper(),fields[0])
    locus   = fields[1]
    gdist   = int(fields[2])      if fields[2] else None
    pdist   = abs(int(fields[3])) if fields[3] else None
    allele1 = fields[4]
    allele2 = fields[5]

    if not locus:
      raise ValueError('Invalid PLINK BIM record %d' % (i+1))

    if chr == '0':
      chr = None

    if allele1 == '0':
      allele1 = None

    if allele2 == '0':
      allele2 = None

    key = allele1,allele2

    model = modelcache.get(key)
    if model is None:
      a,b = allele1,allele2
      if a and b:
        genos = [(a,a),(b,b),(a,b)]
      elif a:
        genos = [(a,a)]
      elif b:
        raise RuntimeError('Invalid BIM locus model: (%s/%s)' % (a,b))
      else:
        genos = []

      model = modelcache[key] = build_model(genotypes=genos,max_alleles=2)

    genome.merge_locus(locus, model, chr, pdist)

    yield locus,model


def _plink_encode(model, shift=0):
  genos = [1<<shift]*4

  _plink_update_encoding(model, genos, shift)

  return genos


def _plink_update_encoding(model, genos, shift):
  n  = len(model.alleles)

  if n > 3:
    raise ValueError('PLINK BIM files support only allelic models')

  a1 = model.alleles[1] if n>1 else None
  a2 = model.alleles[2] if n>2 else None

  if a1:
    genos[ model[a1,a1].index ] = 0
  if a2:
    genos[ model[a2,a2].index ] = 3<<shift
  if a1 and a2:
    genos[ model[a1,a2].index ] = 2<<shift


def load_plink_bed(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a PLINK BED format genotype data file.

  See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

  Use undocumented "--make-bed --ind-major" PLINK options to get sdat flavor.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param      phenome: phenome descriptor
  @type       phenome: Phenome instance
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream

  PLINK defines its binary genotype encoding based on an in-memory bit
  pattern, which is written in reverse to disk.  Due to GLU's requirement
  that the missing genotype be encoded as 0, we simply need to xor the
  on-disk representation by 0b01010101 (0x55) to obtain a valid encoding
  with homozygote 1 encoded as 1, homozygote 2 encoded as 2, and
  heterozygotes as 3.  The resulting homogeneous encoding can then be
  remapped to the correct descriptor.

                   geno-  in-    on-   01
                   type  memory disk  xor'd
  ---------------- ----  ------ ----  -----
  Missing genotype  __    10     01     00
  Homozygote 1      AA    00     00     01
  Homozygote 2      BB    11     11     10
  Heterozygote      AB    01     10     11
  '''
  import numpy  as np

  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  gfile = autofile(filename)

  magic = map(ord,gfile.read(2))
  mode  = ord(gfile.read(1))

  if magic != [0x6c,0x1b]:
    raise ValueError('Invalid PLINK BED file magic number')

  if mode not in [0,1]:
    raise ValueError('Invalid PLINK BED file mode')

  loc = get_arg(args, ['loci'])
  bim = get_arg(args, ['map','bim' ]) or guess_related_file(filename,['bim'])
  fam = get_arg(args, ['fam','tfam']) or guess_related_file(filename,['fam','tfam'])

  if bim is None:
    raise ValueError('BIM file must be specified when loading PLINK BED data')

  if fam is None:
    raise ValueError('A FAM file must be specified when loading PLINK BED data')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if phenome is None:
    phenome = Phenome()

  file_genome = Genome()

  loci,models = zip(*load_plink_bim(bim,file_genome))
  models      = list(models)
  samples     = list(load_plink_tfam(fam,phenome))

  if loc and isinstance(loc,basestring):
    loc = list(load_locus_records(loc)[2])
    # Merge map data into genome
    populate_genome(file_genome,loc)

  unique = len(set(samples))==len(samples) and len(set(loci))==len(loci)

  def _decode(a):
    '''
    Reverse pairs of bits and xor result by 0b0101010101
    '''
    return ((a&0x03)<<6 | (a&0x0C)<<2 | (a&0x30)>>2 | (a&0xC0)>>6)^0x55

  decode = dict( (chr(i),_decode(i)) for i in range(256) )

  if mode == 0:
    format = 'sdat'

    def _load_plink():
      descr = GenotypeArrayDescriptor(models)
      rowbytes = (len(loci)*2+7)//8
      decoder = [decode]*rowbytes

      for sample in samples:
        genos = GenotypeArray(descr)
        genos.data = np.array(map(getitem, decoder, gfile.read(rowbytes)),dtype=np.uint8)
        yield sample,genos

  elif mode == 1:
    format = 'ldat'

    def _load_plink():
      n = len(samples)
      rowbytes = (n*2+7)//8
      decoder = [decode]*rowbytes

      for lname,model in izip(loci,models):
        descr = build_descr(model,n)
        genos = GenotypeArray(descr)
        genos.data = np.array(map(getitem, decoder, gfile.read(rowbytes)),dtype=np.uint8)
        yield lname,genos

  genos = GenomatrixStream(_load_plink(),format,loci=loci,samples=samples,models=models,
                                         genome=file_genome,phenome=phenome,unique=unique,
                                         packed=True)

  if genome:
    genos = genos.transformed(recode_models=genome)

  return genos


class PlinkBedWriter(object):
  '''
  Object to write a PLINK BED file

  See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

  Example of writing an sdat file:

  >>> loci =         (    'l1',       'l2',        'l3'  )
  >>> rows = [('s1', ( ('A', 'A'), (None,None),  ('T','T'))),
  ...         ('s2', ((None,None),  ('C','T'),   ('G','T'))),
  ...         ('s3', ( ('A', 'T'),  ('T','C'),   ('G','G')))]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> import tempfile
  >>> bed = tempfile.NamedTemporaryFile()
  >>> bim = tempfile.NamedTemporaryFile()
  >>> fam = tempfile.NamedTemporaryFile()
  >>> with PlinkBedWriter(bed.name,'sbed',genos.columns,genos.genome,genos.phenome,
  ...                     bim=bim.name,fam=fam.name) as writer:
  ...   genos=iter(genos)
  ...   writer.writerow(*genos.next())
  ...   writer.writerow(*genos.next())
  ...   writer.writerows(genos)
  >>> genos = load_plink_bed(bed.name,'bed',bim=bim.name,fam=fam.name)
  >>> genos.format
  'sdat'
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1:s1', [('A', 'A'), (None, None), ('T', 'T')])
  ('s2:s2', [(None, None), ('C', 'T'), ('G', 'T')])
  ('s3:s3', [('A', 'T'), ('C', 'T'), ('G', 'G')])

  Example of writing an ldat file:

  >>> samples =         (    's1',       's2',       's3'   )
  >>> rows    = [('l1', ( ('A', 'A'), (None,None),  ('T','T'))),
  ...            ('l2', ((None,None),  ('T','T'),   ('G','T'))),
  ...            ('l3', ( ('A', 'T'),  ('T','A'),   ('T','T')))]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> with PlinkBedWriter(bed.name,'lbed',genos.columns,genos.genome,genos.phenome,
  ...                     bim=bim.name,fam=fam.name) as writer:
  ...   genos=iter(genos)
  ...   writer.writerow(*genos.next())
  ...   writer.writerow(*genos.next())
  ...   writer.writerows(genos)
  >>> genos = load_plink_bed(bed.name,'bed',bim=bim.name,fam=fam.name)
  >>> genos.format
  'ldat'
  >>> genos.samples
  ('s1:s1', 's2:s2', 's3:s3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('T', 'T')])
  ('l2', [(None, None), ('T', 'T'), ('G', 'T')])
  ('l3', [('A', 'T'), ('A', 'T'), ('T', 'T')])
  '''
  def __init__(self,filename,format,columns,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param      samples: sample names
    @type       samples: list of str
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    '''
    if format not in ('bed','sbed','lbed'):
      raise ValueError('format must be either bed, sbed, or lbed, found %s' % format)
      raise IOError('format must be either ldat or sdat')

    if format=='bed':
      format=='lbed'

    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    bimfile = get_arg(args, ['bim'])
    famfile = get_arg(args, ['fam'])

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    # Careful: file=<blank> is intended to suppress output
    if famfile is None:
      famfile = related_file(filename,'fam')
    if bimfile is None:
      bimfile = related_file(filename,'bim')

    if phenome is None:
      phenome = Phenome()

    self.out     = autofile(filename,'wb')
    self.format  = format
    self.columns = columns
    self.rowkeys = []
    self.genome  = genome
    self.phenome = phenome
    self.famfile = famfile
    self.bimfile = bimfile

    # Write magic number
    self.out.write( ''.join( map(chr,[0x6c,0x1b]) ) )

    # Write magic number and mode=1
    if format=='lbed':
      # Write BED normal mode (ldat)
      self.out.write( chr(1) )
    else:
      # Write BED transposed-mode (sdat)
      self.out.write( chr(0) )

  def writerow(self, rowkey, genos):
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

    n = len(self.columns)
    if len(genos) != n:
      raise ValueError('[ERROR] Internal error: Genotypes do not match columns')

    rowbytes = (n*2+7)//8

    if self.format == 'lbed':
      # Build encodings for each of the 4 shift values for each model
      model  = self.genome.get_model(rowkey)
      values = [ _plink_encode(model,shift) for shift in [0,2,4,6] ]

      row = [0]*rowbytes
      for i,g in enumerate(genos):
        row[i//4] |= values[i%4][g.index]
    else:
      # Build encodings for each locus and store the set of cached encodings
      # so they can be updated as new alleles are detected.
      #
      # FIXME: This approach is awful and should be replaced by something
      #        much more intelligent
      valuecache = {}
      genovalues = []

      for i,locus in enumerate(self.columns):
        model  = self.genome.get_model(locus)
        shift  = 2*(i%4)
        key    = model,shift
        values = valuecache.get(key)
        if values is None:
          values = valuecache[key] = _plink_encode(model,shift)
        genovalues.append(values)

      row = [0]*rowbytes
      for i,g in enumerate(genos):
        row[i//4] |= genovalues[i][g.index]

    # Write row
    out.write( ''.join(map(chr,row)) )
    self.rowkeys.append(rowkey)

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

    n = len(self.columns)
    rowbytes = (n*2+7)//8

    if self.format=='lbed':
      for lname,genos in rows:
        if len(genos) != n:
          raise ValueError('[ERROR] Internal error: Genotypes do not match columns')

        # Build encodings for each of the 4 shift values for each model
        model  = self.genome.get_model(lname)
        values = [ _plink_encode(model,shift) for shift in [0,2,4,6] ]

        row = [0]*rowbytes
        for i,g in enumerate(genos):
          row[i//4] |= values[i%4][g.index]

        # Write row
        out.write( ''.join(map(chr,row)) )
        self.rowkeys.append(lname)

    else:
      valuecache = {}

      for sample,genos in rows:
        if len(genos) != n:
          raise ValueError('[ERROR] Internal error: Genotypes do not match columns')

        # Build encodings for each locus and store the set of cached encodings
        # so they can be updated as new alleles are detected.
        #
        # FIXME: This approach is awful and should be replaced by something
        #        much more intelligent
        genovalues = []
        for i,locus in enumerate(self.columns):
          model  = self.genome.get_model(locus)
          shift  = 2*(i%4)
          key    = model,shift
          values = valuecache.get(key)
          if values is None:
            values = valuecache[key] = _plink_encode(model,shift)
          genovalues.append(values)

        row = [0]*rowbytes
        for i,g in enumerate(genos):
          row[i//4] |= genovalues[i][g.index]

        # Write row
        out.write( ''.join(map(chr,row)) )
        self.rowkeys.append(sample)

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')

    self.out.close()
    self.out = None

    if self.format=='lbed':
      loci,samples = self.rowkeys,self.columns
    else:
      loci,samples = self.columns,self.rowkeys

    # FIXME: fam file writer should be refactored
    if self.famfile:
      phenome = self.phenome
      out = autofile(self.famfile,'wb')
      for sample in samples:
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

        family     = escape_ident_ws(family or individual)
        individual = escape_ident_ws(individual)
        parent1    = escape_ident_ws(parent1 or '0')
        parent2    = escape_ident_ws(parent2 or '0')

        sex   = SEX_RMAP[phenos.sex]
        pheno = PHENO_RMAP[phenos.phenoclass]

        row = [family,individual,parent1,parent2,sex,pheno]
        out.write( ' '.join(row))
        out.write('\r\n')

    # FIXME: bim file writer should be refactored
    if self.bimfile:
      out = autofile(self.bimfile,'wb')
      for locus in loci:
        loc   = self.genome.get_locus(locus)
        chrom = loc.chromosome or ''
        chrom = CHR_RMAP.get(chrom.upper(),chrom)
        a1,a2 = (loc.model.alleles[1:]+[0,0])[:2]

        out.write( ' '.join( map(str,[chrom,locus,0,loc.location or 0, a1, a2] )) )
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


def save_plink_bed(filename,genos,format,extra_args=None,**kwargs):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence

  Example of writing an sdat file:

  >>> loci =         (    'l1',       'l2',        'l3'  )
  >>> rows = [('s1', ( ('A', 'A'), (None,None),  ('T','T'))),
  ...         ('s2', ((None,None),  ('C','T'),   ('G','T'))),
  ...         ('s3', ( ('A', 'T'),  ('T','C'),   ('G','G')))]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> import tempfile
  >>> bed = tempfile.NamedTemporaryFile()
  >>> bim = tempfile.NamedTemporaryFile()
  >>> fam = tempfile.NamedTemporaryFile()
  >>> save_plink_bed(bed.name,genos,'bed',bim=bim.name,fam=fam.name)
  >>> genos = load_plink_bed(bed.name,'bed',bim=bim.name,fam=fam.name)
  >>> genos.format
  'sdat'
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1:s1', [('A', 'A'), (None, None), ('T', 'T')])
  ('s2:s2', [(None, None), ('C', 'T'), ('G', 'T')])
  ('s3:s3', [('A', 'T'), ('C', 'T'), ('G', 'G')])

  Example of writing an ldat file:

  >>> samples =         (    's1',       's2',       's3'   )
  >>> rows    = [('l1', ( ('A', 'A'), (None,None),  ('T','T'))),
  ...            ('l2', ((None,None),  ('T','T'),   ('G','T'))),
  ...            ('l3', ( ('A', 'T'),  ('T','A'),   ('T','T')))]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> save_plink_bed(bed.name,genos,'bed',bim=bim.name,fam=fam.name)
  >>> genos = load_plink_bed(bed.name,'bed',bim=bim.name,fam=fam.name)
  >>> genos.format
  'ldat'
  >>> genos.samples
  ('s1:s1', 's2:s2', 's3:s3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('T', 'T')])
  ('l2', [(None, None), ('T', 'T'), ('G', 'T')])
  ('l3', [('A', 'T'), ('A', 'T'), ('T', 'T')])
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  if format=='bed':
    if genos.format=='genotriple':
      genos = genos.as_ldat(mergefunc)
      format = 'lbed'
    elif genos.format=='sdat':
      format = 'sbed'
    elif genos.format=='ldat':
      format = 'lbed'
    else:
      raise ValueError('Invalid genotype format')

    if mergefunc is not None:
      genos = genos.merged(mergefunc)

  elif format=='lbed':
    genos = genos.as_ldat(mergefunc)
  elif format=='sbed':
    genos = genos.as_sdat(mergefunc)
  else:
    raise ValueError('Invalid genotype format')

  with PlinkBedWriter(filename, format, genos.columns, genos.genome, genos.phenome,
                                extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


###############################################################################################


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
