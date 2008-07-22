# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'PLINK genotype format input/output objects'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import string

from   itertools                 import islice,izip

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,               \
                                        guess_related_file,related_file, \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import model_from_alleles
from   glu.lib.genolib.locus     import Genome,load_locus_records,populate_genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN, \
                                        PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


__all__ = ['PlinkPedWriter',  'save_plink_ped',  'load_plink_ped',
           'PlinkTPedWriter', 'save_plink_tped', 'load_plink_tped',
           'PlinkBedWriter',  'save_plink_bed',  'load_plink_bed']


ALLELE_MAP  = {'0':None,'1':'1','2':'2','A':'A','C':'C','G':'G','T':'T','B':'B','+':'+','-':'-'}
ALLELE_RMAP = {None:'0'}

SEX_MAP    = {'1':SEX_MALE,'2':SEX_FEMALE}
SEX_RMAP   = {SEX_UNKNOWN:'0', SEX_MALE:'1', SEX_FEMALE:'2'}

PHENO_MAP  = {'1':PHENO_UNAFFECTED, '2':PHENO_AFFECTED}
PHENO_RMAP = {PHENO_UNKNOWN:'0',PHENO_UNAFFECTED:'1',PHENO_AFFECTED:'2'}

PARENT_MAP = {'0':None}

CHR_MAP    = {'0':None,'23':'X','24':'Y','25':'XY','26':'M'}
CHR_RMAP   = {None:'0','':'0','X':'23','Y':'24','XY':'25','M':'26','MT':26}


# Only ASCII is currently supported
_ident_ws_trans = string.maketrans('\t\n\x0b\x0c\r ','______')


# FIXME: Merlin also needs similar escaping (gack, I hate whitespace delimited formats)
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


def load_plink_ped(filename,genome=None,phenome=None,extra_args=None,**kwargs):
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
  >>> with PlinkPedWriter(o,genos.loci,genos.genome,genos.phenome,mapfile=m) as w:
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
  def __init__(self,filename,loci,genome,phenome,extra_args=None,**kwargs):
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

    # Careful: mapfile=<blank> is intended to supress output
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


def save_plink_ped(filename,genos,extra_args=None,**kwargs):
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

  with PlinkPedWriter(filename, genos.loci, genos.genome, genos.phenome,
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


def load_plink_tped(filename,genome=None,phenome=None,extra_args=None,**kwargs):
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
  >>> with PlinkTPedWriter(o,genos.samples,genos.genome,genos.phenome,tfamfile=m) as w:
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
  def __init__(self,filename,samples,genome,phenome,extra_args=None,**kwargs):
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

    # Careful: mapfile=<blank> is intended to supress output
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


def save_plink_tped(filename,genos,extra_args=None,**kwargs):
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
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_ldat(mergefunc)

  with PlinkTPedWriter(filename, genos.samples, genos.genome, genos.phenome,
                                 extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


##############################################################################################


def load_plink_bim(filename,genome):
  mfile = autofile(filename)

  # FIXME: Add support for hemizygote models
  m = genome.max_alleles+1
  modelcache = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                           if locus.model is not None and len(locus.model.alleles)==m )

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

    alleles = []

    if allele1 == '0':
      allele1 = None
    else:
      alleles.append(allele1)

    if allele2 == '0':
      allele2 = None
    else:
      alleles.append(allele2)

    alleles = tuple(sorted(alleles))

    model = genome.get_locus(locus).model

    if not model:
      model = modelcache.get(alleles)

    if not model:
      model = model_from_alleles(alleles,max_alleles=2)

      if genome.default_model is None and len(alleles) == genome.max_alleles:
        modelcache[alleles] = model

    genome.merge_locus(locus, model, True, chr, pdist)

    yield locus,allele1,allele2


def _plink_decode(model,allele1,allele2):
  genos = [model.genotypes[0]]*4

  if allele1:
    genos[0] = model.add_genotype( (allele1,allele1) )
  if allele2:
    genos[3] = model.add_genotype( (allele2,allele2) )
  if allele1 and allele2:
    genos[2] = model.add_genotype( (allele1,allele2) )

  return genos


def _plink_encode(model, shift=0):
  genos = [1<<shift]*4

  _plink_update_encoding(model, genos, shift)

  return genos


def _plink_update_encoding(model, genos, shift):
  n  = len(model.alleles)

  if n > 3:
    raise ValueError('PLINK BIM files support only biallic models')

  a1 = model.alleles[1] if n>1 else None
  a2 = model.alleles[2] if n>2 else None

  if a1:
    genos[ model.add_genotype( (a1,a1) ).index ] = 0
  if a2:
    genos[ model.add_genotype( (a2,a2) ).index ] = 3<<shift
  if a1 and a2:
    genos[ model.add_genotype( (a1,a2) ).index ] = 2<<shift


def load_plink_bed(filename,genome=None,phenome=None,extra_args=None,**kwargs):
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
  '''
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

  if genome is None:
    genome = Genome()

  if phenome is None:
    phenome = Phenome()

  bim_loci = list(load_plink_bim(bim,genome))
  loci     = [ l[0] for l in bim_loci ]
  samples  = list(load_plink_tfam(fam,phenome))
  models   = [ genome.get_model(locus) for locus in loci ]

  if loc and isinstance(loc,basestring):
    loc = list(load_locus_records(loc)[2])
    # Merge map data into genome
    populate_genome(genome,loc)

  unique = len(set(samples))==len(samples) and len(set(loci))==len(loci)

  if mode == 0:
    format = 'sdat'

    def _load_plink():
      valuecache = {}
      genovalues = []

      for i,(locus,allele1,allele2) in enumerate(bim_loci):
        # Cache must key off model and allele order
        model  = models[i]
        key    = model,allele1,allele2

        values = valuecache.get(key)

        if values is None:
          values = valuecache[key] = _plink_decode(*key)

        byte  = i//4
        shift = 2*(i%4)

        genovalues.append( (byte,shift,values) )

      rowbytes = (len(loci)*2+7)//8

      for sample in samples:
        data  = map(ord,gfile.read(rowbytes))
        genos = [ values[(data[byte]>>shift)&3] for byte,shift,values in genovalues ]
        yield sample,genos

  elif mode == 1:
    format = 'ldat'

    def _load_plink():
      genovalues = []

      for i,sample in enumerate(samples):
        byte  = i//4
        shift = 2*(i%4)
        genovalues.append( (byte,shift) )

      valuecache = {}
      rowbytes = (len(samples)*2+7)//8

      for (locus,allele1,allele2),model in izip(bim_loci,models):
        # Cache must key off model and allele order
        key    = model,allele1,allele2
        values = valuecache.get(key)

        if values is None:
          values = valuecache[key] = _plink_decode(*key)

        data = map(ord,gfile.read(rowbytes))
        genos = [ values[(data[byte]>>shift)&3] for byte,shift in genovalues ]

        yield locus,genos

  return GenomatrixStream(_load_plink(),format,loci=loci,samples=samples,models=models,
                                        genome=genome,phenome=phenome,unique=unique)


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
  >>> with PlinkBedWriter(bed.name,genos.format,genos.columns,genos.genome,genos.phenome,
  ...                     bim=bim.name,fam=fam.name) as writer:
  ...   genos=iter(genos)
  ...   writer.writerow(*genos.next())
  ...   writer.writerow(*genos.next())
  ...   writer.writerows(genos)
  >>> genos = load_plink_bed(bed.name,bim=bim.name,fam=fam.name)
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
  >>> with PlinkBedWriter(bed.name,genos.format,genos.columns,genos.genome,genos.phenome,
  ...                     bim=bim.name,fam=fam.name) as writer:
  ...   genos=iter(genos)
  ...   writer.writerow(*genos.next())
  ...   writer.writerow(*genos.next())
  ...   writer.writerows(genos)
  >>> genos = load_plink_bed(bed.name,bim=bim.name,fam=fam.name)
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
    if format not in ('ldat','sdat'):
      raise IOError('format must be either ldat or sdat')

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

    # Careful: file=<blank> is intended to supress output
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
    if format=='ldat':
      # Write BED mode
      self.out.write( chr(1) )
    else:
      # Write BED mode
      self.out.write( chr(0) )
      models = [ self.genome.get_model(locus) for locus in self.columns ]

      # Build encodings for each locus and store the set of cached encodings
      # so they can be updated as new alleles are detected.
      #
      # FIXME: This approach is optimal if all alleles are known, but falls
      #        apart otherwise.  A better approach would be to use an
      #        optimistic dict based approach that resolves new genotypes
      #        only when encoding results in a KeyError.
      #
      # FIXME: Even using this appoach, we should be able to remove fixed
      #        models from the cache to avoid unnecessary updates.
      self.valuecache = valuecache = {}
      self.genovalues = genovalues = []

      for i,model in enumerate(models):
        shift  = 2*(i%4)
        key    = model,shift
        values = valuecache.get(key)
        if values is None:
          values = valuecache[key] = _plink_encode(model,shift)
        genovalues.append(values)

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
    row      = [0]*rowbytes

    if self.format == 'ldat':
      # Build encodings for each of the 4 shift values for each model
      model  = self.genome.get_model(rowkey)
      values = [ _plink_encode(model,shift) for shift in [0,2,4,6] ]

      for i,g in enumerate(genos):
        row[i//4] |= values[i%4][g.index]
    else:
      # FIXME: This update encodings once per model,shift value, but it is
      #        still fairly blecherously slow.
      for (model,shift),values in self.valuecache.iteritems():
        _plink_update_encoding(model, values, shift)

      values = self.genovalues

      for i,g in enumerate(genos):
        row[i//4] |= values[i][g.index]

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

    if self.format=='sdat':
      values = self.genovalues

    for rowkey,genos in rows:
      if len(genos) != n:
        raise ValueError('[ERROR] Internal error: Genotypes do not match columns')

      row = [0]*rowbytes

      if self.format=='ldat':
        # Build encodings for each of the 4 shift values for each model
        model  = self.genome.get_model(rowkey)
        values = [ _plink_encode(model,shift) for shift in [0,2,4,6] ]

        for i,g in enumerate(genos):
          row[i//4] |= values[i%4][g.index]
      else:
        # FIXME: This updates once per model,shift value, but it is still
        #        fairly blecherously slow.
        for (model,shift),values in self.valuecache.iteritems():
          _plink_update_encoding(model, values, shift)

        values = self.genovalues

        for i,g in enumerate(genos):
          row[i//4] |= values[i][g.index]

      # Write row
      out.write( ''.join(map(chr,row)) )
      self.rowkeys.append(rowkey)

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

    if self.format=='ldat':
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


def save_plink_bed(filename,genos,extra_args=None,**kwargs):
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
  >>> save_plink_bed(bed.name,genos,bim=bim.name,fam=fam.name)
  >>> genos = load_plink_bed(bed.name,bim=bim.name,fam=fam.name)
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
  >>> save_plink_bed(bed.name,genos,bim=bim.name,fam=fam.name)
  >>> genos = load_plink_bed(bed.name,bim=bim.name,fam=fam.name)
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

  if genos.format not in ('sdat','ldat'):
    genos = genos.as_ldat(mergefunc)
  elif mergefunc is not None:
    genos = genos.merged(mergefunc)

  with PlinkBedWriter(filename, genos.format, genos.columns, genos.genome, genos.phenome,
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
