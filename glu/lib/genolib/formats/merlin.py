# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'Merlin/MACH genotype format input/output objects'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   itertools                 import islice

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,               \
                                        guess_related_file,related_file, \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.locus     import Genome, load_locus_records, populate_genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN


__all__ = ['MerlinWriter', 'save_merlin', 'load_merlin']


ALLELE_MAP  = {'0':None}
ALLELE_RMAP = {None:'0'}

SEX_MAP     = {'1':SEX_MALE,'2':SEX_FEMALE}
SEX_RMAP    = {SEX_UNKNOWN:'0', SEX_MALE:'M', SEX_FEMALE:'F'}

PARENT_MAP  = {'0':None}


def load_merlin_dat(filename,genome):
  mfile = autofile(filename)

  for i,line in enumerate(mfile):
    line = line.rstrip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    if len(fields) != 2:
      raise ValueError('Invalid Merlin data record %d' % (i+1))

    if fields[0] != 'M':
      raise ValueError('Merlin file reader only supports marker (M) records')

    lname = fields[1]
    genome.merge_locus(lname)

    yield lname


def load_merlin(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a Merlin format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
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
  dat    = get_arg(args, ['dat','data']) or guess_related_file(filename,['dat'])

  if loci is None and dat is None:
    raise ValueError('Dat file or locus file must be specified when loading Merlin files')

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

  if dat:
    map_loci = list(load_merlin_dat(dat,genome))
    if not loci:
      loci = map_loci
    elif loci != map_loci:
      raise ValueError('Locus list and MERLIN DAT file are not identical')

  loci = loci or []

  gfile = autofile(filename)
  n     = 5 + 2*len(loci)

  def _load_merlin():
    aget = ALLELE_MAP.get

    for line_num,line in enumerate(gfile):
      if not line or line.startswith('#'):
        continue

      with gcdisabled():
        fields = line.split()

        if len(fields) != n:
          print len(fields)
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        family,name,father,mother,sex = [ s.strip() for s in fields[:5] ]

        if name == '0':
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        father  = PARENT_MAP.get(father,father)
        mother  = PARENT_MAP.get(mother,mother)

        ename   = '%s:%s' % (family,name)
        efather = '%s:%s' % (family,father) if father else None
        emother = '%s:%s' % (family,mother) if mother else None
        sex     = SEX_MAP.get(sex,SEX_UNKNOWN)

        if father:
          phenome.merge_phenos(efather, family, father, sex=SEX_MALE)
        if mother:
          phenome.merge_phenos(emother, family, mother, sex=SEX_FEMALE)

        phenome.merge_phenos(ename, family, name, efather, emother, sex)

        a1 = [ aget(a,a) for a in islice(fields,5,None,2) ]
        a2 = [ aget(a,a) for a in islice(fields,6,None,2) ]
        genos = zip(a1,a2)

      yield ename,genos

  genos = GenomatrixStream.from_tuples(_load_merlin(),'sdat',loci=loci,genome=genome,phenome=phenome,unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos


class MerlinWriter(object):
  '''
  Object to write Merlin data

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> m = StringIO()
  >>> with MerlinWriter(o,'merlin',genos.loci,genos.genome,genos.phenome,datfile=m) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 0 A A 0 0 C T
  s2 s2 0 0 0 A G C G C C
  s3 s3 0 0 0 G G 0 0 C T
  >>> print m.getvalue() # doctest: +NORMALIZE_WHITESPACE
  M l1
  M l2
  M l3
  '''
  def __init__(self,filename,format,loci,genome,phenome,extra_args=None,**kwargs):
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

    datfile  = get_arg(args, ['datfile','dat'])

    # Careful: mapfile=<blank> is intended to supress output
    if datfile is None:
      datfile = related_file(filename,'dat')

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if phenome is None:
      phenome = Phenome()

    self.out       = autofile(filename,'wb')
    self.loci      = loci
    self.genome    = genome
    self.phenome   = phenome
    self.datfile   = datfile

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

    sex = SEX_RMAP[phenos.sex]
    row = [family or individual,individual,parent1 or '0',parent2 or '0',sex]

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

      sex = SEX_RMAP[phenos.sex]

      row = [family or individual,individual,parent1 or '0',parent2 or '0',sex]

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

    if self.datfile:
      out = autofile(self.datfile,'wb')
      for locus in self.loci:
        out.write('M %s\r\n' % locus)

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


def save_merlin(filename,genos,format,extra_args=None,**kwargs):
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
  >>> save_merlin(o,genos,'merlin')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 s1 0 0 0 A A 0 0 C T
  s2 s2 0 0 0 A G C G C C
  s3 s3 0 0 0 G G 0 0 C T
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_sdat(mergefunc)

  with MerlinWriter(filename, format, genos.loci, genos.genome, genos.phenome,
                                      extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
