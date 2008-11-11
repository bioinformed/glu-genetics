# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'EIGENSOFT (smartpca) genotype format input/output'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   itertools                 import izip

from   glu.lib.utils             import izip_exact
from   glu.lib.fileutils         import autofile,namefile,parse_augmented_filename,get_arg, \
                                        related_file,guess_related_file

from   glu.lib.genolib.streams   import GenotripleStream,GenomatrixStream,NonUniqueError
from   glu.lib.genolib.genoarray import count_genotypes, count_alleles_from_genocounts, \
                                        major_allele_from_allelecounts, \
                                        GenotypeArrayDescriptor, GenotypeArray, UnphasedMarkerModel
from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN


__all__ = ['load_eigensoft_smartpca','save_eigensoft_smartpca','EigensoftSmartPCAWriter']


SEX_MAP  = {'M':SEX_MALE,'m':SEX_MALE,'F':SEX_FEMALE,'f':SEX_FEMALE,'U':SEX_UNKNOWN}
SEX_RMAP = {SEX_UNKNOWN:'U', SEX_MALE:'M', SEX_FEMALE:'F'}

CHR_MAP  = {'0':None,'23':'X','24':'Y','90':'M','91':'XY'}
CHR_MAP.update( (str(i),)*2 for i in range(1,23) )

CHR_RMAP = {None:'1','X':'23','Y':'24','M':'90','XY':'91'}
CHR_RMAP.update( (str(i),)*2 for i in range(1,23) )

DEFAULT_ALLELES = 'A','B'
UNKNOWN = 'UNKNOWN'


def load_eigensoft_snps(filename,genome):
  loci = []
  models = []
  modelcache = {}
  for i,line in enumerate(autofile(filename)):
    line = line.strip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    m = len(fields)

    lname = fields[0]

    if not lname:
      raise ValueError('Invalid SMARTPCA locus record %d' % (i+1))

    chr = location = None
    if m>1:
      chr = CHR_MAP[fields[1]]
    # m==2 is genetic map location
    if m>3:
      location = int(fields[3]) or None

    if m>=6:
      alleles = (intern(fields[4]),intern(fields[5]))
    else:
      alleles = DEFAULT_ALLELES

    model = modelcache.get(alleles)
    if model is None:
      model = modelcache[alleles] = UnphasedMarkerModel(max_alleles=2)

      # Assign genotypes in the appropriate order
      a,b = alleles
      model.add_genotype( (a,a) )
      model.add_genotype( (a,b) )
      model.add_genotype( (b,b) )

    loci.append(lname)
    models.append(model)
    genome.set_locus(lname,model=model,chromosome=chr,location=location)

  return loci,models


def load_eigensoft_inds(filename,phenome):
  samples = []

  for i,line in enumerate(autofile(filename)):
    line = line.strip()

    if not line or line.startswith('#'):
      continue

    fields = line.split()

    m = len(fields)

    if m not in (2,3):
      raise ValueError('Invalid SMARTPCA ind record %d' % (i+1))

    name = fields[0]

    if not name:
      raise ValueError('Invalid SMARTPCA ind record %d' % (i+1))

    sex = SEX_MAP[fields[1]]
    phenoclass = intern(fields[2]) if m==3 else UNKNOWN

    samples.append(name)
    phenome.merge_phenos(name, sex=sex, phenoclass=phenoclass)

  return samples


def load_eigensoft_smartpca(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a SMARTPCA format genotype data.

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

  >>> from StringIO import StringIO
  >>> genos = StringIO('012\\n919\\n101\\n')
  >>> snp   = StringIO('l1 1 0 0 G A\\nl2 1 0 0 G C\\nl3 1 0 0 C T\\n')
  >>> ind   = StringIO('s1 U\\ns2 U\\ns3 U\\n')
  >>> genos = load_eigensoft_smartpca(genos,'eigensoft',ind=ind,snp=snp)
  >>> genos.format
  'ldat'
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), ('A', 'G'), ('A', 'A')])
  ('l2', [(None, None), ('C', 'G'), (None, None)])
  ('l3', [('C', 'T'), ('C', 'C'), ('C', 'T')])
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  unique = get_arg(args, ['unique'], True)
  loci   = get_arg(args, ['snp']) or guess_related_file(filename,['snp'])
  ind    = get_arg(args, ['ind']) or guess_related_file(filename,['ind'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if not loci:
    raise ValueError('Eigenstrat/SmartPCA loader requires a SNP file')

  if not ind:
    raise ValueError('Eigenstrat/SmartPCA loader requires an IND file')

  if phenome is None:
    phenome = Phenome()

  file_genome = Genome()
  loci,models = load_eigensoft_snps(loci,file_genome)
  samples     = load_eigensoft_inds(ind,phenome)

  is_unique = len(set(loci)) == len(loci) and len(set(samples)) == len(samples)
  if unique and not is_unique:
    raise NonUniqueError('EIGENSOFT/SMARTPCA data contains non-unique snps or individuals')

  rows = autofile(filename)
  def _load():
    n = len(samples)

    descrcache = {}
    for locus,model,row in izip_exact(loci,models,rows):
      row = row.rstrip()
      if len(row) != n:
        raise ValueError('Invalid genotype row on line %d of %s' % (rows.line_num+1,namefile(filename)))

      # FIXME: This can be done in load_smartpca_snps
      dg = descrcache.get(model)
      if dg is None:
        gmap  = { '9' : model.genotypes[0],
                  '0' : model.genotypes[1],
                  '1' : model.genotypes[2],
                  '2' : model.genotypes[3] }
        descr = GenotypeArrayDescriptor([model]*n)
        descrcache[model] = descr,gmap
      else:
        descr,gmap = dg

      genos = GenotypeArray(descr, (gmap[g] for g in row))

      yield locus,genos

  genos = GenomatrixStream(_load(),'ldat',samples=samples,loci=loci,models=models,
                                          genome=file_genome,phenome=phenome,
                                          unique=unique,packed=True)

  if genome:
    genos = genos.transformed(recode_models=genome)

  return genos


class EigensoftSmartPCAWriter(object):
  '''
  Object to write Eigensoft SMARTPCA genotype data

  >>> loci =           ('l1',      'l2',      'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_ldat()
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> s = StringIO()
  >>> i = StringIO()
  >>> with EigensoftSmartPCAWriter(o,'eigensoft',genos.samples,genos.genome,genos.phenome,
  ...                                ind=i,snp=s) as w:
  ...   genos=iter(genos)
  ...   w.writerow(*genos.next())
  ...   w.writerow(*genos.next())
  ...   w.writerows(genos)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  210
  919
  101
  >>> print s.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 1 0 0 G A
  l2 1 0 0 G C
  l3 1 0 0 C T
  >>> print i.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 U UNKNOWN
  s2 U UNKNOWN
  s3 U UNKNOWN
  '''
  def __init__(self,filename,format,samples,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param       format: data format string
    @type        format: str
    @param       header: column headings
    @type        header: list or str
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    snpfile = get_arg(args, ['snp'])
    indfile = get_arg(args, ['ind'])

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    # Careful: file=<blank> is intended to supress output
    if snpfile is None:
      snpfile = related_file(filename,'snp')
    if indfile is None:
      indfile = related_file(filename,'ind')

    self.samples = samples
    self.genome  = genome
    self.phenome = phenome
    self.out     = autofile(filename,'w')
    self.snpout  = autofile(snpfile,'w') if snpfile else None

    if indfile:
      indout = autofile(indfile,'w')
      for sample in samples:
        phenos = self.phenome.get_phenos(sample)
        sex = SEX_RMAP[phenos.sex]
        phenoclass = phenos.phenoclass or ''
        indout.write(' '.join([sample,sex,phenoclass or UNKNOWN]))
        indout.write('\n')

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

    model        = genos[0].model
    genocounts   = count_genotypes(genos)
    allelecounts = count_alleles_from_genocounts(model,genocounts)

    # FIXME: Finding major seems much more robust, since it is better defined
    try:
      major,freq = major_allele_from_allelecounts(model,allelecounts)
    except ValueError:
      major = None

    # FIXME: Uninformative locus, so emit a warning
    if major is None:
      return

    other = [ a for a,n in izip(model.alleles[1:],allelecounts[1:]) if a!=major and n ]

    # Non-biallelic locus
    if len(other) > 1:
      # FIXME: Non-bi-alleleic locus, so emit a warning
      return
    if len(other) == 1:
      other = other[0]
    else:
      # Monomorphic, pick a default representation that does not conflict with major
      other = [ a for a,n in izip(model.alleles[1:],allelecounts[1:]) if a!=major ]
      if other:
        other = other[0]
      elif major!='X':
        other = 'X'
      else:
        other = 'N'

    gmap = { model[major,major]:'0',
             model[other,major]:'1',
             model[other,other]:'2',
             model[None,None]  :'9' }

    row = ''.join( gmap[g] for g in genos )

    snpout = self.snpout
    if snpout:
      loc = self.genome.get_locus(locus)
      chr = CHR_RMAP[loc.chromosome]
      pos = str(loc.location or '0')

      snpout.write(' '.join([locus,chr,'0',pos,major,other]))
      snpout.write('\n')

    out.write(row)
    out.write('\n')

  def writerows(self, rows):
    '''
    Write rows of genotypes given pairs of row key and list of genotypes

    @param rows: sequence of pairs of row key and sequence of genotypes in
                 an internal representation
                 class.
    @type  rows: sequence of (str,sequence)
    '''
    for locus,genos in rows:
      self.writerow(locus,genos)

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


def save_eigensoft_smartpca(filename,genos,format,extra_args=None,**kwargs):
  '''
  Write Eigensoft SMARTPCA genotype data

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> i = StringIO()
  >>> s = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...         ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...         ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> save_eigensoft_smartpca(o,genos,'eigensoft',ind=i,snp=s)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  210
  919
  101
  >>> print s.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 1 0 0 G A
  l2 1 0 0 G C
  l3 1 0 0 C T
  >>> print i.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1 U UNKNOWN
  s2 U UNKNOWN
  s3 U UNKNOWN
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)
  mergefunc = get_arg(args, ['mergefunc'])

  genos = genos.as_ldat(mergefunc)

  with EigensoftSmartPCAWriter(filename,format,genos.samples,genos.genome,genos.phenome,
                                        extra_args=args) as writer:

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    writer.writerows(genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
