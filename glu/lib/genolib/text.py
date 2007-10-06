# -*- coding: utf-8 -*-
'''
File:          text.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU text genotype format input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

from __future__ import with_statement

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import re
import csv

from   itertools                 import islice,dropwhile

from   glu.lib.utils             import tally
from   glu.lib.fileutils         import autofile,namefile

from   glu.lib.genolib.streams   import GenotripleStream, GenomatrixStream
from   glu.lib.genolib.genoarray import model_from_alleles
from   glu.lib.genolib.reprs     import snp,hapmap


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']


class NonUniqueError(ValueError): pass


def unique_check_genomatrixstream(genos):
  '''
  Check that all row and column labels of a genomatrix are unique.  Raises
  a NonUniqueError if they are not.

  @param rows: genotype matrix data with the first row
               being the column meta-data
  @type rows: sequence

  >>> genos = GenomatrixStream([],'sdat',loci=['L1','L2','L3','L1'],models=[])
  >>> unique_check_genomatrixstream(genos)
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique column identifiers: L1:2

  >>> loci=('L1','L2')
  >>> rows=[('R1',['AA','AC']),
  ...       ('R1',['AA','AC'])]
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> genos = unique_check_genomatrixstream(genos)
  >>> list(genos)
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique row identifier: R1
  '''
  dcols  = [ (k,n) for k,n in tally(genos.columns).iteritems() if n>1 ]
  if dcols:
    dcols = ','.join( '%s:%d' % kv for kv in dcols )
    raise NonUniqueError,'Non-unique column identifiers: %s' % dcols

  def _check():
    drows = set()
    for label,row in genos:
      if label in drows:
        raise NonUniqueError,'Non-unique row identifier: %s' % label
      else:
        drows.add(label)

      yield label,row

  return genos.clone(_check(),unique=True)


def load_genomatrix_hapmap(filename,limit=None):
  '''
  Load a HapMap genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of samples loaded. Default is None
  @type         limit: int or None
  @rtype             : GenomatrixStream
  '''
  gfile = autofile(filename)
  gfile = dropwhile(lambda s: s.startswith('#'), gfile)

  try:
    header = gfile.next()
  except StopIteration:
    header = []

  if not any(header.startswith(h) for h in HAPMAP_HEADERS):
    raise ValueError, "Input file '%s' does not appear to be in HapMap format." % namefile(filename)

  if limit is not None:
    limit += 11

  columns = [ intern(h.strip()) for h in islice(header.split(),11,limit) ]
  modelcache = {}
  modelmap   = {}

  def _load():
    n = len(columns)
    for line in gfile:
      fields  = line.split()
      locus   = intern(fields[0].strip())
      alleles = tuple(sorted(fields[1].split('/')))
      genos   = fields[11:limit]
      if len(alleles) != 2 or any(a not in 'ACGT' for a in alleles):
        alleles = tuple(set(a for g in genos for a in g if a!='N'))

      # FIXME: Add error recovery and detection
      assert len(alleles)<=2
      assert len(genos) == n

      model = modelcache.get(alleles)
      if model is None:
        model = modelcache[alleles] = model_from_alleles(alleles,max_alleles=2)
      modelmap[locus] = model

      yield locus,genos

  return GenomatrixStream.from_strings(_load(),'ldat',hapmap,samples=columns,modelmap=modelmap)


def load_genomatrix_linkage(filename,limit=None,modelmap=None):
  '''
  Load a Linkage format genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of samples loaded. Default is None
  @type         limit: int or None
  @rtype             : GenomatrixStream
  '''
  gfile = autofile(filename)

  line,gfile = peekfirst(gfile)

  n    = len(split(line))
  loci = [ '%s' % i for i in range(n-6) ]

  def _load():
    split = re_spaces.split
    amap  = {'0':None}

    for line_num,line in enumerate(gfile):
      fields = split(line)

      if len(fields) != n:
        raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

      sample = '%s_%s' % (field[0],field[1])
      a1s    = islice(fields,6,None,2)
      a2s    = islice(fields,7,None,2)
      genos  = [ (amap.get(a1,a1),amap.get(a2,a2)) for a1,a2 in izip(a1s,a2s) ]

      yield sample,genos

  return GenomatrixStream.from_tuples(_load(),'sdat',loci=loci,modelmap=modelmap)


def load_genomatrix_text(filename,format,genorepr,limit=None,unique=True,modelmap=None):
  '''
  Load the genotype matrix data from file.
  Note that the first row is header and the rest rows are genotypes,
  and the file is tab delimited.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param        limit: limit the number of columms loaded
  @type         limit: int or None
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @param       unique: verify that rows and columns are uniquely labeled
                       (default is True)
  @type        unique: bool
  @param     modelmap: map between a locus and an new internal representation of genotypes. Default is None
  @type      modelmap: dict
  @rtype             : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> genos = load_genomatrix_text(data,'ldat',snp)
  >>> genos.format
  'ldat'
  >>> genos.columns
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  '''
  if genorepr is None:
    raise ValueError('genotype representation must be specified when reading a text format')

  gfile = autofile(filename)
  rows = csv.reader(gfile,dialect='tsv')

  if limit is not None:
    limit += 1

  try:
    columns = iter(rows.next())
  except StopIteration:
    raise ValueError, 'Input file "%s" is empty' % namefile(filename)

  format_found = columns.next()

  # 'key' and blank are always allowed for backward compatibility
  if format_found == 'key':
    format_found = ''

  if format is not None and format_found not in ('',format):
    raise ValueError, 'Input file "%s" does not appear to be in %s format.  Found %s.' \
                        % (namefile(filename),format,format_found)

  columns = tuple(intern(h.strip()) for h in columns)
  format = format_found or format

  def _load(rows):
    n = len(columns)

    # Micro-optimization
    local_intern = intern
    local_strip  = str.strip

    for row in rows:
      label = local_intern(local_strip(row[0]))
      genos = row[1:]

      if len(genos) != n:
        raise ValueError('Invalid genotype matrix row on line %d of %s' % (rows.line_num+1,namefile(filename)))

      yield label,genos

  if format=='ldat':
    genos = GenomatrixStream.from_strings(_load(rows),format,genorepr,samples=columns,modelmap=modelmap)
  else:
    genos = GenomatrixStream.from_strings(_load(rows),format,genorepr,loci=columns,modelmap=modelmap)

  if unique:
    genos = unique_check_genomatrixstream(genos)

  return genos


class TextGenomatrixWriter(object):
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
  >>> with TextGenomatrixWriter(o,genos.format,genos.columns,snp) as w:
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
  def __init__(self,filename,format,header,genorepr,dialect='tsv'):
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
    self.header    = header
    self.headerlen = len(header)
    self.genorepr  = genorepr

    self.out.writerow( [format]+[h.strip() for h in self.header] )

  def writerow(self, rowkey, genos):
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

    if len(genos) != self.headerlen:
      raise ValueError('[ERROR] Internal error: Genotypes do not match header')

    out.writerow( [rowkey]+self.genorepr.to_strings(genos) )

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

    headerlen = self.headerlen
    writerow  = out.writerow
    repr      = self.genorepr.to_strings

    for rowkey,genos in rows:
      if len(genos) != headerlen:
        raise ValueError('[ERROR] Internal error: Genotypes do not match header')
      writerow( [rowkey]+repr(genos) )

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


def save_genomatrix_text(filename,genos,genorepr):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        genos: genomatrix stream
  @type         genos: sequence
  @param     genorepr: function to convert internal genotype representation
                       to the desired string representation
  @type      genorepr: unary function

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> rows = [('s1', [('A','A'),(None,None),('C','T')]),
  ...           ('s2', [('A','G'), ('C','G'), ('C','C')]),
  ...           ('s3', [('G','G'),(None,None),('C','T')]) ]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> save_genomatrix_text(o,genos,snp)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  sdat	l1	l2	l3
  s1	AA	  	CT
  s2	AG	CG	CC
  s3	GG	  	CT
  '''
  with TextGenomatrixWriter(filename, genos.format, genos.columns, genorepr) as writer:
    writer.writerows(genos)


def load_genotriples_text(filename,genorepr,unique=True,limit=None,modelmap=None):
  '''
  Load genotype triples from file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @param       unique: verify that rows and columns are uniquely labeled
                       (default is True)
  @type        unique: bool
  @param        limit: limit the number of genotypes loaded
  @type         limit: int or None
  @param     modelmap: map between a locus and an new internal representation of genotypes. Default is None
  @type      modelmap: dict
  @rtype             : GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO('s1\\tl1\\tAA\\ns1\\tl2\\tGG\\ns2\\tl1\\tAG\\ns2\\tl2\\tCC\\n')
  >>> triples = load_genotriples_text(data,snp)
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  if genorepr is None:
    raise ValueError('genotype representation must be specified when reading a text genotype format')

  rows = csv.reader(autofile(filename),dialect='tsv')

  if limit:
    rows = islice(rows,limit)

  def _load():
    # Micro-optimization
    local_intern = intern
    local_strip  = str.strip
    repr         = genorepr.from_string

    for row in rows:
      if not row:
        continue
      elif len(row) != 3:
        raise ValueError('Invalid genotriple on line %d of %s' % (rows.line_num+1,namefile(filename)))

      sample = local_intern(local_strip(row[0]))
      locus  = local_intern(local_strip(row[1]))
      geno   = repr(row[2])

      yield sample,locus,geno

  return GenotripleStream.from_tuples(_load(),unique=unique,modelmap=modelmap)


class TextGenotripleWriter(object):
  '''
  Object to write genotype triple data to a delimited ASCII file

  Genotype triple files must be supplied as and are output to delimited
  ASCII files as a sequence of three items:

    1. Sample name
    2. Locus name
    3. Genotype

  All rows output have exactly these three columns and no file header is
  output. Sample and locus names are arbitrary and user-specified strings.

  >>> triples = [('s1','l1',('C','T')), ('s1','l2',(None,None)),
  ...            ('s1','l3',('A','A')), ('s2','l2', ('C','C'))]
  >>> triples = iter(GenotripleStream.from_tuples(triples))
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with TextGenotripleWriter(o,snp) as w:
  ...   w.writerow(*triples.next())
  ...   w.writerow(*triples.next())
  ...   w.writerows(triples)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1	l1      CT
  s1	l2
  s1	l3	AA
  s2	l2	CC
  '''
  def __init__(self,filename,genorepr,dialect='tsv'):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    @param      dialect: csv module dialect name ('csv' or 'tsv', default is 'tsv')
    @type       dialect: str or csv dialect object
    '''
    if genorepr is None:
      raise ValueError('genotype representation must be specified when reading a text genotype format')

    self.out       = csv.writer(autofile(filename,'w'),dialect=dialect)
    self.genorepr  = genorepr

  def writerow(self, sample, locus, geno):
    '''
    Write a genotype triple (sample,locus,genotype)

    @param sample: sample identifier
    @type  sample: str
    @param  locus: locus identifier
    @type   locus: str
    @param   geno: genotypes internal representation, to be converted to
                   the appropiate string representation by the supplied
                   genorepr class
    @type    geno: genotype representation
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    out.writerow( (sample,locus,self.genorepr.to_string(geno)) )

  def writerows(self, triples):
    '''
    Write a genotype sequence of triples (sample,locus,genotype)

    @param  triples: sequence of (sample,locus,genotype)
    @type   triples: sequence of (str,str,genotype representation)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    writerow = out.writerow
    repr     = self.genorepr.to_string

    for sample,locus,geno in triples:
      writerow( (sample,locus,repr(geno)) )

  def close(self):
    '''
    Close the writer.

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


def save_genotriples_text(filename,triples,genorepr):
  '''
  Write the genotype triple data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      triples: genotype triple data
  @type       triples: sequence
  @param     genorepr: function to convert internal genotype representation
                       to the desired string representation
  @type      genorepr: unary function

  >>> triples = [ ('s1', 'l1',  ('C','T')),
  ...             ('s1', 'l2', (None,None)),
  ...             ('s1', 'l3',  ('A','A')) ]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> save_genotriples_text(o,triples,snp)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1	l1      CT
  s1	l2
  s1	l3	AA
  '''
  with TextGenotripleWriter(filename,genorepr) as w:
    w.writerows(triples)


def load_genotriples_prettybase(filename,unique=True,limit=None,modelmap=None):
  '''
  Load genotype triples from file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @param       unique: verify that rows and columns are uniquely labeled
                       (default is True)
  @type        unique: bool
  @param        limit: limit the number of genotypes loaded
  @type         limit: int or None
  @param     modelmap: map between a locus and an new internal representation of genotypes. Default is None
  @type      modelmap: dict
  @rtype             : GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO('l1 s1 A A\\nl2 s1 G G\\nl1 s2 N N\\nl2 s2 C C\\n')
  >>> triples = load_genotriples_prettybase(data)
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', (None, None))
  ('s2', 'l2', ('C', 'C'))
  '''
  re_spaces = re.compile('[\t ,]+')
  gfile = autofile(filename)

  if limit:
    gfile = islice(gfile,limit)

  def _load():
    # Micro-optimization
    split        = re_spaces.split
    local_intern = intern
    local_strip  = str.strip
    amap         = {'N':None,'n':None}

    for line_num,line in enumerate(gfile):
      row = split(local_strip(line))
      if not row:
        continue
      elif len(row) != 4:
        raise ValueError('Invalid prettybase row on line %d of %s' % (line_num+1,namefile(filename)))

      locus  = local_intern(local_strip(row[0]))
      sample = local_intern(local_strip(row[1]))
      a1,a2  = row[2],row[3]
      geno   = amap.get(a1,a1),amap.get(a2,a2)

      yield sample,locus,geno

  return GenotripleStream.from_tuples(_load(),unique=unique,modelmap=modelmap)


class PrettybaseGenotripleWriter(object):
  '''
  Object to write genotype triple data to a Prettybase format file

  Genotype triple files must be supplied as and are output to whitespace
  delimited ASCII files as a sequence of four items:

    1. Locus name
    2. Sample name
    3. Allele 1, N for missing
    4. Allele 2, N for missing

  All rows output have exactly these four columns and no file header is
  output. Sample and locus names are arbitrary and user-specified strings.

  >>> triples = [('s1','l1',('C','T')), ('s1','l2',(None,None)),
  ...            ('s1','l3',('A','A')), ('s2','l2', ('C','C'))]
  >>> triples = iter(GenotripleStream.from_tuples(triples))
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with PrettybaseGenotripleWriter(o) as w:
  ...   w.writerow(*triples.next())
  ...   w.writerow(*triples.next())
  ...   w.writerows(triples)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 s1 C T
  l2 s1 N N
  l3 s1 A A
  l2 s2 C C
  '''
  def __init__(self,filename):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    '''
    self.out = autofile(filename,'w')

  def writerow(self, sample, locus, geno):
    '''
    Write a genotype triple (sample,locus,genotype)

    @param sample: sample identifier
    @type  sample: str
    @param  locus: locus identifier
    @type   locus: str
    @param   geno: genotypes internal representation
    @type    geno: genotype representation
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    out.write( ' '.join( [locus,sample,geno[0] or 'N',geno[1] or 'N'] ) )
    out.write('\n')

  def writerows(self, triples):
    '''
    Write a genotype sequence of triples (sample,locus,genotype)

    @param  triples: sequence of (sample,locus,genotype)
    @type   triples: sequence of (str,str,genotype representation)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    write = out.write
    join  = ' '.join

    for sample,locus,geno in triples:
      write( join( [locus,sample,geno[0] or 'N',geno[1] or 'N'] ) )
      write('\n')

  def close(self):
    '''
    Close the writer.

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


def save_genotriples_prettybase(filename,triples):
  '''
  Write the genotype triple data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      triples: genotype triple data
  @type       triples: sequence

  >>> triples = [ ('s1', 'l1',  ('C','T')),
  ...             ('s1', 'l2', (None,None)),
  ...             ('s1', 'l3',  ('A','A')) ]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> save_genotriples_prettybase(o,triples)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 s1 C T
  l2 s1 N N
  l3 s1 A A
  '''
  with PrettybaseGenotripleWriter(filename) as w:
    w.writerows(triples)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
