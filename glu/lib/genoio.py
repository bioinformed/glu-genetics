# -*- coding: utf-8 -*-
'''
File:          genoio.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU genotype data input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

from __future__ import with_statement

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import csv

from   itertools   import izip,islice,dropwhile,chain

from   utils       import tally,deprecated_by
from   fileutils   import autofile,namefile,load_table,guess_format
from   genodata    import GenotripleStream, GenomatrixStream
from   genoarray   import snp_marker


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']

INPUT_FORMATS  = ('ldat','hapmap','sdat','trip','genotriple')
OUTPUT_FORMATS = ('ldat','sdat','trip','genotriple')


class NonUniqueError(ValueError): pass


def unique_check_genomatrixstream(columns,rows):
  '''
  Check that all row and column labels of a genomatrix are unique.  Raises
  a NonUniqueError if they are not.

  @param rows: genotype matrix data with the first row
               being the column meta-data
  @type rows: sequence

  >>> columns,rows = unique_check_genomatrixstream(['L1','L2','L3','L1'],[])
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique column identifiers: L1:2

  >>> columns,rows = unique_check_genomatrixstream(['L1','L2'],[('R1',['AA','AB']),('R1',['AA','AB'])])
  >>> list(rows)
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique row identifier: R1
  '''
  dcols  = [ (k,n) for k,n in tally(columns).iteritems() if n>1 ]
  if dcols:
    dcols = ','.join( '%s:%d' % kv for kv in dcols )
    raise NonUniqueError,'Non-unique column identifiers: %s' % dcols

  def _check():
    drows = set()
    for label,row in rows:
      if label in drows:
        raise NonUniqueError,'Non-unique row identifier: %s' % label
      else:
        drows.add(label)

      yield label,row

  return columns,_check()


def guess_informat(filename):
  return guess_format(filename, INPUT_FORMATS)


def guess_informat_list(filenames):
  formats = set( guess_informat(f) for f in filenames )
  formats.discard(None)
  if len(formats) == 1:
    return formats.pop()
  return None


def guess_outformat(filename):
  return guess_format(filename, OUTPUT_FORMATS)


def load_rename_alleles_file(filename):
  '''
  Load an allele renameing file

  >>> from StringIO import StringIO
  >>> data = StringIO('l1\\tA,C,G,T\\tT,G,C,A\\nl3\\tA\\tC\\nl5\\tA,B\\tC,T')
  >>> for lname,alleles in sorted(load_rename_alleles_file(data).iteritems()):
  ...   print lname,sorted(alleles.iteritems())
  l1 [(None, None), ('A', 'T'), ('C', 'G'), ('G', 'C'), ('T', 'A')]
  l3 [(None, None), ('A', 'C')]
  l5 [(None, None), ('A', 'C'), ('B', 'T')]
  '''
  rows = load_table(filename)

  rename = {}
  for i,row in enumerate(rows):
    if not row:
      continue
    if len(row) != 3:
      raise ValueError('Invalid allele rename record %d in %s' % (i+1,namefile(filename)))

    lname,old_alleles,new_alleles = row

    lname       = intern(lname.strip())
    old_alleles = [ intern(a.strip()) for a in old_alleles.split(',') ]
    new_alleles = [ intern(a.strip()) for a in new_alleles.split(',') ]

    if len(old_alleles) != len(new_alleles):
      raise ValueError('Invalid allele rename record %d in %s' % (i+1,namefile(filename)))

    locus_rename = dict( izip(old_alleles,new_alleles) )
    locus_rename[None] = None

    if lname in rename and rename[lname] != locus_rename:
      raise ValueError('Inconsistent rename record %d in %s' % (i+1,namefile(filename)))

    rename[lname] = locus_rename

  return rename


def load_hapmap_genotypes(filename,limit=None,genorepr=snp_marker.pack_strs):
  '''
  Load the hampmap genotype data from file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of samples loaded
  @type         limit: int or None
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @return            : rows of genotypes with the first row being the sample names
  @rtype             : generator
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

  header = list(islice(header.split(),11,limit))
  header = [ intern(h.strip()) for h in header ]

  def _load_hapmap_genotypes():
    n = len(header)
    for line in gfile:
      fields = line.split()
      genos = genorepr(g.replace('N','') for g in islice(fields,11,limit))
      assert len(genos) == n
      yield intern(fields[0].strip()),genos

  return header,_load_hapmap_genotypes()


def load_genomatrix(filename,format=None,limit=None,genorepr=snp_marker.pack_strs,unique=True):
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
  @return            : format and sequence of column names followed by
                       tuples of row label and row data
  @rtype             : tuple of string and generator

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> format,columns,rows = load_genomatrix(data,'ldat')
  >>> format
  'ldat'
  >>> columns
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  '''
  format,columns,rows = _load_genomatrix(filename,format=format,limit=limit,genorepr=genorepr)

  if unique:
    columns,rows = unique_check_genomatrixstream(columns,rows)

  return format,columns,rows


def _load_genomatrix(filename,format,limit,genorepr):
  '''See load_genomatrix'''

  gfile = autofile(filename)
  rows = csv.reader(gfile,dialect='tsv')

  if limit is not None:
    limit += 1

  try:
    header = iter(rows.next())
  except StopIteration:
    raise ValueError, 'Input file "%s" is empty' % namefile(filename)

  format_found = header.next()

  # 'key' and blank are always allowed for backward compatibility
  if format_found == 'key':
    format_found = ''

  if format is not None and format_found not in ('',format):
    raise ValueError, 'Input file "%s" does not appear to be in %s format.  Found %s.' \
                        % (namefile(filename),format,format_found)

  header = tuple(intern(h.strip()) for h in header)

  def _gen_load_genomatrix():
    n = len(header)

    # Micro-optimization
    ir = iter
    it = intern
    st = str.strip
    gr = genorepr

    for i,row in enumerate(rows):
      row   = ir(row)
      label = it(st(row.next()))
      genos = gr(row)

      if len(genos) != n:
        raise ValueError('Invalid genotype matrix row on line %d of %s' % (i+1,namefile(filename)))

      yield label,genos

  return format_found or format,header,_gen_load_genomatrix()


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

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with TextGenomatrixWriter(o,'ldat',('l1','l2','l3'),snp_marker) as w:
  ...   w.writerow('s1', [('A','A'),  None,   ('C','T')])
  ...   w.writerow('s2', [('A','G'),('C','G'),('C','C')])
  ...   w.writerows( [('s3', [('G','G'),  None,   ('C','T')]),
  ...                 ('s4', [('A','A'),  None,     None   ])])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  ldat	l1	l2	l3
  s1	AA	  	CT
  s2	AG	CG	CC
  s3	GG	  	CT
  s4	AA
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
    self.out       = csv.writer(autofile(filename,'w'),dialect=dialect)
    self.header    = header
    self.headerlen = len(header)
    self.genorepr  = genorepr

    self.out.writerow( [format]+list(self.header) )

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

    out.writerow( [rowkey]+self.genorepr.strs_from_reps(genos) )

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
    repr      = self.genorepr.strs_from_reps

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

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with TextGenotripleWriter(o,snp_marker) as w:
  ...   w.writerow('s1', 'l1', ('C','T'))
  ...   w.writerow('s1', 'l2',   None   )
  ...   w.writerows( [('s1', 'l3', ('A','A')),
  ...                 ('s2', 'l2', ('C','C'))])
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

    out.writerow( (sample,locus,self.genorepr.str_from_rep(geno)) )

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
    repr     = self.genorepr.str_from_rep

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


def save_genomatrix(filename,header,matrix,format=None,genorepr=snp_marker.strs_from_reps):
  '''
  Write the genotype matrix data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       matrix: genotype matrix data
  @type        matrix: sequence
  @param       format: text string output in the first header field to
                       indicate data format (default is blank)
  @type        format: string
  @param     genorepr: function to convert internal genotype representation
                       to the desired string representation
  @type      genorepr: unary function

  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> loci =              ('l1',     'l2',    'l3')
  >>> matrix = [('s1', [('A','A'),  None,   ('C','T')]),
  ...           ('s2', [('A','G'),('C','G'),('C','C')]),
  ...           ('s3', [('G','G'),  None,   ('C','T')]) ]
  >>> save_genomatrix(o,loci,matrix,'ldat',snp_marker.strs_from_reps)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  ldat	l1	l2	l3
  s1	AA	  	CT
  s2	AG	CG	CC
  s3	GG	  	CT
  '''
  # Create temporary representation class -- ugly, but necessary for
  # backward compatibility
  class TmpRepr(object):
    strs_from_reps = genorepr

  with TextGenomatrixWriter(filename, format or '', header, TmpRepr()) as writer:
    writer.writerows(matrix)


def load_genotriples(filename,limit=None,genorepr=snp_marker.rep_from_str):
  '''
  Load genotype triples from file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of genotypes loaded
  @type         limit: int or None
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @return            : sequence of tuples of sample name, locus name, and genotype representation
  @rtype             : generator

  >>> from StringIO import StringIO
  >>> data = StringIO('s1\\tl1\\tAA\\ns1\\tl2\\tGG\\ns2\\tl1\\tAG\\ns2\\tl2\\tCC\\n')
  >>> triples = load_genotriples(data,genorepr=snp_marker.rep_from_str)
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  rows = csv.reader(autofile(filename),dialect='tsv')

  if limit:
    rows = islice(rows,limit)

  # Micro-optimization
  it = intern
  st = str.strip

  for i,row in enumerate(rows):
    if not row:
      continue
    elif len(row) != 3:
      raise ValueError('Invalid genotriple on line %d of %s' % (i+1,namefile(filename)))

    sample = it(st(row[0]))
    locus  = it(st(row[1]))
    geno   = genorepr(row[2])

    yield sample,locus,geno


def save_genotriples(filename,triples,genorepr=snp_marker.str_from_rep):
  '''
  Write the genotype triple data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      triples: genotype triple data
  @type       triples: sequence
  @param     genorepr: function to convert internal genotype representation
                       to the desired string representation
  @type      genorepr: unary function

  >>> triples = [ ('s1', 'l1', ('C','T')),
  ...             ('s1', 'l2',   None   ),
  ...             ('s1', 'l3', ('A','A')) ]
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> save_genotriples(o,triples,snp_marker.str_from_rep)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1	l1      CT
  s1	l2
  s1	l3	AA
  '''
  # Create temporary representation class -- ugly, but necessary for
  # backward compatibility
  class TmpRepr(object):
    str_from_rep = genorepr

  with TextGenotripleWriter(filename,TmpRepr()) as w:
    w.writerows(triples)


def load_genotriplestream(filename, limit=None, unique=False):
  '''
  Load genotriple file and return a GenotripleStream object

  @param filename: file name or file object
  @type  filename: str or file object
  @param    limit: limit the number of samples loaded
  @type     limit: int or None
  @param genorepr: object representing the input/output encoding and
                   internal representation of genotypes
  @type  genorepr: UnphasedMarkerRepresentation or similar object
  @param   unique: flag indicating if repeated triplets do not exist
  @type    unique: bool
  @return        : genotriple stream
  @rtype         : GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO('s1\\tl1\\tAA\\ns1\\tl2\\tGG\\ns2\\tl1\\tAG\\ns2\\tl2\\tCC\\n')
  >>> triples = load_genotriplestream(data)
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  triples = load_genotriples(filename,limit=None)
  return GenotripleStream(triples, unique=unique)


def load_genomatrixstream(filename, format=None, limit=None, unique=True):
  '''
  Load genomatrix file depending on matrix format and return a GenotripleMatrix object

  @param filename: file name or file object
  @type  filename: str or file object
  @param   format: format of input genomatrix, 'hapmap', 'ldat' or 'sdat'
  @type    format: str
  @param    limit: limit the number of samples loaded
  @type     limit: int or None
  @param genorepr: object representing the input/output encoding and
                   internal representation of genotypes
  @type  genorepr: UnphasedMarkerRepresentation or similar object
  @param   unique: flag indicating if repeated row or column elements do not exist
  @type    unique: bool
  @return        : loaded genomatrix stream
  @rtype         : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> ldat = load_genomatrixstream(data,'ldat')
  >>> ldat.columns
  ('s1', 's2', 's3')
  >>> for row in ldat:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  '''
  if format is None:
    format = guess_informat(filename)

  samples = loci = None

  if format == 'hapmap':
    samples,genos = load_hapmap_genotypes(filename,limit=limit)
    format = 'ldat'
  elif format == 'ldat':
    format,samples,genos = load_genomatrix(filename,format=format,limit=limit,unique=unique)
  elif format == 'sdat':
    format,loci,genos = load_genomatrix(filename,format=format,limit=limit,unique=unique)
  elif not format:
    raise ValueError, "Input file format for '%s' must be specified" % namefile(filename)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % format

  return GenomatrixStream(genos,format,samples=samples,loci=loci,unique=unique,packed=True)


def load_genostream(filename, format=None, limit=None, unique=None):
  '''
  Load genotype data in the format of (ldat, sdat, hapmap, trip, genotriple) and return a
  GenomatrixStream or GenotripleStream object

  @param filename: file name or file object
  @type  filename: str or file object
  @param   format: format of input
  @type    format: str
  @param    limit: limit the number of samples loaded
  @type     limit: int or None
  @param genorepr: function to convert a list of genotype strings to desired
                   internal representation
  @type  genorepr: unary function
  @param   unique: flag indicating if repeated row or column elements do not exist
  @type    unique: bool
  @return        : new genomatrix stream
  @rtype         : GenomatrixStream or GenotripleStream

  >>> from StringIO import StringIO
  >>> loaded = load_genostream(StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n"),'ldat')
  >>> loaded.columns
  ('s1', 's2', 's3')
  >>> for row in loaded:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  >>> loaded.samples
  ('s1', 's2', 's3')
  >>> loaded.loci
  >>> loaded.unique
  True
  '''
  if format is None:
    format = guess_informat(filename)

  if format in ('ldat','sdat','hapmap'):
    if unique is None:
      unique = True
    return load_genomatrixstream(filename,format=format,limit=limit,unique=unique)

  elif format in ('trip','genotriple'):
    if unique is None:
      unique = False
    return load_genotriplestream(filename, limit=limit, unique=unique)

  elif not format:
    raise ValueError, "Input file format for '%s' must be specified" % namefile(filename)

  else:
    raise NotImplementedError,"File format '%s' is not supported" % format


def save_genostream(filename, genos, format=None, mergefunc=None):
  '''
  Write genotype data to file in one of the specified formats (ldat, sdat, trip, genotriple).

  @param  filename: file name or file object
  @type   filename: str or file object
  @param     genos: genomatrix/genotriple
  @type      genos: genomatrix/genotriple generator
  @param    format: format of input
  @type     format: str
  @param mergefunc: optional function to merge multiple genotypes into a consensus genotype
  @type  mergefunc: function or None
  '''
  if format is None:
    format = guess_outformat(filename)

  if mergefunc is not None:
    genos = genos.merged(mergefunc)

  if format == 'ldat':
    genos = genos.as_ldat(mergefunc)
    save_genomatrix(filename, genos.columns, genos, format)
  elif format == 'sdat':
    genos = genos.as_sdat(mergefunc)
    save_genomatrix(filename, genos.columns, genos, format)
  elif format in ('trip','genotriple'):
    save_genotriples(filename, genos.as_genotriples())
  elif not format:
    raise ValueError, "Output file format for '%s' must be specified" % namefile(filename)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % format


def transform_files(infiles,informat,ingenorepr,
                    outfile,outformat,outgenorepr,
                    transform=None,
                    mergefunc=None,limit=None):
  '''
  The driver for transforming multiple genodata files into different formats
  (ldat, sdat, trip, or genotriples), representations (snp_acgt, snp_ab,
  snp_marker, generic_marker) and, depending on the presence and attributes
  of the transform object, performing operations on samples and loci such as
  exclude, include, and rename.

  @param     infiles: list of input file names or file objects
  @type      infiles: str or file objects
  @param    informat: input file format for all input files
  @type     informat: str
  @param  ingenorepr: internal genotype representation for the input
  @type   ingenorepr: UnphasedMarkerRepresentation or similar object
  @param    outfiles: output file name or file object
  @type     outfiles: str or file object
  @param   outformat: output file format
  @type    outformat: str
  @param outgenorepr: internal genotype representation for the output
  @type  outgenorepr: UnphasedMarkerRepresentation or similar object
  @param   transform: transformation object (optional)
  @type    transform: GenoTransform object
  @param       limit: limit the number of samples loaded
  @type        limit: int or None
  @return           : transformed genodata
  @rtype            : GenomatrixStream or GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> out  = StringIO()
  >>> transform_files([data],'ldat',snp_marker,out,'trip',snp_marker)
  >>> print out.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1	l1	AA
  s2	l1	AG
  s3	l1	GG
  s1	l2	CC
  s2	l2	CT
  s3	l2	TT
  '''
  if informat is None:
    informat = guess_informat_list(infiles)

  genos = [ load_genostream(f,informat,limit=limit).transformed(transform) for f in infiles ]
  n = len(genos)

  if outformat is None:
    outformat = guess_outformat(outfile)

  # Guess output format based on input format if it is unique
  if outformat is None:
    outformat = informat

  if outformat == 'ldat':
    genos = GenomatrixStream.from_streams(genos,'ldat',mergefunc=mergefunc)
  elif outformat == 'sdat':
    genos = GenomatrixStream.from_streams(genos,'sdat',mergefunc=mergefunc)
  elif outformat in ('trip','genotriple'):
    genos = GenotripleStream.from_streams(genos,mergefunc=mergefunc)
  elif not outformat:
    raise ValueError, "Output file format for '%s' must be specified" % namefile(outfile)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % outformat

  # Order again after merging, if necessary
  if n>1 and (transform.loci.order or transform.samples.order):
    genos = genos.transformed(order_loci=transform.loci.order,
                              order_samples=transform.samples.order)

  save_genostream(outfile,genos,outformat)


def _test_genoio():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test_genoio()
