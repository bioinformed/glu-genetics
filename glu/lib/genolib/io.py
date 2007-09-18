# -*- coding: utf-8 -*-
'''
File:          io.py

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

from   operator          import getitem, itemgetter
from   collections       import defaultdict
from   itertools         import izip,islice,dropwhile,imap,repeat

from   glu.lib.utils     import tally
from   glu.lib.fileutils import autofile,namefile,guess_format

from   streams           import GenotripleStream, GenomatrixStream
from   genoarray         import model_from_alleles
from   reprs             import snp,hapmap,marker
from   text              import TextGenomatrixWriter,    TextGenotripleWriter,    \
                                save_genotriples_text,   load_genotriples_text,   \
                                save_genomatrix_text,    load_genomatrix_text,    \
                                load_genomatrix_hapmap
from   binary            import BinaryGenomatrixWriter,  BinaryGenotripleWriter,  \
                                save_genotriples_binary, load_genotriples_binary, \
                                save_genomatrix_binary,  load_genomatrix_binary


INPUT_FORMATS  = ('ldat','hapmap','sdat','trip','genotriple','lbat','sbat','tbat')
OUTPUT_FORMATS = ('ldat','sdat','trip','genotriple','lbat','sbat','tbat')


def guess_informat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''

  return guess_format(filename, INPUT_FORMATS)


def guess_informat_list(filenames):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''

  formats = set( guess_informat(f) for f in filenames )
  formats.discard(None)
  if len(formats) == 1:
    return formats.pop()
  return None


def guess_outformat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''

  return guess_format(filename, OUTPUT_FORMATS)


def load_genostream(filename, format=None, genorepr=None, limit=None, unique=True, modelmap=None):
  '''
  Load genomatrix file depending on matrix format and return a GenotripleMatrix object

  @param filename: a file name or file object
  @type  filename: str or file object
  @param   format: format of input file , 'hapmap', 'ldat', 'sdat', 'trip', 'genotriple'
                   'lbat', 'sbat', or 'tbat'. Default is None
  @type    format: str
  @param genorepr: internal representation of genotypes for the input/output. Default is None
  @type  genorepr: UnphasedMarkerRepresentation or similar object
  @param    limit: limit the number of samples loaded. Default is None
  @type     limit: int or None
  @param   unique: flag indicating if repeated row or column elements do not exist. Default is None
  @type    unique: bool
  @para  modelmap: map between a locus and an new internal representation of genotypes. Default is None
  @type  modelmap: dict
  @return        : loaded genomatrix stream
  @rtype         : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> ldat = load_genostream(data,'ldat',snp)
  >>> ldat.columns
  ('s1', 's2', 's3')
  >>> for row in ldat:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  >>> ldat.loci
  >>> ldat.unique
  True

  >>> from StringIO import StringIO
  >>> data = StringIO('s1\\tl1\\tAA\\ns1\\tl2\\tGG\\ns2\\tl1\\tAG\\ns2\\tl2\\tCC\\n')
  >>> triples = load_genostream(data,'trip',genorepr=snp)
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  if format is None:
    format = guess_informat(filename)

  samples = loci = None

  if format == 'hapmap':
    genos = load_genomatrix_hapmap(filename,limit=limit)
  elif format == 'ldat':
    genos = load_genomatrix_text(filename,format,genorepr,limit=limit,unique=unique,modelmap=modelmap)
  elif format == 'sdat':
    genos = load_genomatrix_text(filename,format,genorepr,limit=limit,unique=unique,modelmap=modelmap)
  elif format == 'lbat':
    genos = load_genomatrix_binary(filename,'ldat',limit=limit,unique=unique,modelmap=modelmap)
  elif format == 'sbat':
    genos = load_genomatrix_binary(filename,'sdat',limit=limit,unique=unique,modelmap=modelmap)
  elif format in ('trip','genotriple'):
    genos = load_genotriples_text(filename,genorepr,limit=limit,unique=unique,modelmap=modelmap)
  elif format=='tbat':
    genos = load_genotriples_binary(filename,limit=limit,unique=unique,modelmap=modelmap)
  elif not format:
    raise ValueError, "Input file format for '%s' must be specified" % namefile(filename)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % format

  return genos


def save_genostream(filename, genos, format=None, genorepr=None, mergefunc=None, compress=True):
  '''
  Write genotype data to file in one of the specified formats (ldat, sdat, trip, genotriple).

  @param  filename: a file name or file object
  @type   filename: str or file object
  @param     genos: genomatrix/genotriple stream
  @type      genos: sequence
  @param    format: format of input file , 'ldat', 'sdat', 'trip', 'genotriple'
                    'lbat', 'sbat', or 'tbat'. Default is None
  @type     format: str
  @param  genorepr: internal representation of genotypes for the input/output. Default is None
  @type   genorepr: UnphasedMarkerRepresentation or similar object
  @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type  mergefunc: callable
  @param  compress: flag indicating if a compressed format is desired. Default is True
  @type   compress: bool
  '''
  if format is None:
    format = guess_outformat(filename)

  if mergefunc is not None:
    genos = genos.merged(mergefunc)

  if format == 'ldat':
    save_genomatrix_text(filename, genos.as_ldat(mergefunc), genorepr)
  elif format == 'sdat':
    save_genomatrix_text(filename, genos.as_sdat(mergefunc), genorepr)
  elif format in ('trip','genotriple'):
    save_genotriples_text(filename, genos.as_genotriples(), genorepr)
  elif format == 'lbat':
    save_genomatrix_binary(filename, genos.as_ldat(mergefunc), compress=compress)
  elif format == 'sbat':
    save_genomatrix_binary(filename, genos.as_sdat(mergefunc), compress=compress)
  elif format == 'tbat':
    save_genotriples_binary(filename, genos.as_genotriples(), compress=compress)
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
  (ldat, sdat, trip, or genotriples), representations (...) and, depending
  on the presence and attributes of the transform object, performing
  operations on samples and loci such as exclude, include, and rename.

  @param     infiles: list of input file names or file objects
  @type      infiles: str or file objects
  @param    informat: format of input file, 'hapmap', 'ldat', 'sdat', 'trip', 'genotriple'
                      'lbat', 'sbat', or 'tbat'
  @type     informat: str
  @param  ingenorepr: internal genotype representation for the input
  @type   ingenorepr: UnphasedMarkerRepresentation or similar object
  @param    outfiles: output file name or file object
  @type     outfiles: str or file object
  @param   outformat: format of output file , 'ldat', 'sdat', 'trip', 'genotriple'
                      'lbat', 'sbat', or 'tbat'
  @type    outformat: str
  @param outgenorepr: internal genotype representation for the output
  @type  outgenorepr: UnphasedMarkerRepresentation or similar object
  @param   transform: transformation object (optional)
  @type    transform: GenoTransform object
  @param   mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type    mergefunc: callable
  @param       limit: limit the number of samples loaded
  @type        limit: int or None
  @return           : transformed genodata
  @rtype            : GenomatrixStream or GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\t\\tCT\\tTT\\n")
  >>> out  = StringIO()
  >>> transform_files([data],'ldat',snp,out,'trip',marker)
  >>> print out.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1  l1      A/A
  s2  l1      A/G
  s3  l1      G/G
  s1  l2
  s2  l2      C/T
  s3  l2      T/T
  '''
  if informat is None:
    informat = guess_informat_list(infiles)

  genos = [ load_genostream(f,informat,ingenorepr,limit=limit).transformed(transform) for f in infiles ]
  n = len(genos)

  if outformat is None:
    outformat = guess_outformat(outfile)

  # Guess output format based on input format if it is unique
  if outformat is None:
    outformat = informat

  if outformat in ('ldat','lbat'):
    genos = GenomatrixStream.from_streams(genos,'ldat',mergefunc=mergefunc)
  elif outformat in ('sdat','sbat'):
    genos = GenomatrixStream.from_streams(genos,'sdat',mergefunc=mergefunc)
  elif outformat in ('trip','genotriple','tbat'):
    genos = GenotripleStream.from_streams(genos,mergefunc=mergefunc)
  elif not outformat:
    raise ValueError, "Output file format for '%s' must be specified" % namefile(outfile)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % outformat

  # Order again after merging, if necessary
  if n>1 and (transform.loci.order or transform.samples.order):
    genos = genos.transformed(order_loci=transform.loci.order,
                              order_samples=transform.samples.order)

  save_genostream(outfile,genos,outformat,outgenorepr)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
