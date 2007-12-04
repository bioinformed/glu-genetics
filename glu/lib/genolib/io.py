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


from   glu.lib.fileutils         import namefile, guess_format, parse_augmented_filename, get_arg

from   glu.lib.genolib.streams   import GenotripleStream, GenomatrixStream
from   glu.lib.genolib.merge     import get_genomerger
from   glu.lib.genolib.locus     import load_genome, Genome
from   glu.lib.genolib.reprs     import get_genorepr

# FIXME: Format support should ultimately be pluggable with a registration protocol
from   glu.lib.genolib.formats   import *


INPUT_FORMATS  = ('ldat','hapmap','sdat','tdat','trip','genotriple','prettybase','pb','lbat','sbat','tbat')
OUTPUT_FORMATS = ('ldat','sdat','tdat','trip','genotriple','prettybase','pb','lbat','sbat','tbat')


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


def load_genostream(filename, extra_args=None, **kwargs):
  '''
  Load genomatrix file depending on matrix format and return a GenotripleMatrix object

  @param filename: a file name or file object
  @type  filename: str or file object
  @param   format: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                   lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type    format: str
  @param  genorepr: string or representation object for text genotypes. Default is None
  @type   genorepr: str, UnphasedMarkerRepresentation or similar object
  @param    limit: limit the number of samples loaded. Default is None
  @type     limit: int or None
  @param   unique: flag indicating if repeated row or column elements do not exist. Default is None
  @type    unique: bool
  @para    genome: map between a locus and an new internal genotype
                   representation. If a string is specified, it is passed to load_genome().
                   Default is None.
  @type    genome: str or Genome instance
  @return        : loaded genomatrix stream
  @rtype         : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> ldat = load_genostream(data,format='ldat',genorepr='snp')
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
  >>> triples = load_genostream(data,format='tdat',genorepr='snp')
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  format   = get_arg(args, ['format'])
  genome   = get_arg(args, ['genome','loci'])
  genorepr = get_arg(args, ['genorepr']) or 'snp'
  unique   = get_arg(args, ['unique'], True)
  limit    = int(get_arg(args, ['limit']) or 0) or None

  if format is None:
    format = guess_informat(filename)

  if isinstance(genorepr,basestring):
    genorepr = get_genorepr(genorepr)

  if isinstance(genome,basestring) or genome is None:
    genome = load_genome(genome,extra_args=args)

  if format == 'hapmap':
    genos = load_genomatrix_hapmap(filename,limit=limit,genome=genome)
  elif format == 'ldat':
    genos = load_genomatrix_text(filename,format,genorepr,limit=limit,unique=unique,genome=genome)
  elif format == 'sdat':
    genos = load_genomatrix_text(filename,format,genorepr,limit=limit,unique=unique,genome=genome)
  elif format == 'lbat':
    genos = load_genomatrix_binary(filename,'ldat',limit=limit,unique=unique,genome=genome)
  elif format == 'sbat':
    genos = load_genomatrix_binary(filename,'sdat',limit=limit,unique=unique,genome=genome)
  elif format in ('tdat','trip','genotriple'):
    genos = load_genotriples_text(filename,genorepr,limit=limit,unique=unique,genome=genome)
  elif format in ('pb','prettybase'):
    genos = load_genotriples_prettybase(filename,limit=limit,unique=unique,genome=genome)
  elif format=='tbat':
    genos = load_genotriples_binary(filename,limit=limit,unique=unique,genome=genome)
  elif not format:
    raise ValueError, "Input file format for '%s' must be specified" % namefile(filename)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % format

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  return genos


def save_genostream(filename, genos, extra_args=None, **kwargs):
  '''
  Write genotype data to file

  @param  filename: a file name or file object
  @type   filename: str or file object
  @param     genos: genomatrix/genotriple stream
  @type      genos: sequence
  @param    format: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                    lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type     format: str
  @param  genorepr: string or representation object for text genotypes. Default is None
  @type   genorepr: str, UnphasedMarkerRepresentation or similar object
  @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type  mergefunc: callable
  @param  compress: flag indicating if a compressed format is desired. Default is True
  @type   compress: bool
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  format    = get_arg(args, ['format'])
  genorepr  = get_arg(args, ['genorepr']) or 'snp'
  mergefunc = get_arg(args, ['mergefunc'])
  compress  = get_arg(args, ['compress'], True)

  if format is None:
    format = guess_outformat(filename)

  if isinstance(genorepr,basestring):
    genorepr = get_genorepr(genorepr)

  if mergefunc is not None:
    genos = genos.merged(mergefunc)

  if format == 'ldat':
    save_genomatrix_text(filename, genos.as_ldat(mergefunc), genorepr)
  elif format == 'sdat':
    save_genomatrix_text(filename, genos.as_sdat(mergefunc), genorepr)
  elif format in ('tdat','trip','genotriple'):
    save_genotriples_text(filename, genos.as_genotriples(), genorepr)
  elif format in ('pb','prettybase'):
    genos = save_genotriples_prettybase(filename, genos.as_genotriples())
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
                    transform=None,genome=None,
                    mergefunc=None,limit=None):
  '''
  A driver for transforming multiple genodata files into different formats
  (ldat, sdat, trip, or genotriples), representations (...) and, depending
  on the presence and attributes of the transform object, performing
  operations on samples and loci such as exclude, include, and rename.  The
  results are then saved to the specified output file.

  @param     infiles: list of input file names or file objects
  @type      infiles: str or file objects
  @param    informat: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                      lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type     informat: str
  @param  ingenorepr: internal genotype representation for the input
  @type   ingenorepr: UnphasedMarkerRepresentation or similar object
  @param    outfiles: output file name or file object
  @type     outfiles: str or file object
  @param   outformat: format of output file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                      lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type    outformat: str
  @param outgenorepr: internal genotype representation for the output
  @type  outgenorepr: UnphasedMarkerRepresentation or similar object
  @param   transform: transformation object (optional)
  @type    transform: GenoTransform object
  @param   mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type    mergefunc: callable
  @param       limit: limit the number of samples loaded
  @type        limit: int or None

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\t\\tCT\\tTT\\n")
  >>> out  = StringIO()
  >>> transform_files([data],'ldat','snp',out,'tdat','marker')
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

  if isinstance(ingenorepr,basestring):
    ingenorepr = get_genorepr(ingenorepr)

  if isinstance(outgenorepr,basestring):
    outgenorepr = get_genorepr(outgenorepr)

  if not outgenorepr:
    outgenorepr = ingenorepr

  if genome is None:
    genome = Genome()
  elif isinstance(genome,basestring):
    genome = load_genome(genome)

  if isinstance(mergefunc,basestring):
    mergefunc = get_genomerger(mergefunc)

  genos = [ load_genostream(f,format=informat,genorepr=ingenorepr,limit=limit,genome=genome).transformed(transform) for f in infiles ]
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
  elif outformat in ('tdat','trip','genotriple','pb','prettybase','tbat'):
    genos = GenotripleStream.from_streams(genos,mergefunc=mergefunc)
  elif not outformat:
    raise ValueError, "Output file format for '%s' must be specified" % namefile(outfile)
  else:
    raise NotImplementedError,"File format '%s' is not supported" % outformat

  # Order again after merging, if necessary
  if n>1 and (transform.loci.order or transform.samples.order):
    genos = genos.transformed(order_loci=transform.loci.order,
                              order_samples=transform.samples.order)

  save_genostream(outfile,genos,format=outformat,genorepr=outgenorepr)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
