# -*- coding: utf-8 -*-
'''
File:          split.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       June 7, 2006

Abstract:      This utility script can split an input matrix into multiple
               output files by row and column groupings

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id$
'''

from __future__ import with_statement

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys

from   glu.lib.utils      import pick
from   glu.lib.fileutils  import hyphen, load_map
from   glu.lib.genolib    import get_genorepr
from   glu.lib.genolib.io import load_genostream, TextGenomatrixWriter, TextGenotripleWriter


def genomatrix_multiplexer(format, header, matrix, samplegroups, locusgroups,
                               defaultsamplegroup, defaultlocusgroup):
  '''
  Sequentially split the contents of a genotype matrix into groups of rows
  and columns based on supplied mappings from row and column labels to group
  identifiers.  No buffering is performed, so partial results are returned
  tagged by the row and column group keys.
  '''
  if format == 'ldat':
    rowgroups          = locusgroups
    columngroups       = samplegroups
    defaultrowgroup    = defaultlocusgroup
    defaultcolumngroup = defaultsamplegroup
  elif format == 'sdat':
    rowgroups          = samplegroups
    columngroups       = locusgroups
    defaultrowgroup    = defaultsamplegroup
    defaultcolumngroup = defaultlocusgroup
  else:
    raise ValueError('Unknown genotype matrix format')

  rdefault = [defaultrowgroup   ] if defaultrowgroup    else []
  cdefault = [defaultcolumngroup] if defaultcolumngroup else []

  if columngroups is not None:
    groupcols = {}
    for i,colkey in enumerate(header):
      for columngroup in columngroups.get(colkey, cdefault):
        if columngroup:
          groupcols.setdefault(columngroup,[]).append(i)

    groupcols = [ (key,indices,pick(header,indices)) for key,indices in groupcols.iteritems() ]

  # This used to be a single loop, but performance is an issue and breaking
  # out the cases helps significantly.  One day Python will have a JIT and
  # take care of this for me...
  if columngroups and rowgroups:
    for rowkey,genos in matrix:
      for rowgroup in rowgroups.get(rowkey, rdefault):
        for columngroup,indices,header in groupcols:
          yield (rowgroup,columngroup),header,rowkey,pick(genos[:],indices)

  elif columngroups:
    for rowkey,genos in matrix:
      for columngroup,indices,header in groupcols:
        yield (None,columngroup),header,rowkey,pick(genos[:],indices)

  elif rowgroups:
    for rowkey,genos in matrix:
      for rowgroup in rowgroups.get(rowkey, rdefault):
        yield (rowgroup,None),header,rowkey,genos

  elif defaultrowgroup or defaultcolumngroup:
    key = (defaultrowgroup,defaultcolumngroup)
    for rowkey,genos in matrix:
      yield key,header,rowkey,genos


def genotriple_multiplexer(triples, samplegroups, locusgroups, sampledefault, locusdefault):
  '''
  Sequentially split the contents of a genotype triple stream into groups
  based on sample and locus using supplied mappings from row and column
  labels to group identifiers.  No buffering is performed, so partial
  results are returned tagged by group keys.
  '''
  sdefault = [sampledefault] if sampledefault else []
  ldefault = [locusdefault ] if locusdefault  else []

  # This used to be be a single loop, but performance is an issue and
  # breaking out the cases helps significantly.  One day Python will have a
  # JIT and take care of this for me...
  if samplegroups and locusgroups:
    for sample,locus,geno in triples:
      for samplegroup in samplegroups.get(sample, sdefault):
        for locusgroup in samplegroups.get(locus, ldefault):
          yield (samplegroup,locusgroup),sample,locus,geno

  elif samplegroups:
    for sample,locus,geno in triples:
      for samplegroup in samplegroups.get(sample, sdefault):
        yield (samplegroup,None),sample,locus,geno

  elif locusgroups:
    for sample,locus,geno in triples:
      for locusgroup in locusgroups.get(locus, ldefault):
        yield (None,locusgroup),sample,locus,geno

  elif sampledefault or locusdefault:
    key = (sampledefault,locusdefault)
    for sample,locus,geno in triples:
      yield key,sample,locus,geno


class RollingTextGenomatrixWriter(object):
  '''
  A wrapper around TextGenomatrixWriter that accepts a maximum number of
  rows per file.  Once that limit is reached, another filed is opened.
  '''
  def __init__(self, filename, format, header, genorepr, maxrows):
    self.filename = filename
    self.format   = format
    self.header   = header
    self.genorepr = genorepr

    self.rows     = sys.maxint
    self.cycles   = 0
    self.maxrows  = maxrows
    self.writer   = None

  def cycle(self):
    self.cycles += 1
    self.rows    = 0

    self.close()

    try:
      filename = self.filename % self.cycles
    except TypeError:
      prefix,suffix = split_fullname(self.filename,'')
      filename = build_filename(prefix + '_part', suffix, (self.cycles,) )

    self.writer = TextGenomatrixWriter(filename,self.format,self.header,self.genorepr)

  def writerow(self,rowkey,genos):
    if self.rows >= self.maxrows:
      self.cycle()
    self.rows += 1
    self.writer.writerow(rowkey,genos)

  def writerows(self,rows):
    for rowkey,genos in rows:
      self.writerow(rowkey, genos)

  def close(self):
    if self.writer is not None:
      self.writer.close()
      self.writer = None

  def __enter__(self):
    return self

  def __exit__(self, *exc_info):
    self.close()

  def __del__(self):
    self.close()


class RollingTextGenotripleWriter(object):
  '''
  A wrapper around TextGenotripleWriter that accepts a maximum number of
  rows per file.  Once that limit is reached, another filed is opened.
  '''
  def __init__(self, filename, genorepr, maxrows):
    self.filename = filename
    self.genorepr = genorepr

    self.rows     = sys.maxint
    self.cycles   = 0
    self.maxrows  = maxrows
    self.writer   = None

  def cycle(self):
    self.cycles += 1
    self.rows    = 0

    self.close()

    try:
      filename = self.filename % self.cycles
    except TypeError:
      prefix,suffix = split_fullname(self.filename,'')
      filename = build_filename(prefix + '_part', suffix, (self.cycles,) )

    self.writer = TextGenotripleWriter(filename,self.genorepr)

  def writerow(self, sample, locus, geno):
    if self.rows >= self.maxrows:
      self.cycle()
    self.rows += 1
    self.writer.writerow(sample,locus,geno)

  def writerows(self,triples):
    for sample,locus,geno in triples:
      self.writerow(sample,locus,geno)

  def close(self):
    if self.writer is not None:
      self.writer.close()
      self.writer = None

  def __enter__(self):
    return self

  def __exit__(self, *exc_info):
    self.close()

  def __del__(self):
    self.close()


class TextGenomatrixFileMap(object):
  '''
  Container for TextGenomatrixWriter and RollingTextGenomatrixWriter objects
  stored by a key tuple.
  '''
  def __init__(self, prefix, suffix, format, header, genorepr, maxrows=None):
    self.writers  = {}
    self.prefix   = prefix
    self.suffix   = suffix
    self.format   = format
    self.header   = header
    self.genorepr = genorepr
    self.maxrows  = maxrows

  def emit(self, keys, header, rowkey, genos):
    self.get_writer(keys, header).writerow(rowkey,genos)

  def emit_sequence(self, seq):
    for keys,header,rowkey,genos in seq:
      self.get_writer(keys, header).writerow(rowkey,genos)

  def get_writer(self, keys, header):
    writer = self.writers.get(keys)

    if writer is None:
      filename = build_filename(self.prefix, self.suffix, keys)

      if self.maxrows:
        writer = RollingTextGenomatrixWriter(filename,self.format,header,self.genorepr,self.maxrows)
      else:
        writer = TextGenomatrixWriter(filename,self.format,header,self.genorepr)

      self.writers[keys] = writer

    return writer

  def close(self):
    for writer in self.writers.itervalues():
      writer.close()
    self.writers = {}

  def __enter__(self):
    return self

  def __exit__(self, *exc_info):
    self.close()

  def __del__(self):
    self.close()


class TextGenotripleFileMap(object):
  '''
  Container for TextGenotripleWriter and RollingTextGenotripleWriter objects stored by a
  key tuple.
  '''
  def __init__(self, prefix, suffix, genorepr, maxrows=None):
    self.writers  = {}
    self.prefix   = prefix
    self.suffix   = suffix
    self.genorepr = genorepr
    self.maxrows  = maxrows

  def emit(self, keys, sample, locus, geno):
    self.get_writer(keys).writerow(sample,locus,geno)

  def emit_sequence(self, seq):
    for keys,sample,locus,geno in seq:
      self.get_writer(keys).writerow(sample,locus,geno)

  def get_writer(self, keys):
    writer = self.writers.get(keys)

    if writer is None:
      filename = build_filename(self.prefix, self.suffix, keys)

      if self.maxrows:
        writer = RollingTextGenotripleWriter(filename,self.genorepr,self.maxrows)
      else:
        writer = TextGenotripleWriter(filename,self.genorepr)

      self.writers[keys] = writer

    return writer

  def close(self):
    for writer in self.writers.itervalues():
      writer.close()
    self.writers = {}

  def __enter__(self):
    return self

  def __exit__(self, *exc_info):
    self.close()

  def __del__(self):
    self.close()


def build_filename(prefix, suffix, keys):
  filename = prefix

  keys = [k for k in keys if k]

  if not keys:
    raise ValueError('Internal error: cannot construct filename from null keys')

  for key in keys:
    if not key:
      continue
    if filename:
      filename += '_'
    filename += str(key)

  if suffix:
    filename += '.%s' % suffix

  return filename


def split_fullname(filename,destdir):
  # Set destination directory
  if destdir:
    dirname = destdir
  else:
    dirname  = os.path.dirname(filename)

  # Get filename
  filename = os.path.basename(filename)

  # Split filename into 1 or 2 parts up to the first '.'
  parts = filename.split('.',1)

  # Combine dirname and up to the first '.' of filename as prefix
  prefix = os.path.join(dirname,parts[0])

  # Suffix the remainder of filename after the first '.'
  suffix = '' if len(parts) == 1 else parts[1]

  return prefix,suffix


def matrix_split(matrix, genorepr, prefix, suffix, options):
  format       = matrix.format
  header       = matrix.columns
  maxrows      = options.maxrows
  locusgroups  = load_map(options.locusgroups, unique=False) if options.locusgroups  else None
  samplegroups = load_map(options.samplegroups,unique=False) if options.samplegroups else None

  if samplegroups is not None or locusgroups is not None:
    with TextGenomatrixFileMap(prefix,suffix,format,header,genorepr,maxrows) as filecache:
      mplx = genomatrix_multiplexer(format,header,matrix,samplegroups,locusgroups,
                                    options.defaultsamplegroup,
                                    options.defaultlocusgroup)

      filecache.emit_sequence(mplx)

  elif maxrows:
    writer = RollingTextGenomatrixWriter('%s.%s' % (prefix,suffix),format,header,genorepr,maxrows)
    writer.writerows(matrix)

  else:
    sys.stderr.write('Terminating: No grouping or splitting specified\n')


def triple_split(triples, genorepr, prefix, suffix, options):
  maxrows      = options.maxrows
  locusgroups  = load_map(options.locusgroups, unique=False) if options.locusgroups  else None
  samplegroups = load_map(options.samplegroups,unique=False) if options.samplegroups else None

  if locusgroups is not None or samplegroups is not None:
    with TextGenotripleFileMap(prefix,suffix,genorepr,maxrows) as filecache:
      mplx = genotriple_multiplexer(triples,samplegroups,locusgroups,
                                    options.defaultsamplegroup,
                                    options.defaultlocusgroup)

      filecache.emit_sequence(mplx)

  elif maxrows:
    writer = RollingTextGenotripleWriter('%s.%s' % (prefix,suffix),genorepr,maxrows)
    writer.writerows(triples)

  else:
    sys.stderr.write('Terminating: No grouping or splitting specified\n')


def option_parser():
  import optparse

  usage = 'usage: %prog [options] matrixfile'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f','--format', dest='format',
                    help='Input format for genotype data. Values=hapmap, ldat, sdat, trip, or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='genotype representation.  Values=snp (default), hapmap, marker')

  parser.add_option('-d', '--destdir', dest='destdir', default='',
                    help='Destination directory for output files.  Write to input file directory by default.')
  parser.add_option('--maxrows', dest='maxrows', metavar='N', type='int',
                    help='Split matrix output so that each contains at most N rows of data')
  parser.add_option('--locusgroups', dest='locusgroups', metavar='FILE',
                    help='map from locus name to locus group')
  parser.add_option('--samplegroups', dest='samplegroups', metavar='FILE',
                    help='map from samples name to sample group')
  parser.add_option('--defaultsamplegroup', dest='defaultsamplegroup', metavar='NAME',
                    help='Default group for any unmapped sample')
  parser.add_option('--defaultlocusgroup', dest='defaultlocusgroup', metavar='NAME',
                    help='Default group for any unmapped sample')
  parser.add_option('--template', dest='template', metavar='NAME', type='str',
                    help='The template for names of the output files')

  return parser


def main():
  parser = option_parser()
  options, args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    return

  filename = options.template if options.template else args[0]
  if filename == '-':
    sys.stderr.write('Error: A filename template must be specified when taking input from stdin')
    return
  prefix,suffix = split_fullname(filename,options.destdir)

  genorepr = get_genorepr(options.genorepr)
  infile   = hyphen(args[0],sys.stdin)
  genos    = load_genostream(infile,format=options.format,genorepr=genorepr)

  if genos.format in ('sdat','ldat'):
    matrix_split(genos, genorepr, prefix, suffix, options)
  elif genos.format == 'genotriple':
    triple_split(genos, genorepr, prefix, suffix, options)
  else:
    raise ValueError('Unsupported input file format %s' % genos.format)


if __name__ == '__main__':
  main()
