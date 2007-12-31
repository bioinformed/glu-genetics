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

from   collections               import defaultdict

from   glu.lib.utils             import pick
from   glu.lib.fileutils         import hyphen, load_map
from   glu.lib.genolib           import get_genorepr
from   glu.lib.genolib.io        import load_genostream, guess_outformat,               \
                                        TextGenomatrixWriter, TextGenotripleWriter,     \
                                        BinaryGenomatrixWriter, BinaryGenotripleWriter
from   glu.lib.genolib.genoarray import GenotypeArrayDescriptor, GenotypeArray


def genomatrix_multiplexer(genos, samplegroups, locusgroups, defaultsamplegroup, defaultlocusgroup):
  '''
  Sequentially split the contents of a genotype matrix into groups of rows
  and columns based on supplied mappings from row and column labels to group
  identifiers.  No buffering is performed, so partial results are returned
  tagged by the row and column group keys.
  '''
  if genos.format == 'ldat':
    rowgroups          = locusgroups
    columngroups       = samplegroups
    defaultrowgroup    = defaultlocusgroup
    defaultcolumngroup = defaultsamplegroup
  elif genos.format == 'sdat':
    rowgroups          = samplegroups
    columngroups       = locusgroups
    defaultrowgroup    = defaultsamplegroup
    defaultcolumngroup = defaultlocusgroup
  else:
    raise ValueError('Unknown genotype matrix format')

  header = genos.columns

  rdefault = [defaultrowgroup   ] if defaultrowgroup    else []
  cdefault = [defaultcolumngroup] if defaultcolumngroup else []

  if columngroups is not None:
    groupcols = defaultdict(list)
    for i,colkey in enumerate(genos.columns):
      for columngroup in columngroups.get(colkey, cdefault):
        if columngroup:
          groupcols[columngroup].append(i)

    if genos.format == 'sdat':
      groupcols = [ (key,indices,
                         GenotypeArrayDescriptor(pick(genos.models,indices)),
                         pick(genos.columns,indices))
                    for key,indices in groupcols.iteritems() ]
    else:
      descrcache = {}
      groupcols = [ (key,indices,
                         None,
                         pick(genos.columns,indices))
                    for key,indices in groupcols.iteritems() ]

  # This used to be a single loop, but performance is an issue and breaking
  # out the cases helps significantly.  One day Python will have a JIT and
  # take care of this for me...
  if columngroups and rowgroups:
    for rowkey,row in genos:
      models = row.descriptor.models
      row = row[:]
      for rowgroup in rowgroups.get(rowkey) or rdefault:
        if not rowgroup:
          continue
        for columngroup,indices,descr,header in groupcols:
          if not descr:
            n     = len(indices)
            model = models[0]
            descr = descrcache.get( (model,n) )
          if not descr:
            descr = descrcache[model,n] = GenotypeArrayDescriptor([model]*n)
          grow = GenotypeArray(descr,pick(row,indices))
          yield (rowgroup,columngroup),header,(rowkey,grow)

  elif columngroups:
    for rowkey,row in genos:
      models = row.descriptor.models
      row = row[:]
      for columngroup,indices,descr,header in groupcols:
        if not descr:
          n     = len(indices)
          model = models[0]
          descr = descrcache.get( (model,n) )
        if not descr:
          descr = descrcache[model,n] = GenotypeArrayDescriptor([model]*n)
        grow = GenotypeArray(descr,pick(row,indices))
        yield (None,columngroup),header,(rowkey,grow)

  elif rowgroups:
    for rowkey,row in genos:
      for rowgroup in rowgroups.get(rowkey) or rdefault:
        if rowgroup:
          yield (rowgroup,None),header,(rowkey,row)

  elif defaultrowgroup or defaultcolumngroup:
    key = (defaultrowgroup,defaultcolumngroup)
    for row in genos:
      yield key,header,row


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
          yield (samplegroup,locusgroup),None,(sample,locus,geno)

  elif samplegroups:
    for sample,locus,geno in triples:
      for samplegroup in samplegroups.get(sample, sdefault):
        yield (samplegroup,None),None,(sample,locus,geno)

  elif locusgroups:
    for sample,locus,geno in triples:
      for locusgroup in locusgroups.get(locus, ldefault):
        yield (None,locusgroup),None,(sample,locus,geno)

  elif sampledefault or locusdefault:
    key = (sampledefault,locusdefault)
    for sample,locus,geno in triples:
      yield key,None,(sample,locus,geno)


def getWriter(filename,format,genome=None,phenome=None,header=None,genorepr=None,maxrows=None):
  if maxrows:
    return RollingWriter(filename,format,genome,phenome,maxrows,header,genorepr)
  elif format in ('ldat','sdat'):
    return TextGenomatrixWriter(filename,format,header,genorepr=genorepr)
  elif format == 'lbat':
    return BinaryGenomatrixWriter(filename,'ldat',header,genome,phenome)
  elif format == 'sbat':
    return BinaryGenomatrixWriter(filename,'sdat',header,genome,phenome)
  elif format in ('tdat','trip','genotriple'):
    return TextGenotripleWriter(filename,genorepr=genorepr)
  elif format == 'tbat':
    return BinaryGenotripleWriter(filename,genome,phenome)
  else:
    raise ValueError('Unknown format')


class RollingWriter(object):
  '''
  A wrapper around Text and Binary Writer objects that accepts a maximum
  number of rows per file.  Once that limit is reached, another file is
  opened.
  '''
  def __init__(self, filename, format, genome, phenome, maxrows, header=None, genorepr=None):
    self.filename = filename
    self.format   = format
    self.genome   = genome
    self.phenome  = phenome
    self.maxrows  = maxrows
    self.header   = header
    self.genorepr = genorepr

    self.rows     = sys.maxint
    self.cycles   = 0
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

    self.writer = getWriter(filename,self.format,self.genome,self.phenome,
                                     header=self.header,genorepr=self.genorepr)

  def writerow(self, *row):
    if self.rows >= self.maxrows:
      self.cycle()
    self.rows += 1
    self.writer.writerow(*row)

  def writerows(self,rows):
    for row in rows:
      self.writerow(*row)

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


class FileMap(object):
  '''
  Container for TextGenomatrixWriter and RollingTextGenomatrixWriter objects
  stored by a key tuple.
  '''
  def __init__(self, prefix, suffix, format, genome, phenome, outformat=None, genorepr=None, maxrows=None):
    self.writers   = {}
    self.prefix    = prefix
    self.suffix    = suffix
    self.format    = format
    self.genome    = genome
    self.phenome   = phenome
    self.outformat = outformat or format
    self.genorepr  = genorepr
    self.maxrows   = maxrows

  def emit(self, keys, header, row):
    self.get_writer(keys, header).writerow(*row)

  def emit_sequence(self, seq):
    for keys,header,row in seq:
      self.get_writer(keys, header).writerow(*row)

  def get_writer(self, keys, header=None):
    writer = self.writers.get(keys)

    if writer is None:
      filename = build_filename(self.prefix, self.suffix, keys)

      writer = getWriter(filename,self.outformat,self.genome,self.phenome,
                                  header=header,genorepr=self.genorepr,
                                  maxrows=self.maxrows)

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


def split(genos, outformat, prefix, suffix, options):
  header       = genos.columns if genos.format in ('sdat','ldat') else None
  locusgroups  = load_map(options.locusgroups, unique=False,default=options.defaultlocusgroup)  if options.locusgroups  else None
  samplegroups = load_map(options.samplegroups,unique=False,default=options.defaultsamplegroup) if options.samplegroups else None

  if samplegroups is not None or locusgroups is not None:
    filecache = FileMap(prefix,suffix,genos.format,genos.genome,genos.phenome,outformat=outformat,
                        genorepr=options.genorepr,maxrows=options.maxrows)

    if genos.format in ('sdat','ldat'):
      mplx = genomatrix_multiplexer(genos,samplegroups,locusgroups,
                                    options.defaultsamplegroup,options.defaultlocusgroup)
    else:
      mplx = genotriple_multiplexer(genos,samplegroups,locusgroups,
                                    options.defaultsamplegroup,options.defaultlocusgroup)

    with filecache:
      filecache.emit_sequence(mplx)

  elif options.maxrows:
    writer = RollingWriter('%s.%s' % (prefix,suffix),outformat,genome=genos.genome,phenome=genos.phenome,
                                     header=header,genorepr=options.genorepr,maxrows=options.maxrows)
    with writer:
      writer.writerows(genos)

  else:
    sys.stderr.write('Terminating: No grouping or splitting specified\n')


def option_parser():
  import optparse

  usage = 'usage: %prog [options] matrixfile'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f','--format', dest='format',
                    help='Input format for genotype data. Values=hapmap, ldat, sdat, trip, or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
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

  # FIXME: Must know about augmented filenames
  if filename == '-':
    sys.stderr.write('Error: A filename template must be specified when taking input from stdin')
    return

  prefix,suffix = split_fullname(filename,options.destdir)

  options.genorepr = get_genorepr(options.genorepr)

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                  genome=options.loci,hyphen=sys.stdin)

  outformat = guess_outformat(args[0]) or options.format or genos.format
  split(genos, outformat, prefix, suffix, options)


if __name__ == '__main__':
  main()
