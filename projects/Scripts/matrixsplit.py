# -*- coding: utf-8 -*-
'''
File:

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       June 7, 2006

Abstract:      This utility script can split an input matrix into multiple
               output files by row and column groupings

Compatibility: Python 2.4 and above

Requires:      biozilla

Version:       0.99

Revision:      $Id: matrixsplit.py 410 2006-10-24 19:55:16Z jacobske $

Copyright (c) 2006 BioInformed Consulting Services.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
'''

__version__ = '0.99'

import csv
import sys
import os

from   operator           import itemgetter
from   itertools          import chain, islice

from   biozilla.utils     import pick
from   biozilla.fileutils import autofile, hyphen, load_map
from   biozilla.genodata  import load_genomatrix, save_genomatrix
from   biozilla.genoarray import snp_acgt


def multiplexer(matrix, format, rmap, cmap, rdefault, cdefault):

  header = matrix.next()

  groupcols = {}
  for i,colkey in enumerate(header):
    colgroup = cmap.get(colkey, cdefault)
    group = groupcols.setdefault(colgroup,([],[format]))
    group[0].append(i)
    group[1].append(colkey)

  groupcols = [ (key,inds,cols) for key,(inds,cols) in groupcols.iteritems() ]

  for rowkey,row in matrix:
    rowgroup = rmap.get(rowkey, rdefault)
    for colgroup, colinds, header in groupcols:
      if not colgroup and not rowgroup:
        continue
      rowdata = pick(row,colinds)
      yield (rowgroup,colgroup),header,rowkey,rowdata


class RollingWriter(object):
  def __init__(self, filename, header, maxrows):
    self.rows     = sys.maxint
    self.cycles   = 0
    self.maxrows  = maxrows
    self.filename = filename
    self.header   = header

  def cycle(self):
    self.cycles += 1
    self.rows    = 0

    try:
      filename = self.filename % self.cycles
    except TypeError:
      prefix,suffix = split_fullname(self.filename,'')
      filename = build_filename(prefix, suffix, (self.cycles,) )

    self.writer = csv.writer(autofile(filename, 'w'), dialect='excel-tab')
    self.writer.writerow(self.header)

  def writerow(self,row):
    if self.rows >= self.maxrows:
      self.cycle()
    self.rows += 1
    self.writer.writerow(row)

  def writerows(self,rows):
    for row in rows:
      self.writerow(row)


class MatrixFileCache(object):
  def __init__(self, prefix, suffix, maxrows=None):
    self.outfiles = {}
    self.prefix  = prefix
    self.suffix  = suffix
    self.maxrows = maxrows

  def emit(self, keys, header, rowkey, data):
    outfile = self.get_file(keys, header)
    data = chain([rowkey], snp_acgt.genos_str(data))
    outfile.writerow(list(data))

  def get_file(self, keys, header):
    outfile = self.outfiles.get(keys)

    if outfile is None:
      outfilename = build_filename(self.prefix, self.suffix, keys)

      if self.maxrows:
        outfile = RollingWriter(outfilename, header, self.maxrows)
      else:
        outfile = csv.writer(autofile(outfilename, 'w'), dialect='excel-tab')
        outfile.writerow(header)

      self.outfiles[keys] = outfile

    return outfile


def build_filename(prefix, suffix, keys):
  filename = prefix

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
  dirname  = os.path.dirname(filename)
  basename = os.path.basename(filename)

  parts = basename.split('.',1)

  if destdir:
    dirname = destdir

  prefix = os.path.join(dirname,parts[0])

  suffix = ''
  if len(parts) == 2:
    suffix = parts[1]

  return prefix,suffix


def option_parser():
  import optparse

  usage = 'usage: %prog [options] matrixfile'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-d', '--destdir', dest='destdir', default='',
                    help='Destination directory for output files.  Default="" to write to input file directory.')
  parser.add_option('--maxrows', dest='maxrows', metavar='N', type='int',
                    help='Split matrix output so that each contains at most N rows of data')
  parser.add_option('--rowgroups', dest='rowgroups', metavar='FILE',
                    help='File containing the map between the row key in the matrix file and the row key group')
  parser.add_option('--colgroups', dest='colgroups', metavar='FILE',
                    help='File containing the map between the col key in the matrix file and the col key group')
  parser.add_option('--defaultrowgroup', dest='defaultrowgroup', metavar='NAME',
                    help='Default group for any row mapped in rowgroups')
  parser.add_option('--defaultcolumngroup', dest='defaultcolumngroup', metavar='NAME',
                    help='Default group for any column mapped in colgroups')
  parser.add_option('--template', dest='template', metavar='NAME', type='str',
                    help='The template for names of the output files')

  return parser


def main():
  parser = option_parser()
  options, args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    return

  if options.template:
    prefix,suffix = options.template.split('.',1)
  else:
    prefix,suffix = split_fullname(args[0],options.destdir)

  format,matrix = load_genomatrix(hyphen(args[0],sys.stdin))

  rmap = {}
  cmap = {}

  if options.rowgroups:
    rmap = load_map(options.rowgroups)
  if options.colgroups:
    cmap = load_map(options.colgroups)

  if rmap or cmap:
    filecache = MatrixFileCache(prefix, suffix, options.maxrows)
    mplx = multiplexer(matrix, format, rmap, cmap, options.defaultrowgroup, options.defaultcolumngroup)
    for m in mplx:
      filecache.emit(*m)
  elif options.maxrows:
    header = list(chain([format],matrix.next()))
    writer = RollingWriter('%s_part.%s' % (prefix,suffix), header, options.maxrows)
    writer.writerows( [rowkey]+list(snp_acgt.genos_str(data)) for rowkey,data in matrix)
  else:
    print >> sys.stderr,'Terminating: No split specified'


if __name__ == '__main__':
  main()
