# -*- coding: utf-8 -*-
'''
File:

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       February 27, 2006

Abstract:

Compatibility: Python 2.4 and above

Requires:      No external dependencies, yet...

Version:       0.99

Revision:      $Id: ill2raw.py 359 2006-09-19 16:51:59Z jacobske $

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

import re
import sys
import csv
from   itertools         import islice, chain, groupby
from   operator          import attrgetter, itemgetter
from   biozilla.utils    import autofile, hyphen
from   biozilla.xtab     import rowsby
from   biozilla.genodata import load_list,load_map,save_genomatrix


DATA_HEADER = ['SNP Name','Sample ID','Allele1 - Forward','Allele2 - Forward','GC Score','X','Y','X Raw','Y Raw']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] sorted_genofile'
  parser = optparse.OptionParser(usage=usage,version='%%prog %s' % __version__)

  #parser.add_option('--gcthreshold', dest='gcthreshold', type='float', metavar='N', default=0)
  parser.add_option('-o','--output', dest='output',                    metavar='FILE')
  parser.add_option('--samplemap',   dest='samplemap',                 metavar='FILE')
  parser.add_option('--locusmap',    dest='locusmap',                  metavar='FILE')
  parser.add_option('--gcthreshold', dest='gcthreshold', type='float', metavar='N')
  parser.add_option('--completion',  dest='completion',  type='float', metavar='N')
  parser.add_option('--limit',       dest='limit',       type='int',   metavar='N')
  parser.add_option('--columns',     dest='columns',                   metavar='FILE')
  parser.add_option('--rowkey',      dest='rowkey',      type='int',   metavar='N')
  parser.add_option('--colkey',      dest='colkey',      type='int',   metavar='N')
  parser.add_option('--datakey',     dest='datakey',     type='int',   metavar='N')

  return parser


def load_illumina_genotypes(filename):
  genofile = csv.reader(autofile(hyphen(filename,sys.stdin)),dialect='excel-tab')

  if 0:
    headings = []
    header = genofile.next()
    if header == ['[Header]']:
      for line in genofile:
        if line == ['[Data]']:
          break
        headings.append(line)
      header = genofile.next()

    print header
    assert header == DATA_HEADER

  for row in genofile:
    if row == DATA_HEADER:
      continue

    lname  = row[0]
    sample = row[1]
    geno   = intern( (row[2]+row[3]).replace('-','') )
    gc     = float(row[4])
    yield lname,sample,geno,gc


def remap_loci(genos,locusmap):
  for lname,sample,geno,gc in genos:
    yield locusmap[lname],sample,geno,gc


def remap_samples(genos,samplemap):
  for lname,sample,geno,gc in genos:
    yield lname,samplemap.get(sample,sample),geno,gc


def datastream(data,rowkey,colkey,datakey):
  for row in data:
    yield row[rowkey],row[colkey],row[datakey]


def merge_genos(ind,locus,genos):
  if not genos:
    return ''

  if len(genos) == 1:
    return genos[0]

  genos = set(genos)
  genos.discard('')

  if not genos:
    return ''
  elif len(genos) > 1:
    print 'BAD GENO',ind,locus,genos
    return ''
  else:
    return genos.pop()


# FIXME: gcthreshold is not honored!!!
def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1 or None in [options.columns,options.colkey,options.rowkey,options.output]:
    parser.print_help(sys.stderr)
    return

  columns = load_list(options.columns)
  genos   = load_illumina_genotypes(args[0])

  if options.locusmap:
    genos = remap_loci(genos, load_map(options.locusmap))

  if options.samplemap:
    smap = load_map(options.samplemap)
    columns = [ smap.get(c,c) for c in columns ]
    genos = remap_samples(genos, smap)

  data = datastream(genos,options.rowkey,options.colkey,options.datakey)

  rows = rowsby(data, columns, itemgetter(0), itemgetter(1), itemgetter(2), merge_genos)

  out = hyphen(options.output,sys.stdout)
  save_genomatrix(out,rows,genorepr=list)


if __name__ == '__main__':
  main()
