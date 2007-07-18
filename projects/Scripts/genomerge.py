# -*- coding: utf-8 -*-
'''
File:          genomerge.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)

Created:       June 13, 2006

Abstract:      This utility script merges the genotypes in .sdat matrix format
               according to the sample key mapping file. In addition to the
               merged genotype matrix file, the bookkeeping concordant statistics
               can also be output to files.

Compatibility: Python 2.4 and above

Requires:      biozilla

Version:       0.99

Revision:      $Id: genomerge.py 476 2007-01-19 19:38:39Z jacobske $

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

import sys
import csv
from   operator               import itemgetter
from   itertools              import izip,chain,islice,imap,repeat
from   biozilla.utils         import tally,hyphen
from   biozilla.genodata      import load_map,load_genostream,save_genomatrix
from   biozilla.genomerge     import get_genomerger, output_merge_statistics


def merge_dups(dupset, ind, locuskeys, mergefunc):
  return list(imap(mergefunc,repeat(ind),locuskeys,izip(*dupset)))


def merge_genotypes_safe(genos, simap, locuskeys, mergefunc):
  '''
  Safe merge function for input data that may not have unique row labels.
  Output order is deterministic, based on the order that input rows are
  observed.

  The major limitation is that the entire dataset must be loaded before
  any output can be generated.
  '''
  dups = {}
  inds = []

  for sample,geno in genos:
    if sample not in simap:
      continue

    ind = simap[sample]

    # Preserve order of individuals
    if ind not in dups:
      inds.append(ind)

    dups.setdefault(ind, []).append(geno)

  # Output rows, merging as necessary
  for ind in inds:
    genos = dups[ind]
    if len(genos) == 1:
      # Fast path for non-dups
      yield ind,genos[0]
    else:
      # Slow path for dups
      yield ind,merge_dups(genos, ind, locuskeys, mergefunc)


def merge_genotypes_unique(genos, simap, locuskeys, mergefunc):
  '''
  Fast and memory efficient merge function for input data with unique row
  labels.  If this assumption is violated, the resulting output may look
  very strange.  Output order is deterministic, based on the order that
  input rows are observed.

  The major optimization is that only the smallest portion of the dataset
  needed must be loaded before output is generated.  If a record has no
  duplicates, it can output without being buffered at all.  Any records with
  duplicates are kept in memory only until all duplicates are seen.
  '''
  ismap = tally(simap.itervalues())

  dups = {}
  inds = {}

  for i,(sample,geno) in enumerate(genos):
    if sample not in simap:
      continue

    ind   = simap[sample]
    edups = ismap[ind]

    # Fast path for non-dups
    if edups == 1:
      yield ind,geno
      continue

    # Slow path for dups
    inds.setdefault(ind, i)
    dups.setdefault(ind, []).append(geno)

    if len(dups[ind]) == edups:
      inds.pop(ind)
      yield ind,merge_dups(dups.pop(ind), ind, locuskeys, mergefunc)

  # Invert order dictionary
  inds = sorted(inds.iteritems(), key=itemgetter(1))

  # Output remaining rows that are missing some of their possible dups
  for i,ind in inds:
    genos = dups.get(ind,[])
    if len(genos) == 1:
      # Fast path for non-dups
      yield ind,genos[0]
    elif genos:
      # Slow path for dups
      yield ind,merge_dups(genos, ind, locuskeys, mergefunc)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-m', '--simap',       dest='simap',       metavar='FILE')
  parser.add_option('-f','--format',  dest='format', metavar='string', default='sdat',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat (default), trip or genotriple')
  parser.add_option('-o', '--outfile',     dest='outfile',     metavar='FILE',default='-')
  parser.add_option('-u', '--unique',      dest='unique',      action='store_true')
  parser.add_option('--merge', dest='merge', metavar='METHOD:T', default='vote:1',
                    help='Genotype merge algorithm and optional consensus threshold used to form a consensus genotypes. '
                         'Values=vote, first.  Value may be optionally followed by a colon and a threshold.  Default=vote:1')
  parser.add_option('-s','--samplemerge', dest='samplemerge', metavar='FILE',
                    help='Output the concordance statistics by sample to FILE (optional)')
  parser.add_option('-l','--locusmerge',  dest='locusmerge',  metavar='FILE',
                    help='Output the concordance statistics by locus to FILE (optional)')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    return

  merger = get_genomerger(options.merge)
  genos  = iter(load_genostream(hyphen(args[0],sys.stdin),options.format).as_sdat(merger))
  header = genos.next()
  simap  = load_map(options.simap)

  if options.unique:
    mergedgenos = merge_genotypes_unique(genos, simap, header, merger)
  else:
    mergedgenos = merge_genotypes_safe(genos, simap, header, merger)

  mergedgenos = chain([header],mergedgenos)

  outfile = hyphen(options.outfile,sys.stdout)
  save_genomatrix(outfile,mergedgenos,options.format)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__=='__main__':
  main()
