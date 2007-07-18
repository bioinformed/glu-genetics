# -*- coding: utf-8 -*-
'''
File:          transform.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-06-29

Abstract:      Performs various transformation on a genotype file

Requires:      Python 2.5, biozilla

Revision:      $Id: transform.py 336 2006-08-09 18:57:08Z staatsb $
'''

__version__ = '0.2'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv
from   itertools              import islice, imap,chain
from   operator               import itemgetter

from   biozilla.xtab          import xtab_list
from   biozilla.utils         import autofile
from   biozilla.genoarray     import snp_acgt
# FIXME: Remove global import once genodata has better uber-functions
from   biozilla.genodata      import *


# FIXME: Merge with genomerge
class Merger(object):
  def __init__(self):
    self.discordances = {}

  def __call__(self,row,col,genos):
    if not genos:
      return 0
    elif len(genos)==1:
      return genos[0]

    genos = set(genos)
    genos.discard(0)

    if not genos:
      return 0
    elif len(genos) == 1:
      return genos.pop()

    self.discordances[row,col] = genos

    return 0


# FIXME: May belong in genodata
def load_genotriples(filename,limit):
  rows = csv.reader(autofile(filename),dialect='excel-tab')

  # FIXME: incorrect implementation of limit for a triple stream
  if limit:
    rows = islice(rows,limit+1)

  for sample,locus,geno in rows:
    yield sample,locus,snp_acgt.byte(tuple(geno))


# FIXME: May belong in genodata
def save_genotriples(filename,genos):
  out = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
  for sample,locus,geno in genos:
    out.writerow( [ sample,locus,snp_acgt.geno_str(geno) ] )


def write_from_ldat(filename, genos, outfmt):
  if outfmt=='ldat':
    save_genomatrix(filename,genos)

  elif outfmt=='sdat':
    columns = genos.next()
    columns,genos = transpose_generator(list(genos), columns)
    save_genomatrix(filename,chain([columns],genos))

  elif outfmt=='trip':
    genos = build_genotriples_by_locus(genos)
    save_genotriples(filename,genos)

  else:
    raise NotImplementedError,'ldat to %s format conversion is not supported' % outfmt


def write_from_sdat(filename, genos, outfmt):
  if outfmt=='sdat':
    save_genomatrix(filename,genos)

  elif outfmt=='ldat':
    columns = genos.next()
    columns,genos = transpose_generator(list(genos), columns)
    save_genomatrix(filename,chain([columns],genos))

  elif outfmt=='trip':
    genos = build_genotriples_by_sample(genos)
    save_genotriples(filename,genos)

  else:
    raise NotImplementedError,'sdat to %s format conversion is not supported' % outfmt


def write_from_trip(filename, genos, outfmt):
  if outfmt=='trip':
    save_genotriples(filename,genos)

  elif outfmt=='sdat':
    merge = Merger()
    genos = xtab_list(genos, itemgetter(0), itemgetter(1), itemgetter(2), merge)
    columns = genos.next()
    save_genomatrix(filename,chain([columns],genos))

  elif outfmt=='ldat':
    merge = Merger()
    genos = xtab_list(genos, itemgetter(1), itemgetter(0), itemgetter(2), merge)
    columns = genos.next()
    save_genomatrix(filename,chain([columns],genos))

  else:
    raise NotImplementedError,'trip to %s format conversion is not supported' % outfmt


# FIXME: Not implemented + difficult to do for missing columns
def filter_genomatrix_missing(genos):
  raise NotImplementedError


def transform_genomatrix_by_locus(genos, options):
  if options.includelocus:
    genos = filter_genomatrix_by_row(genos,load_map(options.includelocus))
  if options.excludelocus:
    genos = filter_genomatrix_by_row(genos,load_map(options.excludelocus),exclude=True)
  if options.renamelocus:
    genos = rename_genomatrix_row(genos, load_map(options.renamelocus))
  if options.excludesamples:
    genos = filter_genomatrix_by_column(genos,load_map(options.excludesamples),exclude=True)
  if options.includesamples:
    genos = filter_genomatrix_by_column(genos,load_map(options.includesamples))
  if options.renamesamples:
    genos = rename_genomatrix_column(genos, load_map(options.renamesamples))
  if options.filtermissing:
    genos = filter_genomatrix_missing(genos)

  return genos


def transform_genomatrix_by_sample(genos, options):
  if options.excludesamples:
    genos = filter_genomatrix_by_row(genos,load_map(options.excludesamples),exclude=True)
  if options.includesamples:
    genos = filter_genomatrix_by_row(genos,load_map(options.includesamples))
  if options.renamesamples:
    genos = rename_genomatrix_row(genos, load_map(options.renamesamples))
  if options.excludelocus:
    genos = filter_genomatrix_by_column(genos,load_map(options.excludelocus),exclude=True)
  if options.includelocus:
    genos = filter_genomatrix_by_column(genos,load_map(options.includelocus))
  if options.renamelocus:
    genos = rename_genomatrix_column(genos, load_map(options.renamelocus))
  if options.filtermissing:
    genos = filter_genomatrix_missing(genos)

  return genos


def transform_genotriples(genos, options):
  if options.filtermissing:
    genos = filter_genotriples_missing(genos)

  if options.excludesamples or options.excludelocus:
    samples = load_map(options.excludesamples) if options.excludesamples else None
    loci    = load_map(options.excludelocus)   if options.excludelocus   else None
    genos   = filter_genotriples(genos,samples,loci,exclude=True)

  if options.includesamples or options.includelocus:
    samples = load_map(options.includesamples) if options.includesamples else None
    loci    = load_map(options.includelocus)   if options.includelocus   else None
    genos   = filter_genotriples(genos,samples,loci)

  if options.samplemap or options.locusmap:
    samples = load_map(options.renamesamples) if options.locusmap    else None
    loci    = load_map(options.renamelocus)   if options.renamelocus else None
    genos   = rename_genotriples(genos,samples,loci)

  return genos


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-i','--inputformat',  dest='informat', metavar='string',
                    help='The file input format for genotype data. Values=ldat, sdat, trip')
  parser.add_option('-o','--outputformat',  dest='outformat', metavar='string',
                    help='The file output format for genotype data. Values=ldat (default), sdat, trip')

  parser.add_option('-c', '--filtermissing', action='store_true', dest='filtermissing',
                    help='Filters out the samples or loci with missing data')

  parser.add_option('-n', '--includesamples', dest='includesamples', metavar='FILE',
                    help='Include list for those samples to only use in the transformation and output')
  parser.add_option('-u', '--includelocus', dest='includelocus', metavar='FILE',
                    help='Include list for those loci to only use in the transformation and output')

  parser.add_option('-x', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='Exclude a list of samples from the transformation and output')
  parser.add_option('-e', '--excludelocus', dest='excludelocus', metavar='FILE',
                    help='Exclude a list of loci from the transformation and output')

  parser.add_option('-m', '--renamesamples', dest='renamesamples', metavar='FILE',
                    help='A list of labels to rename samples')
  parser.add_option('-r', '--renamelocus', dest='renamelocus', metavar='FILE',
                    help='A list of labels to rename loci')

  parser.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows in a triple representation considered to N for testing purposes')

  # parser.add_option('-l', '--samplelimit', dest='samplelimit', metavar='N', type='int', default=0,
  #                   help='Limit the number of samples considered to N for testing purposes (default=0 for unlimited)')
  # parser.add_option('-L', '--locuslimit', dest='locuslimit', metavar='N', type='int', default=0,
  #                   help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  if not options.informat:
    print >> sys.stderr, 'Error: Must specify input data format'
    return
  if not options.outformat:
    print >> sys.stderr, 'Error: Must specify output data format'
    return

  if options.informat=='ldat':
    genos = load_genomatrix(args[0],options.limit)
    genos = transform_genomatrix_by_locus(genos,options)
    write_from_ldat(args[1], genos, options.outformat)
  elif options.informat=='sdat':
    genos = load_genomatrix(args[0],options.limit)
    genos = transform_genomatrix_by_sample(genos,options)
    write_from_sdat(args[1], genos, options.outformat)
  elif options.informat=='trip':
    genos = load_genotriples(args[0],options.limit)
    genos = transform_genotriples(genos,options)
    write_from_trip(args[1], genos, options.outformat)
  else:
    raise NotImplementedError,'Input file format "%s" is not supported' % options.informat


if __name__ == '__main__':
  main()
