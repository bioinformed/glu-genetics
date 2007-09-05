# -*- coding: utf-8 -*-
'''
File:          ginfo.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       September 5, 2007

Abstract:      This utility script can extract the meta data information if any from the genotype data file

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id: $
'''

import sys
import csv
from itertools import izip
from glu.lib.fileutils  import autofile
from glu.lib.genolib.io import load_genostream


def emit(filename,rows):
  out = autofile(filename,'w')
  for row in rows:
    out.write('%s\n' % row)


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] genofile'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--loci',   dest='loci',    metavar='FILE',
                     help='Output the list of loci to FILE')
  parser.add_option('--samples',dest='samples', metavar='FILE',
                     help='Output the list of samples to FILE')
  parser.add_option('--models',  dest='models', metavar='FILE',
                     help='Output the genotype models (alleles for each locus) to FILE')
  parser.add_option('--format',  dest='format', type='str',
                     help='The input genotype file format, possible values=hapmap,ldat,sdat,lbat,sbat,tbat,trip,genotriple')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) < 1:
    parser.print_help()
    return

  genos  = load_genostream(args[0],options.format)

  if genos.samples is not None:
    print 'sample count=', len(genos.samples)
    if options.samples:
      emit(options.samples,genos.samples)

  if genos.loci is not None:
    print 'locus count=', len(genos.columns)
    if options.loci:
      emit(options.loci,genos.loci)

  if genos.models is not None:
    print 'model count=', len(genos.models)
    if options.models:
      rows = []
      for locus,model in izip(genos.loci,genos.models):
        row = [locus]
        if None in model.alleles:
          model.alleles.remove(None)
        row.extend(model.alleles)
        rows.append('\t'.join(row))
      emit(options.models,rows)


if __name__=='__main__':
  main()
