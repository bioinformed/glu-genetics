# -*- coding: utf-8 -*-
'''
File:          ginfo.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       September 5, 2007

Abstract:      Extract available metadata from a GLU genotype file

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id: $
'''

import sys
import csv

from   itertools import izip

from   glu.lib.fileutils     import autofile, namefile, hyphen
from   glu.lib.genolib.io    import load_genostream
from   glu.lib.genolib.reprs import get_genorepr


def emit(filename,rows):
  outfile = autofile(hyphen(filename,sys.stdout),'w')
  csv.writer(outfile,dialect='tsv').writerows(rows)


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] genofile'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--format',  dest='format', type='str',
                     help='The input genotype file format, possible values=hapmap,ldat,sdat,lbat,sbat,tbat,trip,genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  parser.add_option('--loci',   dest='loci',    metavar='FILE',
                     help='Output the list of loci to FILE')
  parser.add_option('--samples',dest='samples', metavar='FILE',
                     help='Output the list of samples to FILE')
  parser.add_option('--models',  dest='models', metavar='FILE',
                     help='Output the genotype models (alleles for each locus) to FILE')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  infile   = hyphen(args[0],sys.stdin)
  genorepr = get_genorepr(options.genorepr)
  genos    = load_genostream(infile,options.format,genorepr=genorepr)

  out = autofile(hyphen(options.output,sys.stdout),'w')

  out.write('Filename    : %s\n' % namefile(infile))
  out.write('Format      : %s\n' % genos.format)

  if genos.samples is not None:
    out.write('sample count: %d\n' % len(genos.samples))
    if options.samples:
      emit(options.samples,([s] for s in genos.samples) )

  if genos.loci is not None:
    out.write('locus  count: %d\n' % len(genos.loci))
    if options.loci:
      emit(options.loci, ([l] for l in genos.loci))

  if genos.samples is not None and genos.loci is not None:
    out.write('model  count: %d\n' % len(set(genos.models)))
    if options.models:
      def _models():
        for locus,model in izip(genos.loci,genos.models):
          yield [locus, ','.join(sorted(model.alleles[1:]))]
      emit(options.models,_models())


if __name__=='__main__':
  main()
