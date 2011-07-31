# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Performs various transformation on a genotype files'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.genolib.merge     import get_genomerger, output_merge_statistics
from   glu.lib.genolib.transform import GenoTransform
from   glu.lib.genolib.io        import transform_files, geno_options


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', nargs='+', help='Input genotype file(s)')

  ioopts = parser.add_argument_group('Input/Output Options')

  ioopts.add_argument('-o', '--output', metavar='FILE', default='-', required=True,
                    help='Output transformed data to FILE (default is "-" for standard out)')

  geno_options(ioopts,input=True,output=True)

  mopts = parser.add_argument_group('Genotype Merging and Reporting')
  geno_options(mopts,merge=True)

  filter = parser.add_argument_group('Filtering')
  geno_options(filter,filter=True)

  trans = parser.add_argument_group('Transformations')
  geno_options(trans,transform=True)

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  infiles = sorted(set(options.genotypes))
  merger  = get_genomerger(options.merge, bool(options.samplemerge or options.locusmerge))

  transform = GenoTransform.from_options(options)

  transform_files(infiles, options.informat,  options.ingenorepr,
                  options.output, options.outformat, options.outgenorepr,
                  transform, mergefunc=merger, genome=options.loci, phenome=options.pedigree,
                  inhyphen=sys.stdin,outhyphen=sys.stdout)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__ == '__main__':
  main()
