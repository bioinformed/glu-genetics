# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Convert GLU GDAT file into any other supported GLU genotype data format'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   glu.lib.genolib.io        import load_genostream, save_genostream, geno_options


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('input', help='GDAT file')

  geno_options(parser,input=True,filter=True,transform=True,output=True)

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output genotype file name')
  #parser.add_argument('-s', '--targetstrand', metavar='T', default='customer',
  #                  help='Target strand based on Illumina manifest file: ab, top, bottom, forward, '
  #                       'reverse, customer (default), anticustomer, design, antidesign')
  parser.add_argument('-t', '--gcthreshold', type=float, metavar='N', default=0,
                    help='Genotypes with GC score less than N set to missing')
  parser.add_argument('--samplestats', metavar='FILE',
                    help='Output per sample average GC statistics to FILE')
  parser.add_argument('--locusstats', metavar='FILE',
                    help='Output per locus average GC statistics to FILE')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if options.informat and options.informat!='gdat':
    raise ValueError('Input must be in gdat format')

  genos = load_genostream(options.input,format='gdat',genorepr=options.ingenorepr,
                                        genome=options.loci, phenome=options.pedigree,
                                        transform=options, hyphen=sys.stdin,
                                        gcthreshold=options.gcthreshold,
                                        samplestats=options.samplestats,
                                        locusstats=options.locusstats)

  save_genostream(options.output,genos,format=options.outformat,genorepr=options.outgenorepr,
                                       hyphen=sys.stdout)


if __name__ == '__main__':
  main()
