# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Wrapper to simplift conversion of a GDAT file into any supported GLU genotype file'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   glu.lib.genolib.io        import load_genostream, save_genostream, geno_options


def option_parser():
  import optparse

  usage = 'usage: %prog [options] infile.gdat'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True,transform=True,output=True)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output genotype file name')
  #parser.add_option('-s', '--targetstrand', dest='targetstrand', metavar='T', default='customer',
  #                  help='Target strand based on Illumina manifest file: ab, top, bottom, forward, '
  #                       'reverse, customer (default), anticustomer, design, antidesign')
  parser.add_option('-t', '--gcthreshold', dest='gcthreshold', type='float', metavar='N', default=0,
                    help='Genotypes with GC score less than N set to missing')
  parser.add_option('--samplestats', dest='samplestats', metavar='FILE',
                    help='Output per sample average GC statistics to FILE')
  parser.add_option('--locusstats',  dest='locusstats',  metavar='FILE',
                    help='Output per locus average GC statistics to FILE')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  if options.informat and options.informat!='gdat':
    raise ValueError('Input must be in gdat format')

  genos = load_genostream(args[0],format='gdat',genorepr=options.ingenorepr,
                                  genome=options.loci,phenome=options.pedigree,
                                  transform=options, hyphen=sys.stdin,
                                  gcthreshold=options.gcthreshold,
                                  samplestats=options.samplestats,
                                  locusstats=options.locusstats)

  save_genostream(options.output,genos,format=options.outformat,genorepr=options.outgenorepr,
                                 hyphen=sys.stdout)


if __name__ == '__main__':
  main()
