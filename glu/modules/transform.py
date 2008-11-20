# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Performs various transformation on a genotype files'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.genolib.merge   import get_genomerger, output_merge_statistics
from   glu.lib.genolib.streams import GenoTransform
from   glu.lib.genolib.io      import transform_files, geno_options


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotypes...'
  parser = optparse.OptionParser(usage=usage)

  ioopts = optparse.OptionGroup(parser, 'Input/Output Options')

  ioopts.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output transformed data to FILE(default is "-" for standard out)')

  geno_options(ioopts,input=True,output=True)

  mopts = optparse.OptionGroup(parser, 'Genotype Merging and Reporting')
  geno_options(mopts,merge=True)

  filter = optparse.OptionGroup(parser, 'Filtering and Renaming')
  geno_options(filter,filter=True)

  trans = optparse.OptionGroup(parser, 'Transformations')
  geno_options(trans,transform=True)

  parser.add_option_group(ioopts)
  parser.add_option_group(mopts)
  parser.add_option_group(filter)
  parser.add_option_group(trans)

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    return

  if not options.output:
    parser.print_help(sys.stderr)
    sys.stderr.write('Error: Must specify output file\n')
    return

  infiles = sorted(set(args))
  merger  = get_genomerger(options.merge, bool(options.samplemerge or options.locusmerge))

  transform = GenoTransform.from_options(options)

  transform_files(infiles, options.informat,  options.ingenorepr,
                  options.output, options.outformat, options.outgenorepr,
                  transform, mergefunc=merger, genome=options.loci, phenome=options.pedigree,
                  inhyphen=sys.stdin,outhyphen=sys.stdout)

  output_merge_statistics(merger, options.samplemerge, options.locusmerge)


if __name__ == '__main__':
  main()
