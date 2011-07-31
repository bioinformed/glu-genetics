# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Generate a matrix of pairwise LD values'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   operator                import itemgetter
from   itertools               import chain, izip

from   glu.lib.fileutils       import table_writer
from   glu.lib.xtab            import xtab
from   glu.lib.genolib         import geno_options


from   glu.modules.ld.tagzilla import generate_ldpairs


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', nargs='+', help='Input genotype file')

  inputgroup = parser.add_argument_group('Input options')

  geno_options(inputgroup,input=True,filter=True)

  inputgroup.add_argument('-S', '--ldsubset', metavar='FILE', default='',
                          help='File containing loci within the region these loci LD will be analyzed (see -d/--maxdist)')
  inputgroup.add_argument('-R', '--range', metavar='S-E,...', default='',
                          help='Ranges of genomic locations to analyze, specified as a comma separated list of start and '
                               'end coordinates "S-E".  If either S or E is not specified, then the ranges are assumed '
                               'to be open.  The end coordinate is exclusive and not included in the range.')

  outputgroup = parser.add_argument_group('Output options')

  outputgroup.add_argument('-o', '--output',   metavar='FILE',   default= '-',
                    help='Output file for formatted data')
  outputgroup.add_argument('-M', '--measure',  default='r2',
                    help="Measure of LD: r2 (default) or D'")

  genoldgroup = parser.add_argument_group('Genotype and LD estimation options')

  genoldgroup.add_argument('-a', '--minmaf', dest='maf', metavar='FREQ', type=float, default=0.05,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_argument('-c', '--mincompletion', metavar='N', default=0, type=int,
                          help='Drop loci with less than N valid genotypes. Default=0')
  genoldgroup.add_argument(      '--mincompletionrate', metavar='N', default=0, type=float,
                          help='Drop loci with completion rate less than N (0-1). Default=0')
  genoldgroup.add_argument('-m', '--maxdist', metavar='D', type=int, default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_argument('-P', '--hwp', metavar='p', default=None, type=float,
                          help='Filter out loci that fail to meet a minimum significance level (pvalue) for a '
                               'test Hardy-Weinberg proportion (no default)')

  bingroup = parser.add_argument_group('LD threshold options')

  bingroup.add_argument('-d', '--dthreshold', dest='d', metavar='DPRIME', type=float, default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_argument('-r', '--rthreshold', dest='r', metavar='N', type=float, default=0,
                          help='Minimum r-squared threshold to output (default=0)')

  return parser


def merge(i,j,x):
  return x[0] if x else ''


def main():
  parser  = option_parser()
  options = parser.parse_args()

  out = table_writer(options.output,hyphen=sys.stdout)

  if options.measure.lower() == 'r2':
    col = 2
  elif options.measure.lower() == "d'":
    col = 3
  else:
    raise ValueError('Unknown or unsupported LD measure specified: %s' % options.measure)

  locusmap = {}
  ldpairs  = generate_ldpairs(options, locusmap, set(), None, None)
  ldpairs  = chain.from_iterable(ldpairs)

  columns,rows,data = xtab(ldpairs,itemgetter(1),itemgetter(0),itemgetter(col),merge)
  out.writerow(['']+columns)

  for label,row in izip(rows,data):
    out.writerow([label]+row)


if __name__=='__main__':
  main()
