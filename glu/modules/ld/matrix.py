# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Generate a matrix of pairwise LD values'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import optparse

from   operator                import itemgetter
from   itertools               import chain, izip

from   glu.lib.fileutils       import table_writer
from   glu.lib.xtab            import xtab
from   glu.lib.genolib         import geno_options


from   glu.modules.ld.tagzilla import check_option01, generate_ldpairs


def option_parser():
  usage = 'usage: %prog [options] genotypes...'
  parser = optparse.OptionParser(usage=usage)

  inputgroup = optparse.OptionGroup(parser, 'Input options')

  geno_options(inputgroup,input=True,filter=True)

  inputgroup.add_option('-S', '--ldsubset', dest='ldsubset', metavar='FILE', default='',
                          help='File containing loci within the region these loci LD will be analyzed (see -d/--maxdist)')
  inputgroup.add_option('-R', '--range', dest='range', metavar='S-E,...', default='',
                          help='Ranges of genomic locations to analyze, specified as a comma separated list of start and '
                               'end coordinates "S-E".  If either S or E is not specified, then the ranges are assumed '
                               'to be open.  The end coordinate is exclusive and not included in the range.')

  outputgroup = optparse.OptionGroup(parser, 'Output options')

  outputgroup.add_option('-o', '--output',   dest='output',   metavar='FILE',   default= '-',
                    help='Output file for formatted data')
  outputgroup.add_option('-M', '--measure',  dest='measure', default='r2',
                    help="Measure of LD: r2 (default) or D'")

  genoldgroup = optparse.OptionGroup(parser, 'Genotype and LD estimation options')

  genoldgroup.add_option('-a', '--minmaf', dest='maf', metavar='FREQ', type='float', default=0.05,
                          action='callback', callback=check_option01,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_option('-c', '--mincompletion', dest='mincompletion', metavar='N', default=0, type='int',
                          help='Drop loci with less than N valid genotypes. Default=0')
  genoldgroup.add_option(      '--mincompletionrate', dest='mincompletionrate', metavar='N', default=0, type='float',
                          action='callback', callback=check_option01,
                          help='Drop loci with completion rate less than N (0-1). Default=0')
  genoldgroup.add_option('-m', '--maxdist', dest='maxdist', metavar='D', type='int', default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_option('-P', '--hwp', dest='hwp', metavar='p', default=None, type='float',
                          action='callback', callback=check_option01,
                          help='Filter out loci that fail to meet a minimum significance level (pvalue) for a '
                               'test Hardy-Weinberg proportion (no default)')

  bingroup = optparse.OptionGroup(parser, 'LD threshold options')

  bingroup.add_option('-d', '--dthreshold', dest='d', metavar='DPRIME', type='float', default=0.,
                          action='callback', callback=check_option01,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_option('-r', '--rthreshold', dest='r', metavar='N', type='float', default=0,
                          action='callback', callback=check_option01,
                          help='Minimum r-squared threshold to output (default=0)')

  parser.add_option_group(inputgroup)
  parser.add_option_group(outputgroup)
  parser.add_option_group(genoldgroup)
  parser.add_option_group(bingroup)

  return parser


def merge(i,j,x):
  return x[0] if x else ''


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    sys.exit(2)

  out = table_writer(options.output,hyphen=sys.stdout)

  if options.measure.lower() == 'r2':
    col = 2
  elif options.measure.lower() == "d'":
    col = 3
  else:
    raise ValueError('Unknown or unsupported LD measure specified: %s' % options.measure)

  args = [(options,arg) for arg in args]
  locusmap = {}
  ldpairs  = generate_ldpairs(args, locusmap, set(), None, None, options)
  ldpairs  = chain.from_iterable(ldpairs)

  columns,rows,data = xtab(ldpairs,itemgetter(1),itemgetter(0),itemgetter(col),merge)
  out.writerow(['']+columns)

  for label,row in izip(rows,data):
    out.writerow([label]+row)


if __name__=='__main__':
  main()
