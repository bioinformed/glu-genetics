# -*- coding: utf-8 -*-

__gluindex__  = True
__program__   = 'TagZilla coverage'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Evaluate maximum coverage of a set of tags'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.fileutils       import list_reader, table_writer
from   glu.lib.genolib         import geno_options

from   glu.modules.ld.tagzilla import epsilon, sfloat, \
                                      build_design_score, generate_ldpairs


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('tagfile',              help='Delimited or tabular file containing tag marker names')
  parser.add_argument('genotypes', nargs='+', help='Genotype file(s)')

  inputgroup = parser.add_argument_group('Input options')

  geno_options(inputgroup,input=True,filter=True)

  inputgroup.add_argument('-e', '--excludetag', dest='exclude', metavar='FILE', default='',
                          help='File containing loci that are excluded from covering another locus')
  inputgroup.add_argument('-R', '--range', metavar='S-E,...', default='',
                          help='Ranges of genomic locations to analyze, specified as a comma separated list of start and '
                               'end coordinates "S-E".  If either S or E is not specified, then the ranges are assumed '
                               'to be open.  The end coordinate is exclusive and not included in the range.')
  inputgroup.add_argument('-D', '--designscores', metavar='FILE', type=str, action='append',
                          help='Read in design scores or other weights to use as criteria to choose the optimal tag for each bin')
  inputgroup.add_argument('--designdefault', metavar='N', type=float, default=0,
                          help='Default design score for any locus not found in a design file')
  inputgroup.add_argument('-L', '--limit', metavar='N', type=int, default=0,
                          help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')

  outputgroup = parser.add_argument_group('Output options')

  outputgroup.add_argument('-o', '--output', metavar='FILE', default='-',
                          help="Output tabular LD information for bins to FILE ('-' for standard out)")

  genoldgroup = parser.add_argument_group('Genotype and LD estimation options')

  genoldgroup.add_argument('-a', '--minmaf', dest='maf', metavar='FREQ', type=float, default=0.05,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_argument('-A', '--minobmaf', dest='obmaf', metavar='FREQ', type=float, default=None,
                          help='Minimum minor allele frequency (MAF) for obligate tags (defaults to -a/--minmaf)')
  genoldgroup.add_argument('-c', '--mincompletion', metavar='N', default=0, type=int,
                          help='Drop loci with less than N valid genotypes. Default=0')
  genoldgroup.add_argument(      '--mincompletionrate', metavar='N', default=0, type=float,
                          help='Drop loci with completion rate less than N (0-1). Default=0')
  genoldgroup.add_argument('-m', '--maxdist', metavar='D', type=int, default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_argument('-P', '--hwp', metavar='p', default=None, type=float,
                          help='Filter out loci that fail to meet a minimum significance level (pvalue) for a '
                               'test Hardy-Weinberg proportion (no default)')

  bingroup = parser.add_argument_group('Binning options')

  bingroup.add_argument('-d', '--dthreshold', dest='d', metavar='DPRIME', type=float, default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_argument('-r', '--rthreshold', dest='r', metavar='N', type=float, default=0,
                          help='Minimum r-squared threshold to output (default=0)')

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  exclude = set()

  tags = set(list_reader(options.tagfile))

  if options.exclude:
    exclude = set(list_reader(options.exclude))
  elif options.designscores:
    exclude = set()

  designscores = build_design_score(options.designscores)

  # Trim undesignable tags
  tags -= exclude
  tags -= set(t for t,score in designscores.iteritems() if score < epsilon)
  if options.designscores and options.designdefault <= epsilon:
    tags -= set(t for t in tags if t not in designscores)

  locusmap = {}
  ldpairs = generate_ldpairs(options, locusmap, set(), None, tags)

  loci    = set()
  besttag = dict( (tag,(tag,1)) for tag in tags )
  missing = '',-1
  for pairs in ldpairs:
    for lname1,lname2,r2,dprime in pairs:
      if lname1==lname2:
        loci.add(lname1)

      for l1,l2 in (lname1,lname2),(lname2,lname1):
        if l2 in tags:
          best_locus,best_r2 = besttag.get(l1,missing)
          if l2 not in exclude and r2 > best_r2:
            besttag[l1] = l2,r2
          break

  outfile = table_writer(options.output, hyphen=sys.stdout)
  outfile.writerow(['LNAME','TAG','BEST RSQUARED'])
  missing = '',sfloat(0)
  for lname in loci:
    tag,r2 = missing
    if lname in besttag:
      tag,r2 = besttag[lname]
      r2 = sfloat(r2)
    outfile.writerow([lname,tag,r2])

  for lname in locusmap:
    if lname not in loci:
      outfile.writerow([lname,'',sfloat(0)])


if __name__ == '__main__':
  main()
