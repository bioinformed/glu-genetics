# -*- coding: utf-8 -*-

__gluindex__  = True
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Find surrogates SNPs'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import heapq

from   collections             import defaultdict

from   glu.lib.fileutils       import list_reader, table_writer
from   glu.lib.genolib         import geno_options

from   glu.modules.ld.tagzilla import epsilon, sfloat, \
                                      build_design_score, generate_ldpairs


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('needles',              help='Tabular or delimited file of markers')
  parser.add_argument('genotypes', nargs='+', help='Input genotype file(s)')

  inputgroup = parser.add_argument_group('Input options')

  geno_options(inputgroup,input=True,filter=True)

  inputgroup.add_argument('--haystack', metavar='FILE',
                          help='List of SNPs eligable to be surrogates (minus any in needle.lst if --indirect is specifed)')
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

  bingroup.add_argument('--indirect', action='store_true',
                          help='Allow only indirect surrogates')
  bingroup.add_argument('-s', '--maxsurrogates', metavar='N', type=int, default=0,
                          help='Maximum number of surrogates (default=0 for unlimited)')
  bingroup.add_argument('-d', '--dthreshold', dest='d', metavar='DPRIME', type=float, default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_argument('-r', '--rthreshold', dest='r', metavar='N', type=float, default=0.15,
                          help='Minimum r-squared threshold to output (default=0.15)')

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  exclude = set()
  if options.designscores:
    designscores = build_design_score(options.designscores)
    exclude.update(lname for lname,d in designscores.iteritems() if d <= epsilon)

  seen     = set()
  needle   = set(list_reader(options.needles))

  # ldsubset=needle
  ldpairs = generate_ldpairs(options, {}, set(), None, needle)

  # FIXME: Obtain default haystack from ldpairs
  if options.haystack:
    haystack = set(list_reader(options.haystack))-exclude

    if options.indirect:
      haystack = haystack-needle
  else:
    direct = needle if options.indirect else set()

  heappushpop    = heapq.heappushpop
  heappush       = heapq.heappush
  maxsurrogates  = options.maxsurrogates
  best_surrogate = defaultdict(list)

  def update_surrogates(lname,surrogate,r2):
    surrogates = best_surrogate[lname]

    # FIXME: very expensive 2-tuple creation to add priorities.  May be
    # better to devise a fixed point encoding and use least-significant bits
    # for priority
    if lname==surrogate:
      priority = r2,2
    elif surrogate in needle:
      priority = r2,1
    else:
      priority = r2,0

    if not maxsurrogates or len(surrogates)<maxsurrogates:
      heappush(surrogates, (priority,surrogate) )
    elif priority>surrogates[0][0]:
      heappushpop(surrogates, (priority,surrogate) )

  # FIXME: this loop is duplicated only because haystack-direct is not
  #        always available?  LD calculation typically requires
  #        pre-materialization, so the universal haystack should always be
  #        known, though perhaps not accessible outside the internal
  #        tagzilla code?
  if haystack is None:
    for pairs in ldpairs:
      for lname1,lname2,r2,dprime in pairs:
        # FIXME: May be able to add to lname1==lname2 branch and reduce to single set op
        seen.add(lname1)
        seen.add(lname2)

        # FIXME: should be able to use lname1 is lname2?
        if lname1==lname2:
          if lname1 in needle and lname1 not in direct:
            update_surrogates(lname1,lname1,r2)
        else:
          if lname1 in needle and lname2 not in direct:
            update_surrogates(lname1,lname2,r2)
          if lname2 in needle and lname1 not in direct:
            update_surrogates(lname2,lname1,r2)
  else:
    for pairs in ldpairs:
      for lname1,lname2,r2,dprime in pairs:
        # FIXME: May be able to add to lname1==lname2 branch and reduce to single set op
        seen.add(lname1)
        seen.add(lname2)

        # FIXME: should be able to use lname1 is lname2?
        if lname1==lname2:
          if lname1 in needle and lname1 in haystack:
            update_surrogates(lname1,lname1,r2)
        else:
          if lname1 in needle and lname2 in haystack:
            update_surrogates(lname1,lname2,r2)
          if lname2 in needle and lname1 in haystack:
            update_surrogates(lname2,lname1,r2)

  outfile = table_writer(options.output, hyphen=sys.stdout)
  outfile.writerow(['LNAME','RANK','SURROGATE','RSQUARED','REASON'])

  priorities = ['AUGMENT','PANEL','SELF']

  for lname in sorted(needle):
    surrogates = best_surrogate.get(lname)

    if not surrogates:
      reason = 'NO SURROGATE' if lname in seen else 'NO DATA'
      outfile.writerow([lname,'','','',reason])
    else:
      surrogates.sort(reverse=True)
      for i,((r2,prior),sname) in enumerate(surrogates):
        outfile.writerow([lname,i+1,sname,sfloat(r2),priorities[prior]])


if __name__ == '__main__':
  main()
