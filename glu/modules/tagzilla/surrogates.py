# -*- coding: utf-8 -*-
'''
File:          surrogates.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       November 8, 2005

Abstract:      Find the best surrogate SNP for any that fail design

Compatibility: Python 2.5 and above

Requires:      No external dependencies, yet...

Revision:      $Id$
'''

__program__   = 'TagZilla surrogates'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

__accelerators__ = ['tagzillac']

from tagzilla import *


def option_parser():
  usage = 'usage: %prog [options] tagfile genofile...'
  parser = TagZillaOptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                        help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  parser.add_option('--profile', dest='profile', metavar='P', help=optparse.SUPPRESS_HELP)

  inputgroup = optparse.OptionGroup(parser, 'Input options')

  inputgroup.add_option('-f', '--format', dest='format', metavar='NAME', default='',
                          help='Format for genotype/pedigree or ld input data.  Values: hapmap (default), linkage, festa, prettybase, raw.')
  inputgroup.add_option(      '--pedformat', dest='pedformat', metavar='NAME', default='',
                          help='Format for pedigree data.  Values: hapmap or linkage.  Defaults to hapmap when '
                               'reading HapMap files and linkage format otherwise.')
  inputgroup.add_option('-l', '--loci', dest='loci', metavar='FILE',
                          help='Locus description file for input in Linkage format')
  inputgroup.add_option('-p', '--pedfile', dest='pedfile', metavar='FILE', action='append',
                          help='Pedigree file for HapMap or PrettyBase data files (optional)')
  inputgroup.add_option('-e', '--excludetag', dest='exclude', metavar='FILE', default='',
                          help='File containing loci that are excluded from being a tag')
  inputgroup.add_option('-s', '--subset', dest='subset', metavar='FILE', default='',
                          help='File containing loci that define the subset to be analyzed of the loci that are read')
  inputgroup.add_option('-R', '--range', dest='range', metavar='S-E,...', default='',
                          help='Ranges of genomic locations to analyze, specified as a comma seperated list of start and '
                               'end coordinates "S-E".  If either S or E is not specified, then the ranges are assumed '
                               'to be open.  The end coordinate is exclusive and not included in the range.')
  inputgroup.add_option('-D', '--designscores', dest='designscores', metavar='FILE', type='str', action='append',
                          help='Read in design scores or other weights to use as criteria to choose the optimal tag for each bin')
  inputgroup.add_option('--designdefault', dest='designdefault', metavar='N', type='float', default=0,
                          help='Default design score for any locus not found in a design file')
  inputgroup.add_option('-L', '--limit', dest='limit', metavar='N', type='int', default=0,
                          help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')

  outputgroup = optparse.OptionGroup(parser, 'Output options')

  outputgroup.add_option('-o', '--output', dest='outfile', metavar='FILE', default='-',
                          help="Output tabular LD information for bins to FILE ('-' for standard out)")

  genoldgroup = optparse.OptionGroup(parser, 'Genotype and LD estimation options')

  genoldgroup.add_option('-a', '--minmaf', dest='maf', metavar='FREQ', type='float', default=0.05,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_option('-A', '--minobmaf', dest='obmaf', metavar='FREQ', type='float', default=None,
                          help='Minimum minor allele frequency (MAF) for obligate tags (defaults to -a/--minmaf)')
  genoldgroup.add_option('-c', '--mincompletion', dest='mincompletion', metavar='N', default=0, type='int',
                          help='Drop loci with less than N valid genotypes. Default=0')
  genoldgroup.add_option(      '--mincompletionrate', dest='mincompletionrate', metavar='N%', default=0, type='float',
                          help='Drop loci with completion rate less than N% (0-100). Default=0')
  genoldgroup.add_option('-m', '--maxdist', dest='maxdist', metavar='D', type='int', default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_option('-P', '--hwp', dest='hwp', metavar='p', default=None, type='float',
                          help='Filter out loci that fail to meet a minimum signficance level (pvalue) for a '
                               'test Hardy-Weignberg proportion (no default)')

  bingroup = optparse.OptionGroup(parser, 'Binning options')

  bingroup.add_option('-d', '--dthreshold', dest='d', metavar='DPRIME', type='float', default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_option('-r', '--rthreshold', dest='r', metavar='N', type='float', default=0,
                          help='Minimum r-squared threshold to output (default=0)')

  parser.add_option_group(inputgroup)
  parser.add_option_group(outputgroup)
  parser.add_option_group(genoldgroup)
  parser.add_option_group(bingroup)

  return parser


def surrogates(options,args):
  subset  = set()
  exclude = set()

  if options.subset:
    read_snp_list(options.subset, subset)

  if options.exclude:
    read_snp_list(options.exclude, exclude)

  designscores = build_design_score(options.designscores)
  if options.designdefault <= epsilon:
    ldsubset = set(lname for lname,d in designscores.iteritems() if d <= epsilon)
  else:
    ldsubset = set()

  locusmap = {}
  options.multipopulation = None
  ldpairs = generate_ldpairs(args, locusmap, set(), subset, ldsubset, options)

  missing = '',0
  best_surrogate = {}
  for pairs in ldpairs:
    for lname1,lname2,r2,dprime in pairs:
      d1=designscores.get(lname1,options.designdefault) > epsilon
      d2=designscores.get(lname2,options.designdefault) > epsilon

      if d1 and not d2 and lname2 in designscores:
        best_locus,best_r2 = best_surrogate.get(lname2,missing)
        if r2 > best_r2:
          best_surrogate[lname2] = lname1,r2
      elif not d1 and d2 and lname1 in designscores:
        best_locus,best_r2 = best_surrogate.get(lname1,missing)
        if r2 > best_r2:
          best_surrogate[lname1] = lname2,r2

  outfile = autofile(options.outfile, 'w', hyphen=sys.stdout)
  outfile.write('LNAME\tSURROGATE\tRSQUARED\n')
  for lname,(surrogate,r2) in best_surrogate.iteritems():
    r2 = sfloat(r2)
    outfile.write('%s\t%s\t%s\n' % (lname,surrogate,r2))


def main():
  launcher(surrogates, option_parser, **globals())


if __name__ == '__main__':
  main()
