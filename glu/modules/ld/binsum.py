# -*- coding: utf-8 -*-

__gluindex__  = True
__program__   = 'TagZilla binsum'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Generate summaries from multiple TagZilla output files'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.imerge          import imerge
from   glu.lib.fileutils       import autofile, hyphen, list_reader

from   glu.modules.ld.tagzilla import BinInfo, NullBinInfo, locus_result_sequence, bin_qualifier


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', nargs='+', help='Input genotype file')

  parser.add_argument('--license', action='store_true',
                          help="show program's copyright and license terms and exit")
  parser.add_argument('-b', '--summary', dest='sumfile', metavar='FILE', default='-',
                          help="Output summary tables FILE (default='-' for standard out)")
  parser.add_argument('-B', '--bininfo', metavar='FILE',
                          help='Output summary information about each bin to FILE')
  parser.add_argument('-H', '--histomax', metavar='N', type=int, default=10,
                          help='Largest bin size output in summary histogram output (default=10)')
  parser.add_argument('-s', '--subset', metavar='FILE', default='',
                          help='File containing loci to filter bins.  Only bins that contain one or more of these loci will be included in the output.')
  parser.add_argument('-t', '--targetbins', metavar='N', type=int, default=0,
                          help='Stop when N bins have been selected (default=0 for unlimited)')
  parser.add_argument('-T', '--targetloci', metavar='N', type=int, default=0,
                          help='Stop when N loci have been tagged (default=0 for unlimited)')

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  infofile = None
  if options.bininfo:
    infofile = autofile(hyphen(options.bininfo,sys.stdout), 'w')

  if options.bininfo or options.sumfile:
    bininfo = BinInfo(infofile, options.histomax+1)
  else:
    bininfo = NullBinInfo()

  sumfile = autofile(hyphen(options.sumfile,sys.stdout), 'w')

  if [infofile,sumfile].count(sys.stdout) > 1:
    sys.stderr.write('ERROR: More than one output file directed to standard out.\n')
    return

  subset  = set()
  exclude = set()

  if options.subset:
    subset = set(list_reader(options.subset))

  locusmap = {}
  results = [ locus_result_sequence(filename, locusmap, exclude) for filename in options.genotypes ]

  # FIXME: Make generic and use a better high-water mark
  # Unix-specific method to determine if we'll overflow the maximum number
  # of files descriptors allowed
  try:
    import resource
    slimit,hlimit = resource.getrlimit(resource.RLIMIT_NOFILE)
    if len(options.genotypes) > slimit*3/4:
      results = [ list(r) for r in results ]
  except ImportError:
    pass

  results = imerge(results)

  popdtags = {}
  populations = set()
  binned_loci = 0
  binnum = None
  for population,bin in results:
    if subset and not subset.intersection(bin):
      continue

    if binnum != bin.binnum:
      binnum = bin.binnum
      popdtags[bin.disposition] = popdtags.get(bin.disposition,0) + 1

    populations.add(population)

    bin_qualifier(bin, binned_loci, options)

    binned_loci += len(bin)
    bininfo.emit_bin(bin, locusmap, exclude, population)

  # Emit useful bin summary table
  for population in sorted(populations):
    bininfo.emit_summary(sumfile, population)

  if len(populations) > 1:
    bininfo.emit_multipop_summary(sumfile,popdtags)


if __name__ == '__main__':
  main()
