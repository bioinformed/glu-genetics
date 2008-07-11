# -*- coding: utf-8 -*-

from   __future__ import division

__gluindex__  = True
__abstract__  = 'Detect duplicate samples based on an all pairs-comparison of genotype concordance rates'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.utils             import pair_generator
from   glu.lib.fileutils         import table_writer
from   glu.lib.union_find        import union_find

from   glu.lib.genolib.io        import load_genostream,geno_options
from   glu.lib.genolib.merge     import get_genomerger
from   glu.lib.genolib.genoarray import genoarray_concordance


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True)

  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='Mapping from sample identifier to subject identifier')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity shared between two individuals (default=85)')
  parser.add_option('-m', '--mingenos', '--mincount', dest='mingenos', metavar='N', type='int', default=20,
                    help='Minimum concordant genotypes to be considered informative for duplicate checking')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  threshold = float(options.threshold)/100

  if len(args) != 1:
    parser.print_help()
    return

  expected_dupset = union_find()
  observed_dupset = union_find()

  if options.duplicates:
    for dupset in load_table(options.duplicates):
      dupset = [ d for d in dupset if d ]
      for dup in dupset:
        expected_dupset.union(dupset[0],dup)

  genos  = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                   genome=options.loci,phenome=options.pedigree)
  genos  = genos.as_sdat().materialize()

  status = { (True, True ): ['EXPECTED',  'CONCORDANT'],
             (True, False): ['EXPECTED',  'DISCORDANT'],
             (False,True ): ['UNEXPECTED','CONCORDANT'] }

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE1','SAMPLE2','CONCORDANT_GENOTYPES','COMPARISONS',
                'EXPECTED_DUPLICATE','OBSERVED_DUPLICATE'])

  # N.B. It is not possible to output duplicates incrementally and enumerate
  # duplicate sets.  This is because sets can merge transitively as
  # pairs are observed, so merges among existing sets can occur at any time.
  # So unless we choose to buffer output, enumeration of duplicate sets is
  # not currently supported.
  for (sample1,genos1),(sample2,genos2) in pair_generator(genos):
    matches,comparisons = genoarray_concordance(genos1, genos2)

    obs_dup = matches>=options.mingenos and comparisons and float(matches)/comparisons >= threshold
    exp_dup = (sample1,sample2) in expected_dupset

    if not obs_dup and not exp_dup:
      continue

    if sample1>sample2:
      sample1,sample2=sample2,sample1

    if obs_dup:
      observed_dupset.union(sample1,sample2)

    out.writerow([sample1,sample2,matches,comparisons]+status[exp_dup,obs_dup])


if __name__ == '__main__':
  main()
