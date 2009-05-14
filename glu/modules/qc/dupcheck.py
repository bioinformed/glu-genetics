# -*- coding: utf-8 -*-

from   __future__ import division

__gluindex__  = True
__abstract__  = 'Detect duplicate samples based on an all pairs-comparison of genotype concordance rates'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools                 import chain

from   glu.lib.utils             import pair_generator
from   glu.lib.fileutils         import table_reader,table_writer
from   glu.lib.union_find        import union_find

from   glu.lib.genolib.io        import load_genostream,geno_options
from   glu.lib.genolib.genoarray import genoarray_concordance


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotypes'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True)

  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='Mapping from sample identifier to subject identifier')
  parser.add_option('--checkexp', dest='checkexp', action='store_true',
                    help='Check only expected duplicate pairs')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-t', '--threshold', dest='threshold', metavar='N', type='float', default=0.80,
                    help='Minimum proportion genotype concordance threshold for duplicates (default=0.8)')
  parser.add_option('-m', '--mingenos', dest='mingenos', metavar='N', type='int', default=20,
                    help='Minimum number of concordant genotypes to be considered informative (default=20)')

  return parser


def squeeze_pairs(sets,universe=None):
  '''
  Squeeze singleton sets after optionally finding the intersection with a
  universal set
  '''
  if universe is not None:
    if not isinstance(universe, (set,dict)):
      universe = set(universe)
    sets = (s&universe for s in sets)
  return [ s for s in sets if len(s)>1 ]


def expected_pairs(genos, dupset):
  '''
  Form pairs by materializing only samples that participate in expected
  dupsets of size>1 and generating pairs from within members of each of the
  resulting sets.
  '''
  sets = squeeze_pairs(dupset.sets(),genos.samples)

  samples = set()
  for s in sets:
    samples.update(s)

  genos = genos.transformed(includesamples=samples).as_sdat().materialize()

  sets = squeeze_pairs(sets,genos.samples)
  genos = dict(genos)
  pairs = [ pair_generator([ (s,genos[s]) for s in dset ]) for dset in sets ]

  return chain(*pairs)


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  threshold = options.threshold

  if threshold < 0 or threshold > 1:
    raise ValueError('Invalid threshold %s given.  Must be between 0 and 1.' % threshold)

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  expected_dupset = union_find()
  observed_dupset = union_find()

  if options.duplicates:
    for dupset in table_reader(options.duplicates):
      dupset = [ d for d in dupset if d ]
      for dup in dupset:
        expected_dupset.union(dupset[0],dup)

  genos = load_genostream(args[0],format=options.informat,genorepr=options.ingenorepr,
                                   genome=options.loci,phenome=options.pedigree,
                                   transform=options)

  if options.checkexp:
    pairs = expected_pairs(genos, expected_dupset)
  else:
    pairs = pair_generator(genos.as_sdat().materialize())

  status = { (True, True ): ['EXPECTED',  'CONCORDANT'],
             (True, False): ['EXPECTED',  'DISCORDANT'],
             (False,True ): ['UNEXPECTED','CONCORDANT'] }

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE1','SAMPLE2','CONCORDANT_GENOTYPES','COMPARISONS','CONCORDANCE_RATE',
                'EXPECTED_DUPLICATE','OBSERVED_DUPLICATE'])

  # N.B. It is not possible to output duplicates incrementally and enumerate
  # duplicate sets.  This is because sets can merge transitively as
  # pairs are observed, so merges among existing sets can occur at any time.
  # So unless we choose to buffer output, enumeration of duplicate sets is
  # not currently supported.
  for (sample1,genos1),(sample2,genos2) in pairs:
    matches,comparisons = genoarray_concordance(genos1, genos2)

    obs_dup = matches>=options.mingenos and comparisons and matches/comparisons >= threshold
    exp_dup = (sample1,sample2) in expected_dupset

    if not obs_dup and not exp_dup:
      continue

    if sample1>sample2:
      sample1,sample2=sample2,sample1

    if obs_dup:
      observed_dupset.union(sample1,sample2)

    rate = matches/comparisons if comparisons else ''

    out.writerow([sample1,sample2,matches,comparisons,rate]+status[exp_dup,obs_dup])


if __name__ == '__main__':
  main()
