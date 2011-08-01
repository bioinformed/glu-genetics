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
from   glu.lib.progressbar       import progress_loop

from   glu.lib.genolib.io        import load_genostream,geno_options
from   glu.lib.genolib.genoarray import genoarray_concordance


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-e', '--duplicates', metavar='FILE',
                    help='Mapping from sample identifier to subject identifier')
  parser.add_argument('--testpairs', metavar='FILE',
                    help='File containing a list of pairs to test')
  parser.add_argument('--testexp', action='store_true',
                    help='Check only expected duplicate pairs')
  parser.add_argument('-t', '--threshold', metavar='N', type=float, default=0.80,
                    help='Minimum proportion genotype concordance threshold for duplicates (default=0.8)')
  parser.add_argument('-m', '--mingenos', metavar='N', type=int, default=20,
                    help='Minimum number of concordant genotypes to be considered informative (default=20)')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output duplicate report')
  parser.add_argument('--dupout', metavar='FILE',
                    help='Output sets of observed duplicate samples')
  parser.add_argument('-P', '--progress', action='store_true',
                    help='Show analysis progress bar, if possible')

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


def flatten_dupsets(sets):
  '''
  Flatten dupsets into a single flat set
  '''
  flat = set()
  for s in sets:
    flat.update(s)
  return flat


def write_dupsets(filename, dupsets):
  '''
  Output a table of duplicates where each row contains all samples that met
  the concordance threshold
  '''
  out = table_writer(filename)
  for dset in dupsets.sets():
    if len(dset) > 1:
      out.writerow(sorted(dset))


def expected_pairs(genos, dupset):
  '''
  Form pairs by materializing only samples that participate in expected
  dupsets of size>1 and generating pairs from within members of each of the
  resulting sets.
  '''
  sets    = squeeze_pairs(dupset.sets(),genos.samples)
  samples = flatten_dupsets(sets)
  genos   = genos.transformed(includesamples=samples).as_sdat().materialize()
  sets    = squeeze_pairs(sets,genos.samples)
  genos   = dict(genos)
  n       = sum( len(dset)*(len(dset)-1)//2 for dset in sets )
  pairs   = [ pair_generator([ (s,genos[s]) for s in dset ]) for dset in sets ]

  return chain(*pairs),n


def file_pairs(genos, pairfile):
  '''
  Form pairs by generating only pairs of samples that appear in a file
  '''
  data = table_reader(pairfile)
  genos = dict(genos)

  # Materialize pairs, since the length is needed for the progress bar
  pairs = []

  # Create a local intern dictionary for sample names, since the list of
  # pairs can be very long
  samples = dict( (sample,sample) for sample in genos )

  for row in data:
    if len(row) < 2:
      continue

    # Intern sample names
    sample1 = samples.get(row[0])
    sample2 = samples.get(row[1])

    if sample1 and sample2:
      pairs.append( ((sample1,genos[sample1]),(sample2,genos[sample2])) )

  # FIXME: We may want to implement a progress bar with a spinner for
  # unknown length, since materializing here may be fatal to performance
  return pairs,len(pairs)


def main():
  parser  = option_parser()
  options = parser.parse_args()

  threshold = options.threshold

  if threshold < 0 or threshold > 1:
    raise ValueError('Invalid threshold %s given.  Must be between 0 and 1.' % threshold)

  expected_dupset = union_find()
  observed_dupset = union_find()

  if options.duplicates:
    for dupset in table_reader(options.duplicates):
      dupset = [ d for d in dupset if d ]
      for dup in dupset:
        expected_dupset.union(dupset[0],dup)

  genos = load_genostream(options.genotypes,format=options.informat,genorepr=options.ingenorepr,
                                   genome=options.loci,phenome=options.pedigree,
                                   transform=options)

  if options.testexp and options.testpairs:
    raise ValueError('Only one of --testexp and --testpairs may be specified')

  if options.testpairs:
    pairs,pair_count = file_pairs(genos, options.testpairs)

  elif options.testexp:
    pairs,pair_count = expected_pairs(genos, expected_dupset)
  else:
    genos = genos.as_sdat().materialize()
    pairs = pair_generator(genos)
    pair_count = len(genos)*(len(genos)-1)//2

  status = { (True, True ): ['EXPECTED',  'CONCORDANT'],
             (True, False): ['EXPECTED',  'DISCORDANT'],
             (False,True ): ['UNEXPECTED','CONCORDANT'] }

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['SAMPLE1','SAMPLE2','CONCORDANT_GENOTYPES','COMPARISONS','CONCORDANCE_RATE',
                'EXPECTED_DUPLICATE','OBSERVED_DUPLICATE'])

  if options.progress:
    pairs = progress_loop(pairs, length=pair_count, units='pairs')

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

  if options.dupout:
    write_dupsets(options.dupout, observed_dupset)


if __name__ == '__main__':
  main()
