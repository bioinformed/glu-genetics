# -*- coding: utf-8 -*-

from   __future__ import division

__gluindex__  = False
__abstract__  = '''\
Detect expected and unexpected duplicate samples based on expected sample
equivalence and empirical genotype concordance rate'''
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import time
import csv

from   itertools                 import chain
from   textwrap                  import fill

from   glu.lib.utils             import percent, pair_generator
from   glu.lib.fileutils         import autofile, hyphen, table_writer
from   glu.lib.union_find        import union_find
from   glu.lib.sections          import save_section, SectionWriter, save_metadata_section

from   glu.lib.genolib.io        import load_genostream
from   glu.lib.genolib.merge     import get_genomerger
from   glu.lib.genolib.genoarray import genoarray_concordance


# Here for comparison purposes to experiment with unoptomized genotype
# comparisons
def block_pair_generator(samples,blocksize=100):
  '''
  Generate all distinct pairs of genotype vectors using a blocking
  algorithm.  Blocking is helpful because genotypes are stored in a
  bit-packed format and unpacking is fairly expensive.  Thus they are
  unpacked in chunks so that the total number of unpacking operations is
  reduced by a factor of the blocksize.

  Example using 6 input samples and a block size of 2:

  >>> samples = [ (i,'ACGT') for i in range(5) ]
  >>> pairs = list(block_pair_generator(samples,blocksize=2))

  We expect (n-1)(n-2)/2 pairs, a total of 10 for 6 inputs

  >>> len(pairs)
  10

  We can see that the pairs are distinct and seem to cover all of the
  possibilities.

  >>> for p1,p2 in pairs:
  ...   print p1,p2
  (0, 'ACGT') (1, 'ACGT')
  (0, 'ACGT') (2, 'ACGT')
  (0, 'ACGT') (3, 'ACGT')
  (1, 'ACGT') (2, 'ACGT')
  (1, 'ACGT') (3, 'ACGT')
  (0, 'ACGT') (4, 'ACGT')
  (1, 'ACGT') (4, 'ACGT')
  (2, 'ACGT') (3, 'ACGT')
  (2, 'ACGT') (4, 'ACGT')
  (3, 'ACGT') (4, 'ACGT')

  This is more easily seen after sorting the pairs.

  >>> for p1,p2 in sorted(pairs):
  ...   print p1,p2
  (0, 'ACGT') (1, 'ACGT')
  (0, 'ACGT') (2, 'ACGT')
  (0, 'ACGT') (3, 'ACGT')
  (0, 'ACGT') (4, 'ACGT')
  (1, 'ACGT') (2, 'ACGT')
  (1, 'ACGT') (3, 'ACGT')
  (1, 'ACGT') (4, 'ACGT')
  (2, 'ACGT') (3, 'ACGT')
  (2, 'ACGT') (4, 'ACGT')
  (3, 'ACGT') (4, 'ACGT')
  '''
  blocks = (len(samples)+blocksize-1)//blocksize

  for i in xrange(blocks):
    block1 = [ (id,row[:]) for id,row in samples[i*blocksize:(i+1)*blocksize] ]
    n = len(block1)

    # Output within-block pairs
    for j in range(n):
      for k in range(j+1,n):
        yield block1[j],block1[k]

    # Output inter-block pairs
    for j in xrange(i+1,blocks):
      block2 = [ (id,row[:]) for id,row in samples[j*blocksize:(j+1)*blocksize] ]
      m = len(block2)

      for j in range(n):
        for k in range(m):
          yield block1[j],block2[k]


def sum_dupsets(data, union=None):
  dups = {}

  if union is None:
    if not isinstance(data, (list,tuple)):
      data = list(data)

    union = union_find()
    for a,b,c,n in data:
      union.union(a,b)

  for a,b,c,n in data:
    exemplar=union.find(a)
    assert exemplar == union.find(b)
    dupdata = dups.setdefault(exemplar, [0,0,0])
    dupdata[0] += c
    dupdata[1] += n
    dupdata[2] += 1.0

  sets = union.setmap()
  return [ (sets[e],c/m,n/m) for e,(c,n,m) in dups.iteritems() ]


def duplicate_output_summary(out,data):
  data = [ (percent(c,n),len(s),s) for s,c,n in data ]
  data.sort(reverse=True)
  out.write('\n')
  out.write('  Rank  %Conc   Size  Members\n')
  out.write('  ----  ------  ----  --------------------------------------------------------\n')

  e = ' '*22
  for r,(p,m,s) in enumerate(data):
    text = fill(', '.join(sorted(s)), initial_indent=e,subsequent_indent=e,width=79).lstrip(' ')
    out.write('  %4d  %5.1f%%  %4d  %s\n' % (r+1,p,m,text))
  out.write('\n\n')


def duplicate_output_details(out,data):
  data = sorted( ((percent(c,n),c,n,a,b) for a,b,c,n in data), reverse=True)
  out.write('\n')
  out.write('  Rank  Individual1                 Individual2                  N / Total   %Conc  \n')
  out.write('  ----  --------------------------  --------------------------  -----------  -------\n')

  for r,(p,c,n,a,b) in enumerate(data):
    out.write('  %4d  %-27s %-27s %4d / %4d  %6.2f%%\n' % (r+1,a,b,c,n,p))

  out.write('\n\n')


def write_dupsets(filename, dupsets):
  out = table_writer(filename)
  for dset in dupsets.sets():
    out.writerow( list(dset) )


def save_results(sw, observed_dupset, expected_dups, unexpected_dups, unexpected_nondups):
  save_section(sw, 'expected_duplicates',             map(list, observed_dupset.sets()) )
  save_section(sw, 'expected_duplicates_detected',    [(i1,i2,m,g) for i1,i2,m,g in expected_dups] )
  save_section(sw, 'unexpected_duplicates_detected',  [(i1,i2,m,g) for i1,i2,m,g in unexpected_dups] )
  save_section(sw, 'expected_duplicates_undetected',  [(i1,i2,m,g) for i1,i2,m,g in unexpected_nondups] )


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f','--format',  dest='format',
                    help='Input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('--merge', dest='merge', metavar='METHOD:T', default='unanimous',
                    help='Genotype merge algorithm and optional consensus threshold used to form a consensus genotypes. '
                         'Values=unique,unanimous,vote,ordered.  Value may be optionally followed by a colon and a threshold.  Default=unanimous')
  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='A delimited file containing expected duplicates')
  parser.add_option('-d', '--dupout', dest='dupout', metavar='FILE',
                    help='Output of duplicate sets')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity shared between two individuals (default=85)')
  parser.add_option('-m', '--mingenos', '--mincount', dest='mingenos', metavar='N', type='int', default=20,
                    help='Minimum concordant genotypes to be considered informative for duplicate checking')
  parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  threshold = float(options.threshold)/100

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  expected_dupset = union_find()
  observed_dupset = union_find()

  if options.duplicates:
    print >> sys.stderr, 'Loading expected duplicate data...',
    for dupset in csv.reader(autofile(options.duplicates), dialect='excel-tab'):
      for dup in dupset:
        expected_dupset.union(dupset[0],dup)
    print >> sys.stderr, 'Done.'

  print >> sys.stderr, 'Loading data...',
  merger   = get_genomerger(options.merge)

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,genome=options.loci)
  genos = genos.as_sdat(merger).materialize()

  print >> sys.stderr, 'Done.'

  print >> sys.stderr, 'Checking for duplicates...',

  out.write('Duplicate sample analysis\n')
  out.write('  Input File:     %s\n' % args[0])
  out.write('  Timestamp:      %s\n' % time.asctime())
  out.write('''
  This analysis will detect putative duplicate samples within the given
  dataset that share at least %.2f%% of their genotypes in common, provided
  that there are at least %d comparable genotypes (the number of genotypes
  that can be compared will differ because of missing data).

Input:
  A CSV file with a row for each locus and a column for each individual's
  genotype.  The first row must contain the names of each column.  The
  locus name must be in the first column.

Output:
  1) A summary table of sets of duplicate samples.  These sets are
     constructed such that at least one member was selected as a
     putative duplicate of another member.  Thus, a low threshold
     can result in very large duplicate sets due to population
     structure or related individuals in the dataset.  The number of
     individuals in each set as well as the average pair-wise genotype
     concordance rate is displayed for each set.

  2) A detailed output table that lists the number of concordant
     genotypes, number of comparable genotypes, and the concordance
     rates for each putative duplicate pair of individuals.

''' % (threshold*100,options.mingenos))

  expected_dups      = []
  unexpected_dups    = []
  unexpected_nondups = []

  unexpected_dupset = union_find()

  for (i1,genos1),(i2,genos2) in pair_generator(genos.use_stream()):
    matches,comparisons = genoarray_concordance(genos1, genos2)

    obs_dup = matches>=options.mingenos and comparisons and float(matches)/comparisons >= threshold
    exp_dup = (i1,i2) in expected_dupset

    if not obs_dup and not exp_dup:
      continue

    if i1>i2:
      i1,i2=i2,i1

    if obs_dup:
      observed_dupset.union(i1,i2)

    d = (i1,i2,matches,comparisons)
    if obs_dup and exp_dup:
      expected_dups.append(d)
      expected_dupset.union(i1,i2)
    elif obs_dup:
      unexpected_dups.append(d)
      unexpected_dupset.union(i1,i2)
    else:
      unexpected_nondups.append(d)

  print >> sys.stderr, 'Done.'

  if options.tabularoutput:
    sw = SectionWriter(options.tabularoutput)
    save_metadata_section(sw, analysis='dupcheck', analysis_version='0.1')
    save_results(sw, observed_dupset, expected_dups, unexpected_dups, unexpected_nondups)

  if options.dupout:
    union_out = autofile(options.dupout, 'w')
    write_dupsets(options.dupout, observed_dupset)

  out.write('All DUPLICATE INDIVIDUALS:\n')
  alldups = chain(expected_dups,unexpected_dups)
  duplicate_output_summary(out,sum_dupsets(alldups, observed_dupset))

  out.write('EXPECTED DUPLICATE INDIVIDUALS:\n')
  duplicate_output_summary(out,sum_dupsets(expected_dups))

  out.write('UNEXPECTED DUPLICATE INDIVIDUALS:\n')
  duplicate_output_summary(out,sum_dupsets(unexpected_dups))

  out.write('EXPECTED DUPLICATE INDIVIDUALS NOT DETECTED:\n')
  duplicate_output_summary(out,sum_dupsets(unexpected_nondups))

  out.write('EXPECTED DUPLICATE INDIVIDUALS DETECTED\n')
  duplicate_output_details(out,expected_dups)

  out.write('UNEXPECTED DUPLICATE INDIVIDUALS DETECTED\n')
  duplicate_output_details(out,unexpected_dups)

  out.write('EXPECTED DUPLICATE INDIVIDUALS NOT DETECTED\n')
  duplicate_output_details(out,unexpected_nondups)


def test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  test()
  main()
