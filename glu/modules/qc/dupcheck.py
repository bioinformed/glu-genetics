# -*- coding: utf-8 -*-
'''
File:          dupcheck.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   itertools          import izip, chain
from   textwrap           import fill

from   glu.lib.utils      import percent
from   glu.lib.fileutils  import autofile, hyphen, load_map
from   glu.lib.union_find import union_find
from   glu.lib.genolib    import load_genomatrixstream, snp
from   glu.lib.sections   import save_section, SectionWriter, save_metadata_section


try:
  from dupcheckc import match_byte_array as match

except ImportError:
  print >> sys.stderr, 'WARNING: Using pure-Python version of dupcheck match function'

  def match(g1, g2):
    i = n = 0
    for a,b in izip(g1,g2):
      # If data are not missing
      if a and b:
        # Check for match
        if a is b:
          i += 1
        # Count total comparisons
        n += 1

    return i,n


def duplicate(genos1, genos2, threshold=0.75):
  matches,genos = match(genos1, genos2)

  dup = genos and float(matches)/genos >= threshold
  return dup,matches,genos


def generate_pairs(data):
  n = len(data)

  for i in xrange(n):
    i1 = data[i]

    for j in xrange(i+1,n):
      i2 = data[j]

      yield i1,i2


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
    text = fill(', '.join(s), initial_indent=e,subsequent_indent=e,width=79).lstrip(' ')
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
  out = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
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
  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='A csv file containing expected duplicates')
  parser.add_option('-d', '--dupout', dest='dupout', metavar='FILE',
                    help='Output of duplicate sets')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity chasred between two individuals (default=85)')
  parser.add_option('-m', '--mincount', dest='mincount', metavar='N', type='int', default=20,
                    help='Minimum concordant genotypes to be considered informative for duplicate checking')
  parser.add_option('-s', '--samplemap', dest='samplemap', metavar='FILE',
                    help='Remap sample names based on map file')
  parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')

  return parser


def main():
  import sys,time

  parser = option_parser()
  options,args = parser.parse_args()

  threshold = float(options.threshold)/100

  if len(args) != 1:
    parser.print_help()
    return

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  if options.samplemap:
    samplemap = load_map(options.samplemap)

  expected_dupset = union_find()
  observed_dupset = union_find()

  if options.duplicates:
    print >> sys.stderr, 'Loading expected duplicate data...',
    for dupset in csv.reader(autofile(options.duplicates), 'excel-tab'):
      if options.samplemap:
        dupset = [ samplemap[dup] for dup in dupset ]
      for dup in dupset:
        expected_dupset.union(dupset[0],dup)
    print >> sys.stderr, 'Done.'

  print >> sys.stderr, 'Loading data...',
  samples = load_genomatrixstream(args[0], options.format, snp).as_sdat()

  if not samples.packed:
    samples = samples.transformed(repack=True)

  if options.samplemap:
    samples = samples.transformed(rename_samples=options.samplemap)

  # Materialize samples
  samples = list(samples)

  print >> sys.stderr, 'Done.'

  print >> sys.stderr, 'Checking for duplicates...',

  out.write('Duplicate sample analysis\n')
  out.write('  Input File:     %s\n' % args[0])
  out.write('  Timestamp:      %s\n' %  time.asctime())
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

''' % (threshold*100,options.mincount))

  pairs = generate_pairs(samples)

  expected_dups      = []
  unexpected_dups    = []
  unexpected_nondups = []

  unexpected_dupset = union_find()

  for (i1,genos1),(i2,genos2) in pairs:
    obs_dup,matches,genos = duplicate(genos1, genos2, threshold)

    if matches < options.mincount:
      obs_dup = False

    exp_dup = (i1,i2) in expected_dupset

    if not obs_dup and not exp_dup:
      continue

    if obs_dup:
      observed_dupset.union(i1,i2)

    d = (i1,i2,matches,genos)
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

  n,m = len(expected_dups)+len(unexpected_dups),len(samples)
  out.write('All DUPLICATE INDIVIDUALS: %d/%d = %5.2f%%\n' % (n,m,percent(n,m)))
  alldups = chain(expected_dups,unexpected_dups)
  duplicate_output_summary(out,sum_dupsets(alldups, observed_dupset))

  n,m = len(expected_dups),len(samples)
  out.write('EXPECTED DUPLICATE INDIVIDUALS: %d/%d = %5.2f%%\n' % (n,m,percent(n,m)))
  duplicate_output_summary(out,sum_dupsets(expected_dups))

  n,m = len(unexpected_dups),len(samples)
  out.write('UNEXPECTED DUPLICATE INDIVIDUALS: %d/%d = %5.2f%%\n' % (n,m,percent(n,m)))
  duplicate_output_summary(out,sum_dupsets(unexpected_dups))

  n,m = len(unexpected_nondups),len(samples)
  out.write('EXPECTED DUPLICATE INDIVIDUALS NOT DETECTED: %d/%d = %5.2f%%\n' % (n,m,percent(n,m)))
  duplicate_output_summary(out,sum_dupsets(unexpected_nondups))

  out.write('EXPECTED DUPLICATE INDIVIDUALS DETECTED\n')
  duplicate_output_details(out,expected_dups)

  out.write('UNEXPECTED DUPLICATE INDIVIDUALS DETECTED\n')
  duplicate_output_details(out,unexpected_dups)

  out.write('EXPECTED DUPLICATE INDIVIDUALS NOT DETECTED\n')
  duplicate_output_details(out,unexpected_nondups)


if __name__ == '__main__':
  main()
