# -*- coding: utf-8 -*-
'''
File:          dupcheck2.py

Authors:       Jun Lu          (lujun@mail.nih.gov)
               Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-10-01

Abstract:      Detect expected and unexpected duplicate samples by genotype
               concordance

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv
import time

from   itertools             import izip, groupby
from   operator              import itemgetter

from   glu.lib.utils         import autofile, hyphen, percent, pair_generator, pick
from   glu.lib.union_find    import union_find
from   glu.lib.genoarray     import get_genorepr
from   glu.lib.genomerge     import get_genomerger
from   glu.lib.genodata      import load_map, load_genostream
from   glu.lib.sections      import save_section, SectionWriter, save_metadata_section


def match_object_generic(g1, g2):
  i = n = 0
  for a,b in izip(g1,g2):
    # If data are not missing
    if a and b:
      # Check for match
      if a == b:
        i += 1
      # Count total comparisons
      n += 1

  return i,n


def get_matchfunc(exemplar):
  import array

  try:
    import dupcheckc

    if isinstance(exemplar, list):
      return dupcheckc.match_object

    elif isinstance(exemplar, array.array):
      if exemplar.typecode == 'B':
        return dupcheckc.match_byte_array
      elif exemplar.typecode == 'H':
        return dupcheckc.match_short_array
      elif exemplar.typecode == 'L':
        return dupcheckc.match_long_array
  except ImportError:
    pass

  print >> sys.stderr, 'WARNING: Using pure-Python/slow version of dupcheck match function'
  return match_object_generic


def write_header(outfile,filename,threshold,mincount):
  outfile.write('Duplicate sample analysis\n')
  outfile.write('  Input File:     %s\n' % filename)
  outfile.write('  Timestamp:      %s\n' %  time.asctime())
  outfile.write('  Threshold: %.2f%%; Mininum Counts: %d (see explanation below)\n' % (threshold,mincount))
  outfile.write('''
  This analysis will detect putative duplicate samples within the given
  dataset that share at least %.2f%% of their genotypes in common, provided
  that there are at least %d matched genotypes.

Output:
  This analysis provides the following output:
    - A side-by-side comparison between expected and observed duplicate sets.
      If any members of an expected duplicate set are not detected, or are found in
      2 or more observed sets, that set will be label with '*' (next to the set_size_number).
      If additional members appear in the same (observed) set with an expected duplicate set,
      an '!' will be added in the first column.
    - If phenotype data is available, this analysis will print out the detected
      duplicate sample-pairs with DIFFERENT phenotypes. These are problematic samples.
    - Concordance results for sample-pairs of three groups:
       1. unexpected
       2. expected, not detectd
       3. expected, detected
''' % (threshold,mincount))


def write_dupsets(filename, dupsets):
  out = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
  for dset in dupsets.sets():
    out.writerow( list(dset) )


def load_sample_phenos(file,phenos):
  phenos = [it.lower() for it in phenos.split(',')]

  rows = csv.reader(autofile(file), 'excel-tab')
  header = rows.next()
  pos = [i for (i,it) in enumerate(header) if it.lower() in phenos]

  if len(pos) != len(phenos):
    raise ValueError,'One or more sample phenotypes not found in %s, please check.' % file

  return dict([('phenos',tuple(pick(header,pos)))]+[(row[0],tuple(pick(row,pos))) for row in rows])


def load_expected_dupsets(file):
  expected_dupset = union_find()

  for dupset in csv.reader(autofile(file), 'excel-tab'):
    dupset = [ d for d in dupset if d ]
    for dup in dupset:
      expected_dupset.union(dupset[0],dup)

  return expected_dupset


def find_dups(samples,threshold=85,mincount=20,exp_dups=None,sample_phenos=None):
  threshold = float(threshold)/100
  samples = iter(samples)
  samples.next()
  samples = list(samples)
  sample_ids = set(map(itemgetter(0),samples))

  dup_pairs = []
  dup_sets  = union_find()

  dup_pairs.append( (None,None,None,None,None,1,1) )
  dup_pairs.append( (None,None,None,None,None,0,1) )
  dup_pairs.append( (None,None,None,None,None,1,0) )

  matchfunc = get_matchfunc(samples[0][1])

  if exp_dups is None:
    exp_dups = set()

  for (s1,genos1),(s2,genos2) in pair_generator(samples):
    matches,total = matchfunc(genos1, genos2)
    if len(genos1)< mincount or len(genos2) < mincount:
      raise ValueError, 'Error: the minimum counts (-m option) can not be greater that the number of loci'

    obs_dup = total and matches >= mincount and float(matches)/total >= threshold

    s1,s2 = sorted([s1,s2])
    is_exp= (s1,s2) in exp_dups

    if not obs_dup and not is_exp:
      continue

    diff_pheno = 'NO|UNKNOWN'
    if sample_phenos is not None:
      s1_pheno,s2_pheno = sample_phenos.get(s1),sample_phenos.get(s2)
      if not s1_pheno or not s2_pheno:
        diff_pheno = 'UNKNOWN'
      elif s1_pheno != s2_pheno:
        diff_pheno = 'YES'
      else:
        diff_pheno = 'NO'

    if obs_dup:
      dup_sets.union(s1,s2)
      dup_pairs.append( (s1,s2,matches,total,diff_pheno,1, is_exp) )
    else:
      dup_pairs.append((s1,s2,matches,total,diff_pheno,0,1))

  return sample_ids,dup_pairs,dup_sets


def compare_dupsets(sample_ids,obs_dupsets,exp_dupsets=None):
  obs_dups   = sorted([sorted(list(i)) for i in obs_dupsets.sets()])
  obs_setmap = {}
  for i,alist in enumerate(obs_dups):
    for s in alist:
      obs_setmap[s] = i+1

  exp      = []
  unexp    = []
  leftover = []
  obs_idx  = set()

  if exp_dupsets is not None:
    exp_dups = [sorted([j for j in s if j in sample_ids]) for s in exp_dupsets.sets()]
    exp_dups = [it for it in exp_dups if len(it)>1]
    exp_dups.sort()
    for i,alist in enumerate(exp_dups):
      exp.append([(s,obs_setmap.get(s,'NOT_FOUND')) for s in alist])
      unexp_sets = set()
      for s in alist:
        if s in obs_setmap:
          unexp_sets |= set(obs_dups[obs_setmap[s]-1])
          obs_idx.add(obs_setmap[s]-1)

      unexp_sets -= set(alist)
      unexp.append(sorted([(s,obs_setmap[s]) for s in unexp_sets]))

  leftover = [(sorted(list(aset)),i+1) for i,aset in enumerate(obs_dups) if i not in obs_idx]

  return exp,unexp,leftover


def write_title_rows(out,h1,h2,sample_phenos=None):
  if sample_phenos is not None:
    for pheno in sample_phenos['phenos']:
      h1 += ' %-20s ' % pheno
      h2 += ' '+'-'*20 + ' '
  out.write(h1 + '\n' + h2 + '\n')


def write_sample_rows(out,samples,first_ind,next_ind,in_obs_set=None,sample_phenos=None):
  if in_obs_set is not None:
    assert len(samples) == len(in_obs_set)

  for i, s in enumerate(samples):
    indent = (first_ind if not i else next_ind)
    row = indent + '%-27s' % s

    if in_obs_set is not None:
      row += '  %15s  ' % str(in_obs_set[i])
    if sample_phenos is not None:
      row += ''.join([' %-20s ' % p for p in sample_phenos.get(s,'')])
    out.write(row+'\n')

  out.write('\n')


def write_unexpected_sets(out,leftover,sample_phenos=None):
  e = 11
  indent2 = ' ' * e

  out.write('    *****  ' + '[OTHER OBSERVED DUPLICATE SETS (i.e. UNEXPECTED)]' + '\n')
  for (alist,i) in leftover:
    indent1 = '  ' + '%7d'% (len(alist)) + '  '
    write_sample_rows(out,alist,indent1,indent2,[i]*len(alist),sample_phenos)


def write_pairs_long(out,indent,samples,sample_phenos=None):
  s1,s2 = samples
  row = indent + '%-27s  %-27s ' % (s1,s2)
  if sample_phenos is not None:
    num_pheno = len(sample_phenos['phenos'])
    pheno1 = sample_phenos.get(s1,['']*num_pheno)
    pheno2 = sample_phenos.get(s2,['']*num_pheno)
    row += ''.join([' %-20s ' % ','.join([p1,p2]) for p1,p2 in izip(pheno1,pheno2)])

  out.write(row+'\n')


def write_dups_diff_phenos(out,dups,sample_phenos=None):
  diff_dups = [(percent(m,n),s1,s2,m,n,exp) for s1,s2,m,n,d,obs,exp in dups if d=='YES']
  diff_dups.sort(reverse=True)

  out.write('\n\n <!>-- Duplicate sample-pairs with DIFFERENT phenotypes --<!>\n\n')
  if sample_phenos is None:
    out.write('   <NOT AVAILABLE, PHENOTYPE INFORMATION NEEDED.> \n')
  else:
    header1 = '  %Concord       N / Total        In_Expected_Set  Duplicate-Pair             '
    header2 = '  ---------  -------------------  ---------------  -------------------------- '
    write_title_rows(out,header1,header2,sample_phenos)
    ind2 = ' '*51
    for p,s1,s2,m,n,e in diff_dups:
      expected = ('YES' if e else 'NO')
      ind1     = '  %8.2f%%  %8d / %8d  %15s  ' % (p,m,n,expected)
      write_sample_rows(out,[s1,s2],ind1,ind2,sample_phenos=sample_phenos)


def write_duplicate_sets(out,sample_ids,obs_dupsets,exp_dupsets=None,sample_phenos=None):
  exp,unexp,leftover = compare_dupsets(sample_ids,obs_dupsets,exp_dupsets)

  out.write('\n\n Comparison of EXPECTED with OBSERVED duplicate-sets\n\n')
  h1 = '  Set_size Expected_duplicate_set      In_observed_set   '
  h2 = '  -------- --------------------------  ----------------  '
  write_title_rows(out,h1,h2,sample_phenos)

  ind2 = ' ' * 11
  if not exp:
    write_unexpected_sets(out,leftover,sample_phenos)
    return

  for i,aset in enumerate(exp):
    samples = [s for s,pos in aset]
    in_obs_set = [pos for s,pos in aset]
    ind1    = '  ' + '%7d'% (len(aset)) + '  '
    if 'NOT_FOUND' in in_obs_set or len(set(in_obs_set))!=1:
      ind1    = '  ' + '%6d*'% (len(aset)) + '  '
    write_sample_rows(out,samples,ind1,ind2,in_obs_set,sample_phenos)
    if not unexp[i]:
      continue

    out.write('       !   '+ '[OTHER MEMBERS IN THE SAME OBSERVED SET]' + '\n')
    samples = [s for s,pos in unexp[i]]
    in_obs_set = [pos for s,pos in unexp[i]]
    write_sample_rows(out,samples,ind2,ind2,in_obs_set,sample_phenos)

  if leftover:
    write_unexpected_sets(out,leftover,sample_phenos)


def write_duplicate_pairs(out,dups,sample_phenos=None,pairformat='long'):
  dups = sorted(dups,key=itemgetter(5,6))
  for key,group in groupby(dups,itemgetter(5,6)):
    group = sorted( ((percent(m,n),s1,s2,m,n) for s1,s2,m,n,d,obs,expect in group),reverse=True )
    if key == (1,0):
      g1 = (group,' I. Unexpected duplicate-pairs\n\n')
    elif key == (0,1):
      g2 = (group,' II. Expected duplicate-pairs, but NOT detected\n\n')
    else:
      g3 = (group,' III. Expected duplicate-pairs, and detected\n\n')

  for group,title in (g1,g2,g3):
    out.write('\n\n%s ' % title)
    if pairformat == 'short':
      h1 = '  %Concord       N / Total        Duplicate-Pair            '
      h2 = '  ---------  -------------------  ------------------------- '
      write_title_rows(out,h1,h2,sample_phenos)
      ind2  = ' '*34
      for p,s1,s2,m,n in group:
        if m is None: continue
        ind1 = '  %8.2f%%  %8d / %8d  ' % (p,m,n)
        write_sample_rows(out,[s1,s2],ind1,ind2,sample_phenos=sample_phenos)
    elif pairformat == 'long':
      h1 = '  %Concord       N / Total        Sample_1                     Sample_2                    '
      h2 = '  ---------  -------------------  ---------------------------  --------------------------- '
      write_title_rows(out,h1,h2,sample_phenos)
      for p,s1,s2,m,n in group:
        if m is None: continue
        ind = '  %8.2f%%  %8d / %8d  ' % (p,m,n)
        write_pairs_long(out,ind,(s1,s2),sample_phenos)
    else:
      raise NotImplementedError,'The pair-format: %s is not supported!' % pairformat


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-e', '--duplicates', dest='duplicates', metavar='FILE',
                    help='A tab-delimited file with expected duplicates on each row, or with columns of sid,pid')
  parser.add_option('-d', '--dupout', dest='dupout', metavar='FILE',
                    help='Output of duplicate sets')
  parser.add_option('-f','--format',  dest='format', metavar='string', default='sdat',
                    help='The file input format for genotype data. Values=hapmap, ldat, sdat (default), trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp_acgt',
                    help='Input genotype representation.  Values=snp_acgt (default), snp_ab, snp_marker, or generic')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of duplicate check report')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity shared between two individuals (default=85)')
  parser.add_option('-m', '--mincount', dest='mincount', metavar='N', type='int', default=20,
                    help='Minimum concordant genotypes to be considered informative for duplicate checking')
  parser.add_option('--merge', dest='merge', metavar='METHOD:T', default='vote:1',
                    help='Genotype merge algorithm and optional consensus threshold used to form a consensus genotypes. '
                         'Values=vote,ordered.  Value may be optionally followed by a colon and a threshold.  Default=vote:1')
  parser.add_option('-l', '--limit', dest='limit', metavar='N', type='int', default=None,
                    help='Limit the number of rows of data to N for testing purposes')
  #parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE',
  #                  help='Generate machine readable tabular output of results')

  parser.add_option('--phenofile', dest='phenofile', metavar='FILE',default=None,
                    help='A file containing sample phenotype information')
  parser.add_option('--phenos', dest='phenos', metavar='string',type='string',default=None,
                    help='Specify sample phenotypes (i.e. column headers in the phenotype file) to be included in output')
  parser.add_option('--pairformat', dest='pairformat', metavar='string',type='string',default='long',
                    help='Specify the output format (options: long or short)')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  print >> sys.stderr, 'Loading genotype data...'
  genorepr  = get_genorepr(options.genorepr)
  limit     = options.limit or None
  merger    = get_genomerger(options.merge,genorepr)
  samples   = load_genostream(args[0], options.format, limit=limit, genorepr=genorepr, unique=False).as_sdat(merger)

  exp_dupsets = None
  if options.duplicates:
    print >> sys.stderr, 'Loading expected duplicates data...'
    exp_dupsets = load_expected_dupsets(options.duplicates)

  sample_phenos = None
  if options.phenos:
    if options.phenofile is not None:
      sample_phenos = load_sample_phenos(options.phenofile,options.phenos)
    else:
      print >> sys.stderr, 'Sample phenotype file is required'
      return

  print >> sys.stderr, 'Checking for duplicates...'
  sample_ids,dups,dup_sets = find_dups(samples,options.threshold,options.mincount,exp_dupsets,sample_phenos)

  if options.dupout:
    union_out = autofile(options.dupout, 'w')
    write_dupsets(options.dupout, dup_sets)

  out = autofile(hyphen(options.output,sys.stdout), 'w')

  write_header(out,args[0],options.threshold,options.mincount)
  write_duplicate_sets(out,sample_ids,dup_sets,exp_dupsets,sample_phenos)
  write_dups_diff_phenos(out,dups,sample_phenos)
  write_duplicate_pairs(out,dups,sample_phenos,options.pairformat)


if __name__ == '__main__':
  main()
