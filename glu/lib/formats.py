# -*- coding: utf-8 -*-
'''
File:          formats.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import re

from fileutils import autofile

GENO_HEADER = 'rs#\tchr\tpos\t'
HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']
LOCUS_HEADER1 = ['LNAME','LOCATION','MAF','BINNUM','DISPOSITION']
LOCUS_HEADER2 = ['LNAME','LOCATION','POPULATION','MAF','BINNUM','DISPOSITION']

re_spaces = re.compile('[\t ,]+')


def hapmap_geno(g):
  return intern(g.strip().replace('N',' '))


def linkage_geno(g):
  return intern(''.join(g).strip().replace('0',' '))


def load_hapmap_genotypes(filename, nonfounders):
  '''
     Loads the RS#, Genome location, and genotypes from a HapMap formatted
     file.  Returns a generator that yields successive lists with the
     following elements:

       1) RS#          (string)
       2) Location     (integer)
       3..n) Genotypes (strings)

     Genotypes strings are 'interned' to save massive amounts of memory.
     i.e. all 'AA' string objects refer to the same immutable string object,
     rather than each 'AA' genotype allocating a distinct string object.
     See the Python builtin 'intern' function for more details.

     Remember: Generators are NOT generally restartable, so materialize the
               results if your code needs to iterate through them more than
               once.
  '''

  gfile = autofile(filename, 'r', hyphen=sys.stdin)

  gfile = dropwhile(lambda s: s.startswith('#'), gfile)
  header = gfile.next()

  if not any(header.startswith(h) for h in HAPMAP_HEADERS):
    if filename == '-':
      filename = '<stdin>'
    raise TagZillaError, "ERROR: HapMap Input file '%s' does not seem to be a HapMap data file." % filename

  header = header.strip().split(' ')

  if nonfounders is None:
    indices = range(11,len(header))
  else:
    indices = [ i for i,name in enumerate(header) if i >= 11 and name not in nonfounders ]

  for line in gfile:
    fields = line.split(' ')
    n = len(fields)
    genos = [ hapmap_geno(fields[i]) for i in indices if i < n ]
    try:
      yield Locus(fields[0], int(fields[3]), genos)
    except ValueError:
      # Ignore invalid loci, with just a warning, for now
      print >> sys.stderr, "WARNING: Invalid locus in file '%s', name '%s'" % (filename,fields[0])


def read_locus_file(filename):

  locfile = autofile(filename)
  locus_info = []

  try:
    for i,line in enumerate(locfile):
      fields = re_spaces.split(line.strip())

      if not fields:
        continue

      if len(fields) == 2:
        locus_info.append( (fields[0], int(fields[1])) )
      else:
        raise ValueError

  except (ValueError,IndexError):
    raise TagZillaError, 'ERROR: Invalid line %d in locus file "%s".' % (i+1,filename)

  return locus_info


def linkage_genotypes(fields):
  '''Return an iterator that yields genotype from a Linkage record'''
  allele1s = islice(fields,6,None,2)
  allele2s = islice(fields,7,None,2)
  return map(linkage_geno, izip(allele1s,allele2s))


def load_linkage_genotypes(filename, loci):
  '''
  Load each individual genotypes from the linkage formatted file
  for all the founders and assign the loaded genotypes to each locus
  in loci and construct a list of Locus objects.
  Each Locus object would have the following elements:
      1) RS#                 (string)
      2) Location            (integer)
      3) MAF                 (float)
      3) Genotypes           (list of strings)
  Note that Genotypes strings are 'interned' to save massive amounts of memory.
  i.e. all 'AA' string objects refer to the same immutable string object,
  rather than each 'AA' genotype allocating a distinct string object.
  See the Python builtin 'intern' function for more details.

  @param filename: name of the linkage formatted genotype data file
  @type  filename: string
  @param     loci: the locus information
  @type      loci: set of tuples with each like (locusname, location)
  @rtype         : list
  @return        : list of Locus objects
  '''

  pedfile = autofile(filename)

  ind_genos = []
  for line in pedfile:
    fields = re_spaces.split(line.strip())

    # Filter out all non-founders
    if len(fields) < 10 or fields[2:4] != ['0','0']:
      continue

    ind_genos.append( linkage_genotypes(fields) )

  n = len(ind_genos)
  missing_geno = intern('  ')
  loci = [ Locus(lname, location, [missing_geno]*n) for lname,location in loci ]

  for i,genos in enumerate(ind_genos):
    for j,g in enumerate(genos):
      loci[j].genos[i] = g

  for locus in loci:
    locus.maf = estimate_maf(locus.genos)

  return loci


def load_raw_genotypes(filename, nonfounders=None):
  '''
     Loads the RS#, Genome location, and genotypes from a native formatted
     genotype file.  Returns a generator that yields successive lists with
     the following elements:

       1) RS#          (string)
       2) Location     (integer)
       3..n) Genotypes (strings)

     Genotypes strings are 'interned' to save massive amounts of memory.
     i.e. all 'AA' string objects refer to the same immutable string object,
     rather than each 'AA' genotype allocating a distinct string object.
     See the Python builtin 'intern' function for more details.

     Remember: Generators are NOT generally restartable, so materialize the
               results if your code needs to iterate through them more than
               once.
  '''

  gfile = autofile(filename, 'r', hyphen=sys.stdin)

  header = gfile.readline()

  if not header.startswith(GENO_HEADER):
    if filename == '-':
      filename = '<stdin>'
    raise TagZillaError, "ERROR: Genotype input file '%s' does not seem to be in the right format." % filename

  header = header.strip().split('\t')

  if nonfounders is None:
    indices = range(3,len(header))
  else:
    indices = [ i for i,name in enumerate(header) if i >= 3 and name not in nonfounders ]

  for line in gfile:
    fields = line.split('\t')
    n = len(fields)
    genos = [ intern(fields[i].strip() or '  ') for i in indices if i < n ]
    yield Locus(fields[0], int(fields[2]), genos)


def load_festa_file(filename, locusmap, subset, rthreshold):
  '''
  Load FESTA formatted file that contain pre-computed LD data for pairs of loci
  '''
  ldfile = autofile(filename)
  header = ldfile.readline()

  for line in ldfile:
    lname1,lname2,ldvalue = re_spaces.split(line.strip())
    ldvalue = float(ldvalue)

    if subset and (lname1 not in subset or lname2 not in subset):
      continue

    if lname1 not in locusmap:
      locusmap[lname1] = Locus(lname1, 0, [])

    if lname2 not in locusmap:
      locusmap[lname2] = Locus(lname2, 0, [])

    if ldvalue >= rthreshold:
      yield lname1,lname2,ldvalue,0


def load_hapmapld_file(filename, locusmap, subset, maxd, rthreshold, dthreshold):
  '''
  Load Hapmap formatted file that contain pre-computed LD data for pairs of loci
  '''
  ldfile = autofile(filename)
  ldfile = dropwhile(lambda s: s.startswith('#'), ldfile)

  for line in ldfile:
    loc1,loc2,pop,lname1,lname2,dprime,r2,lod = line.strip().split(' ')

    if subset and (lname1 not in subset or lname2 not in subset):
      continue

    loc1 = int(loc1)
    loc2 = int(loc2)

    if lname1 not in locusmap:
      locusmap[lname1] = Locus(lname1, loc1, [])

    if lname2 not in locusmap:
      locusmap[lname2] = Locus(lname2, loc2, [])

    if abs(loc1-loc2) > maxd:
      continue

    dprime = float(dprime)
    r2     = float(r2)

    if r2 >= rthreshold and abs(dprime) >= dthreshold:
      yield lname1,lname2,r2,dprime


def read_hapmap_nonfounders(filename, nonfounders):
  pedfile = autofile(filename)

  founder = ['0','0']
  for line in pedfile:
    fields = line.strip().split('\t')

    # Filter out all founders
    if fields[2:4] == founder or len(fields) < 7:
      continue

    sample = fields[6].split(':')[-2]
    nonfounders.add(sample)

  return nonfounders


def read_linkage_nonfounders(filename, nonfounders):
  pedfile = autofile(filename)

  founder = ['0','0']
  for line in pedfile:
    fields = re_spaces.split(line.strip())

    # Filter out all non-founders
    if fields[2:4] == founder:
      continue

    nonfounders.add(tuple(fields[0:2]))

  return nonfounders


def read_snp_list(name, sset):
  if name.startswith(':'):
    sset.update(name[1:].split(','))
    return

  sfile = autofile(name)
  for line in sfile:
    fields = re_spaces.split(line.strip())
    if fields:
      sset.add(fields[0])

  return sset


def read_design_score(filename):
  sf = autofile(filename)
  for line in sf:
    fields = re_spaces.split(line.strip())
    lname = fields[0]
    try:
      score = float(fields[1])
      yield lname,score

    except ValueError:
      pass


def read_illumina_design_score(filename):
  sf = autofile(filename)
  header = sf.next.split(',')
  design_index = header.index('SNP_Score')
  for line in sf:
    fields = line.split(',')
    lname = fields[0]
    try:
      score = float(fields[design_index])
      yield lname,score

    except ValueError:
      pass

