# -*- coding: utf-8 -*-

__program__   = 'TagZilla'
__gluindex__  = True
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)',
                 'Zhaoming Wang (wangzha@mail.nih.gov)']
__abstract__  = 'Robust and fast SNP tagging program'
__descr__     = '''\
A robust and fast SNP binning and tagging program that takes many forms of
input data, can be tuned by over a dozen meaningful parameters, and produces
several useful human and machine readable outputs.  The heart of the program
takes genotype data, haplotype frequencies, from which pairwise r-squared or
D' linkage disequilibrium statistics are computed.  Those LD statistics are
used to compute bins using a greedy maximal algorithm (similar to that of
Carlson et al, 2004), and reports detailed information on bins and tags.
Many useful extensions are also implemented, including sex-linked analysis,
efficient multi-population tagging, incorporation of design scores and bin
informativity measures, and calculation of detailed bin and locus coverage
statistics by type of bin.  Please consult the accompanying manual for more
information.'''
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import re
import sys
import copy
import time

from   math                      import log, ceil
from   operator                  import itemgetter
from   collections               import defaultdict
from   itertools                 import chain, groupby, izip, dropwhile, count

from   glu.lib.hwp               import hwp_biallelic
from   glu.lib.stats             import mean, median
from   glu.lib.utils             import pair_generator, percent
from   glu.lib.fileutils         import autofile, hyphen, list_reader, table_reader, table_writer
from   glu.lib.glu_launcher      import GLUError
from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.ld        import estimate_ld, count_haplotypes, bound_ld
from   glu.lib.genolib.genoarray import minor_allele_from_genos


epsilon = 10e-10

LOCUS_HEADER1   = ['LNAME','LOCATION','MAF','BINNUM','DISPOSITION']
LOCUS_HEADER2   = ['LNAME','LOCATION','POPULATION','MAF','BINNUM','DISPOSITION']
LOCUS_HEADER3   = ['LNAME','CHROMOSOME','LOCATION','POPULATION','MAF','BINNUM','DISPOSITION']
PAIR_HEADER     = ['BIN','LNAME1','LNAME2','POPULATION','RSQUARED','DPRIME','DISPOSITION']
re_spaces = re.compile('[\t ,]+')


class TagZillaError(GLUError): pass


class Locus(object):
  __slots__ = ('name','chromosome','location','maf','genos')
  def __init__(self, name, chromosome, location, genos):
    a,maf = minor_allele_from_genos(genos)

    self.name       = name
    self.chromosome = chromosome
    self.location   = location or 0
    self.maf        = maf
    self.genos      = genos


def scan_ldpairs(loci, maxd, rthreshold, dthreshold):
  '''
  A generator for pairs of loci within a specified genomic distance.
  Loci are assumed to be sorted by genomic location.

  Split loci into non-communicating regions based on chromosome boundaries
  or gaps larger than maxd and generates pairwise ld using
  scan_ldpairs_region
  '''
  regions  = []
  region   = []
  last_chr = None
  last_loc = 0

  # Collect all regions
  for locus in loci:
    chr = locus.chromosome
    loc = locus.location

    # Detect boundary and split the current region
    if (chr!=last_chr or loc-last_loc>maxd) and region:
      regions.append(region)
      region = []

    last_chr = chr
    last_loc = loc
    region.append(locus)

  # Pick up tailings
  if region:
    regions.append(region)

  # Generate ld and chain results together
  sys.stderr.write('[%s] Generating LD in %d region(s)\n' % (time.asctime(),len(regions)))
  return [ iter(scan_ldpairs_region(region, maxd, rthreshold, dthreshold)) for region in regions ]


def scan_ldpairs_region(loci, maxd, rthreshold, dthreshold):
  '''
  A generator for pairs of loci within a specified genomic distance.
  Loci are assumed to be sorted by genomic location.
  '''
  # Scan each locus
  n = len(loci)
  for i in xrange(n):
    locus1    = loci[i]
    name1     = locus1.name
    genos1    = locus1.genos
    location1 = locus1.location

    yield name1,name1,1.0,1.0

    # And up to maxd distance beyond it
    for j in xrange(i+1,n):
      locus2    = loci[j]
      location2 = locus2.location

      if location2 - location1 > maxd:
        break

      counts = count_haplotypes(genos1, locus2.genos)
      r2,dprime = estimate_ld(*counts)

      if r2 >= rthreshold and abs(dprime) >= dthreshold:
        yield name1,locus2.name,r2,dprime


def filter_loci_by_maf(loci, minmaf, minobmaf, include):
  '''
  Generator that filters loci by a minimum MAF

  Loci come in two flavors, each with a distinct minimum MAF.
  If the locus.name is not in the provided include set, then
  the minmaf parameter is used as a threshold.  Otherwise, the
  minobmaf (minimum obligate MAF) threshold is applied.
  '''

  mafs = (minmaf,minobmaf)
  for locus in loci:
    # For mafs[locus.name in include], the index evaluates to:
    #    False == 0: Choose mafs[0] == minmaf
    #    True  == 1: Choose mafs[1] == minobmaf
    if locus.maf >= mafs[locus.name in include]:
      yield locus


def filter_loci_by_inclusion(loci, include):
  '''Generator that filters loci based on an inclusion set'''

  for locus in loci:
    if locus.name in include:
      yield locus


def filter_loci_by_hwp(loci, pvalue):
  '''
  Generator that filters loci based on significance of deviation from
  Hardy-Weinberg proportions
  '''
  for locus in loci:
    p = hwp_biallelic(locus.genos)
    if p >= pvalue:
      yield locus


range_all = (-sys.maxint,sys.maxint)

def filter_loci_by_range(loci, rangestring):
  '''Generator that filters loci based on an inclusion range'''

  ranges = []
  for range in rangestring.split(','):
    try:
      start,stop = range.split('-')
      start = int(start or -sys.maxint)
      stop  = int(stop  or  sys.maxint)

      if stop < start:
        start,stop = stop,start

    except (ValueError,TypeError):
      raise TagZillaError('ERROR: Invalid genomic range: %s' % range)

    ranges.append( (start,stop) )

  if range_all in ranges:
    ranges = [range_all]

  for locus in loci:
    for start,stop in ranges:
      if start <= locus.location < stop:
        yield locus
        break


def completion(locus):
  return sum(1 for g in locus if g),len(locus)


def filter_loci_by_completion(loci, mincompletion, mincompletionrate):
  '''Generator that filters loci by a minimum completion rate'''

  for locus in loci:
    m,n = completion(locus.genos)

    rate = 0
    if n:
      rate = float(m)/n

    if m >= mincompletion and rate >= mincompletionrate:
      yield locus


class Bin(set):
  INCLUDE_UNTYPED = -2
  INCLUDE_TYPED   = -1
  NORMAL          =  0
  EXCLUDE         =  2

  __slots__ = ('maf','disposition','maxcovered')
  def __init__(self, iterable=None, maf=None, disposition=NORMAL, maxcovered=None):
    if iterable is not None:
      set.__init__(self,iterable)
    else:
      set.__init__(self)
    self.maf = maf or 0.
    self.disposition = disposition
    self.maxcovered = max(maxcovered,len(self))

  def add(self, lname, maf):
    set.add(self, lname)
    self.maxcovered = max(self.maxcovered,len(self))
    self.maf += maf

  def remove(self, lname, maf):
    set.remove(self, lname)
    self.maf -= maf

  def discard(self, lname, maf):
    if lname in self:
      set.remove(self, lname)
      self.maf -= maf

  def average_maf(self):
    return float(self.maf)/len(self)

  def priority(self):
    return (self.disposition, -len(self), -self.maf)

  def __reduce__(self):
    return (Bin,(list(self),self.maf,self.disposition,self.maxcovered))

  def __repr__(self):
    return 'Bin(%s,%f,%d,%d)' % (list(self),self.maf,self.disposition,self.maxcovered)


def binldcmp(x,y):
  if x[0] == x[1] or y[0] == y[1]:
    return 1
  return -cmp(x[2],y[2]) or cmp(x[0],y[0]) or cmp(x[1],y[1])


result_priority_map = { 'obligate-untyped'  : -2,
                        'obligate-typed'    : -1,
                        'maximal-bin'       :  0,
                        'residual'          :  1,
                        'obligate-exclude'  :  2 }


class BinResult(object):
  __slots__ = ('binnum','tags','others','tags_required','average_maf','include',
               'ld','disposition','maxcovered','recommended_tags','include_typed')

  def sort(self):
    self.ld.sort(cmp=binldcmp)

  def priority(self):
    return (result_priority_map[self.disposition], -len(self), -self.average_maf)

  def __len__(self):
    return len(self.tags) + len(self.others)

  def __iter__(self):
    return chain(self.tags,self.others)

  def __le__(self, other):
    return self.priority() <= other.priority()


class BinStat(object):
  def __init__(self):
    self.count         = 0
    self.tags_required = 0
    self.loci          = 0
    self.width         = 0
    self.spacing       = 0
    self.total_tags    = 0
    self.others        = 0
    self.includes      = 0
    self.excludes      = 0

  def update(self, required, tags, others, width, spacing, include, excludes):
    self.count         += 1
    self.tags_required += required
    self.loci          += tags + others
    self.width         += width
    self.spacing       += spacing
    self.total_tags    += tags
    self.others        += others
    if include:
      self.includes += 1
    self.excludes += excludes

  def __add__(self, other):
    new = BinStat()
    new.count         = self.count         + other.count
    new.tags_required = self.tags_required + other.tags_required
    new.loci          = self.loci          + other.loci
    new.width         = self.width         + other.width
    new.spacing       = self.spacing       + other.spacing
    new.total_tags    = self.total_tags    + other.total_tags
    new.others        = self.others        + other.others
    new.includes      = self.includes      + other.includes
    new.excludes      = self.excludes      + other.excludes
    return new


class NullPairwiseBinOutput(object):
  def emit_bin(self, bin, qualifier, population, options):
    pass

  def emit_extra(self, lddata, tags, population):
    pass


class PairwiseBinOutput(NullPairwiseBinOutput):
  def __init__(self, outfile, exclude):
    self.outfile = outfile
    self.exclude = exclude
    outfile.writerow(['BIN','LNAME1','LNAME2','POPULATION','RSQUARED','DPRIME','DISPOSITION'])

  def emit_bin(self, bin, qualifier, population, options):
    outfile = self.outfile
    exclude = self.exclude
    bin.sort()

    for lname1,lname2,r2,dprime in bin.ld:
      if options.skip and (bin.disposition in ('obligate-exclude','residual')
                        or lname1 in exclude or lname2 in exclude):
        continue

      r2 = sfloat(r2)
      dprime = sfloat(dprime)
      disposition = pair_disposition(lname1, lname2, bin, qualifier)
      outfile.writerow([bin.binnum,lname1,lname2,population,r2,dprime,disposition])

  def emit_extra(self, lddata, tags, population):
    outfile = self.outfile
    bin = BinResult()
    bin.tags = set(tags)
    for (lname1,lname2),(r2,dprime) in lddata.iteritems():
      disposition = pair_disposition(lname1,lname2,bin,qualifier='interbin')
      r2     = sfloat(r2)
      dprime = sfloat(dprime)
      outfile.writerow(['',lname1,lname2,population,r2,dprime,disposition])


def save_ldpairs(filename, ldpairs):
  out = table_writer(filename,hyphen=sys.stdout)
  out.writerow(['LNAME1','LNAME2','RSQUARED','DPRIME'])

  def _write_pairs(pairs):
    for p in pairs:
      out.writerow(p)
      yield p

  return (_write_pairs(pairs) for pairs in ldpairs)


class NullLocusOutput(object):
  def emit_bin(self, bin, locusmap, qualifier, population):
    pass


class LocusOutput(NullLocusOutput):
  def __init__(self, locusinfofile, exclude):
    self.locusinfofile = locusinfofile
    self.exclude = exclude
    locusinfofile.writerow(LOCUS_HEADER3)

  def emit_bin(self, bin, locusmap, qualifier, population):
    locusinfofile = self.locusinfofile
    exclude = self.exclude
    for lname in chain(bin.tags,bin.others):
      disposition = locus_disposition(lname, bin, exclude, qualifier)
      l = locusmap[lname]
      maf = sfloat(l.maf)
      locusinfofile.writerow([l.name, l.chromosome, l.location, population,
                                   maf, bin.binnum, disposition])


class NullBinInfo(object):
  def emit_bin(self, bin, loci, exclude, population):
    pass

  def emit_summary(self, sumfile, population):
    pass

  def emit_multipop_summary(self, sumfile, ptags):
    pass


class BinInfo(NullBinInfo):
  dispositions = ['obligate-untyped','obligate-typed','maximal-bin','residual','obligate-exclude']

  def __init__(self, outfile, histomax):
    self.outfile = outfile
    self.stats = {}
    self.histomax = histomax

  def emit_bin(self, bin, loci, exclude, population):
    out = self.outfile

    binnum  = bin.binnum
    binsize = len(bin)
    amaf    = bin.average_maf*100
    locs    = sorted([ loci[lname].location for lname in bin ])
    spacing = sorted([ locs[i+1]-locs[i] for i in xrange(len(locs)-1) ])
    width   = locs[-1]-locs[0]
    excls   = exclude.intersection(bin) if exclude is not None else set()

    aspacing = 0
    if len(spacing) > 1:
      aspacing = mean(spacing)

    if bin.maxcovered == 1:
      hlen = 0
    else:
      hlen = min(self.histomax,binsize)

    stats = self.stats.get(population,None)
    if stats is None:
      stats = self.stats[population] = {}

    d = bin.disposition
    if d not in stats:
      stats[d] = [ BinStat() for i in xrange(self.histomax+1) ]

    stats[d][hlen].update(bin.tags_required, len(bin.tags), len(bin.others),
                          width, aspacing, bin.include is not None, len(excls))

    if not out:
      return

    population = population or 'default'
    out.write('Bin %-4d population: %s, sites: %d, tags %d, other %d, tags required %d, width %d, avg. MAF %.1f%%\n' \
                   % (binnum,population,binsize,len(bin.tags),len(bin.others),bin.tags_required,width,amaf))
    out.write('Bin %-4d Location: min %d, median %d, mean %d, max %d\n' \
                  % (binnum,locs[0],median(locs),mean(locs),locs[-1]))
    if len(spacing) > 1:
      out.write('Bin %-4d Spacing: min %d, median %d, mean %d, max %d\n' \
                    % (binnum,spacing[0],median(spacing),mean(spacing),spacing[-1]))
    out.write('Bin %-4d TagSnps: %s\n' % (binnum,' '.join(sorted(bin.tags))))
    if bin.recommended_tags:
      out.write('Bin %-4d RecommendedTags: %s\n' % (binnum, ' '.join(bin.recommended_tags)))
    out.write('Bin %-4d other_snps: %s\n' % (binnum,' '.join(sorted(bin.others))))

    if bin.include is not None:
      if bin.disposition == 'obligate-untyped':
        out.write('Bin %-4d Obligate_tag: %s, untyped\n' % (binnum,bin.include))
      else:
        out.write('Bin %-4d Obligate_tag: %s, typed\n' % (binnum,bin.include))

    if excls:
      out.write('Bin %-4d Excluded_as_tags: %s\n' % (binnum,' '.join(sorted(excls))))

    out.write('Bin %-4d Bin_disposition: %s\n' % (binnum,bin.disposition))
    out.write('Bin %-4d Loci_covered: %s\n' % (binnum,bin.maxcovered))
    out.write('\n')


  def emit_summary(self, sumfile, population):
    out = sumfile
    stats = self.stats.get(population,{})

    tstats = {}
    for d in self.dispositions:
      if d in stats:
        self.emit_summary_stats(out, stats[d], d, population)
        tstats[d] = sum(stats[d], BinStat())

    if not population:
      out.write('\nBin statistics by disposition:\n')
    else:
      out.write('\nBin statistics by disposition for population %s:\n' % population)

    out.write('                      tags                                total   non-     avg    avg\n')
    out.write(' disposition          req.   bins     %    loci      %    tags    tags    tags  width\n')
    out.write(' -------------------- ------ ------ ------ ------- ------ ------- ------- ---- ------\n')

    total_bins = sum(s.count for s in tstats.values())
    total_loci = sum(s.loci  for s in tstats.values())

    for d in self.dispositions:
      self.emit_summary_line(out, '%-20s' % d, tstats.get(d,BinStat()), total_bins, total_loci)

    self.emit_summary_line(out, '              Total ', sum(tstats.values(), BinStat()), total_bins, total_loci)
    out.write('\n')
    out.flush()


  def emit_multipop_summary(self, sumfile, tags):
    n = sum(tags.itervalues())

    sumfile.write('\nTags required by disposition for all populations:\n')

    sumfile.write('                      tags         \n')
    sumfile.write(' disposition          req.     %   \n')
    sumfile.write(' -------------------- ------ ------\n')

    for d in self.dispositions:
      m = tags.get(d,0)
      sumfile.write(' %-20s %6d %6.2f\n' % (d,m,percent(m,n)))

    sumfile.write('              Total   %6d %6.2f\n\n' % (n, 100))
    sumfile.flush()


  def emit_summary_stats(self, out, stats, disposition, population):
    if not population:
      out.write('\nBin statistics by bin size for %s:\n\n' % disposition)
    else:
      out.write('\nBin statistics by bin size for %s for population %s:\n\n' % (disposition,population))

    out.write(' bin   tags                                total   non-     avg    avg\n')
    out.write(' size  req.   bins     %    loci      %    tags    tags    tags  width\n')
    out.write(' ----- ------ ------ ------ ------- ------ ------- ------- ---- ------\n')
    total_bins = sum(s.count for s in stats)
    total_loci = sum(s.loci  for s in stats)

    hlist = [ i for i,s in enumerate(stats) if s.count ]
    hmin = min(hlist)
    hmax = max(hlist)

    for i in xrange(hmin,hmax+1):
      if not i:
        label = 'singl'
      elif i == self.histomax:
        label = '>%2d  ' % (i-1)
      else:
        label = '%3d  ' % i

      self.emit_summary_line(out, label, stats[i], total_bins, total_loci)

    self.emit_summary_line(out, 'Total', sum(stats, BinStat()), total_bins, total_loci)
    out.write('\n')


  def emit_summary_line(self, out, label, stats, total_bins, total_loci):
    n = stats.count
    m = stats.loci
    if n:
      t = float(stats.total_tags) / n
      w = float(stats.width) / n
    else:
      t,w = 0,0
    out.write(' %s %6d %6d %6.2f %7d %6.2f %7d %7d %4.1f %6d\n' % (label,
              stats.tags_required,n,percent(n,total_bins),
              m,percent(m,total_loci),stats.total_tags,stats.others,t,w))


def locus_result_sequence(filename, locusmap, exclude):
  '''
  Returns a generator of BinResult objects for the tagzilla locus output
  file name.

  The locusmap dictionary and the exclude set are filled in incrementally as
  the result stream is processed.

  NOTE: This function is not currently used by tagzilla -- rather it exists
        to unparse tagzilla output and is included as a utility function
        for when tagzilla is used as a module.
  '''
  locusfile = table_reader(filename)

  header = locusfile.next()

  if header == PAIR_HEADER:
    version = 0
    grouper = itemgetter(0,3)
  elif header == LOCUS_HEADER1:
    version = 1
    grouper = itemgetter(3)
  elif header == LOCUS_HEADER2:
    version = 2
    grouper = itemgetter(4,2)
  elif header == LOCUS_HEADER3:
    version = 3
    grouper = itemgetter(5,3)
  else:
    raise TagZillaError('ERROR: Invalid input format for file %s.' % filename)


  for binnum,(_,loci) in enumerate(groupby(locusfile,grouper)):
    bin = BinResult()
    bin.binnum = binnum

    bin.tags = []
    bin.others = []
    bin.ld = []
    bin.include = None
    bin.average_maf = 0
    bin.maxcovered = 0
    bin.recommended_tags = []
    bin.disposition = 'maximal-bin'
    bin.tags_required = 1

    for locus in loci:
      if version == 0:
        if locus[1] != locus[2]:
          continue
        lname = locus[1]
        location = 0
        population = locus[3]
        maf = 0.5
        disposition = locus[6]
      elif version == 1:
        lname,location,maf,binnum,disposition = locus
        population = chr = ''
      elif version == 2:
        lname,location,population,maf,binnum,disposition = locus
        chr = ''
      elif version == 3:
        lname,chr,location,population,maf,binnum,disposition = locus

      bin.binnum = binnum
      locus = locusmap[lname] = Locus(lname, chr, int(location), [])

      maf = float(maf)
      locus.maf = maf
      bin.average_maf += maf

      disposition = disposition.split(',')

      if 'other' in disposition:
        bin.others.append(lname)
      elif 'exclude' in disposition:
        bin.others.append(lname)
        exclude.add(lname)
      elif 'excluded-tag' in disposition:
        bin.tags.append(lname)
        bin.disposition = 'obligate-exclude'
      elif 'obligate-tag' in disposition:
        bin.tags.append(lname)
        bin.disposition = 'obligate-include'
        bin.include = lname
      elif 'untyped-tag' in disposition:
        bin.tags.append(lname)
        bin.disposition = 'obligate-untyped'
        bin.include = lname
      elif 'typed-tag' in disposition:
        bin.tags.append(lname)
        bin.disposition = 'obligate-typed'
        bin.include = lname
      elif 'alternate-tag' in disposition:
        bin.tags.append(lname)
      elif 'candidate-tag' in disposition or \
           'necessary-tag' in disposition:
        bin.tags.append(lname)
      elif 'lonely-tag' in disposition:
        bin.tags.append(lname)
        bin.maxcovered = 2
      elif 'singleton-tag' in disposition:
        bin.tags.append(lname)

      if 'recommended' in disposition:
        bin.recommended_tags.append(lname)

      if 'untyped_bin' in disposition:
        bin.disposition = 'obligate-untyped'

      if 'typed_bin' in disposition:
        bin.disposition = 'obligate-typed'

      if 'residual' in disposition:
        bin.disposition = 'residual'

    bin.maxcovered = max(bin.maxcovered,len(bin))
    bin.average_maf /= len(bin)

    yield population,bin


def must_split_bin(bin, binsets, get_tags_required):
  if not get_tags_required:
    return False

  tags_required = get_tags_required(len(bin))

  if tags_required == 1:
    return False

  tags = [ lname for lname in bin if can_tag(binsets[lname],bin) ]

  return len(tags) < tags_required <= len(bin)


class NaiveBinSequence(object):
  def __init__(self, loci, binsets, lddata, get_tags_required):
    self.loci = loci
    self.binsets = binsets
    self.lddata = lddata
    self.get_tags_required = get_tags_required

  def __iter__(self):
    return self

  def pop(self):
    if not self.binsets:
      raise StopIteration

    while 1:
      ref_lname = self.peek()
      largest = self.binsets[ref_lname]

      if not must_split_bin(largest, self.binsets, self.get_tags_required):
        break

      self.split_bin(ref_lname,largest)

    # Remove all references to this locus from any other
    # binset not captured
    bins = {}
    for lname in largest:
      bin = self.pop_bin(lname)
      bins[lname] = bin
      maf = self.loci[lname].maf
      for lname2 in bin - largest:
        self.reduce_bin(lname2, lname, maf)

    return ref_lname,largest,bins

  next = pop

  def peek(self):
    '''
    Find the largest bin among all the sets, selecting bins ordered by
    inclusion status, size (descending), and MAF (descending).  This
    implementation is mainly for demonstration and testing, as it uses a
    naive and potentially very slow linear search.  See the
    FastBinSequence descendant class for a more efficient solution based
    on a priority queue.
    '''
    bins = self.binsets.iteritems()
    # First set the first bin as the best
    lname, bin = bins.next()
    prio = bin.priority()
    # Then iterate through the remaining items to refine the result
    for current_lname,bin in bins:
      current_prio = bin.priority()
      if current_prio < prio:
        lname = current_lname
        prio = current_prio

    return lname

  def pop_bin(self, lname):
    return self.binsets.pop(lname)

  def reduce_bin(self, other_locus, taken_locus, maf):
    if other_locus in self.binsets:
      self.binsets[other_locus].discard(taken_locus, maf)

  def split_bin(self, ref_lname, bin):
    ld = []
    for lname in bin:
      if lname == ref_lname:
        continue

      covered = len(self.binsets.get(lname,[]))

      lname1,lname2 = ref_lname,lname
      if (lname1,lname2) not in self.lddata:
        lname1,lname2=lname2,lname1

      r2,dprime = self.lddata[lname1,lname2]
      ld.append( (-covered,r2,lname) )

    ld.sort()

    # Remove smallest ld value
    covered,r2,lname = ld[0]
    self.reduce_bin(ref_lname, lname, self.loci[lname].maf)
    self.reduce_bin(lname, ref_lname, self.loci[ref_lname].maf)


class FastBinSequence(NaiveBinSequence):
  def __init__(self, loci, binsets, lddata, get_tags_required):
    NaiveBinSequence.__init__(self, loci, binsets, lddata, get_tags_required)

    import pqueue
    self.pq = pq = pqueue.PQueue()

    for lname,bin in binsets.iteritems():
      pq[lname] = bin.priority()

  def peek(self):
    priority,ref_lname = self.pq.peek()
    return ref_lname

  def pop_bin(self, lname):
    del self.pq[lname]
    return NaiveBinSequence.pop_bin(self, lname)

  def reduce_bin(self, other_locus, taken_locus, maf):
    NaiveBinSequence.reduce_bin(self, other_locus, taken_locus, maf)
    if other_locus in self.binsets:
      self.pq[other_locus] = self.binsets[other_locus].priority()


def BinSequence(loci, binsets, lddata, get_tags_required):
  try:
    return FastBinSequence(loci, binsets, lddata, get_tags_required)
  except ImportError:
    pass

  return NaiveBinSequence(loci, binsets, lddata, get_tags_required)


class NaiveMultiBinSequence(object):
  def __init__(self, loci, binsets, lddata, get_tags_required):
    self.loci    = loci
    self.binsets = binsets
    self.lddata  = lddata
    self.get_tags_required = get_tags_required

    self.lnames = set()
    for pop_binsets in binsets:
      self.lnames.update(pop_binsets)

  def __iter__(self):
    return self

  def pop(self):
    if not any(self.binsets):
      raise StopIteration

    while 1:
      ref_lname = self.peek()

      split = False
      for pop_binsets,pop_loci in izip(self.binsets,self.loci):
        bin = pop_binsets.get(ref_lname,None)
        if bin and must_split_bin(bin, pop_binsets, self.get_tags_required):
          self.split_bin(pop_binsets, pop_loci, ref_lname,bin)
          split = True
          break

      if not split:
        break

    # Remove all references to this locus from any other
    # binset not captured
    largest = []
    bins = []
    touched = set()
    for pop_binsets,pop_loci in izip(self.binsets,self.loci):
      used_bins = {}
      lbin = pop_binsets.get(ref_lname,None)
      if lbin:
        touched.update(lbin)
        for lname in lbin:
          bin = pop_binsets.pop(lname)
          used_bins[lname] = bin
          maf = pop_loci[lname].maf
          for lname2 in bin - lbin:
            self.reduce_bin(pop_binsets, lname2, lname, maf)

      largest.append(lbin)
      bins.append(used_bins)

    self.update_bins(touched)
    return ref_lname,largest,bins

  next = pop

  def peek(self):
    lnames = iter(self.lnames)

    prio = None
    while prio is None:
      lname = lnames.next()
      prio  = self.priority(lname)

    # Then iterate through the remaining items to refine the result
    for current_lname in lnames:
      current_prio = self.priority(current_lname)
      if current_prio is not None and current_prio < prio:
        lname = current_lname
        prio  = current_prio

    return lname

  def reduce_bin(self, binsets, other_locus, taken_locus, maf):
    if other_locus in binsets:
      binsets[other_locus].discard(taken_locus, maf)

  def priority(self, lname):
    disposition = 1000
    binlen = maf = pops = 0
    minlen = sys.maxint

    for pop_binsets in self.binsets:
      bin = pop_binsets.get(lname)
      if bin:
        disposition = min(disposition, bin.disposition)
        minlen  = min(minlen,len(bin))
        binlen += len(bin)
        maf    += bin.maf
        pops   += 1

    if minlen == 1:
      pops   *= 2
      binlen *= 2

    if disposition < 1000:
      return (disposition,-pops,-binlen,-maf)
    else:
      return None

  def split_bin(self, binsets, loci, ref_lname, bin):
    ld = []
    for lname in bin:
      if lname == ref_lname:
        continue

      covered = len(binsets.get(lname,[]))

      lname1,lname2 = ref_lname,lname
      if (lname1,lname2) not in self.lddata:
        lname1,lname2=lname2,lname1

      r2,dprime = self.lddata[lname1,lname2]
      ld.append( (-covered,r2,lname) )

    ld.sort()

    # Remove smallest ld value
    covered,r2,lname = ld[0]
    self.reduce_bin(binsets, ref_lname, lname, loci[lname].maf)
    self.reduce_bin(binsets, lname, ref_lname, loci[ref_lname].maf)
    self.update_bins([ref_lname,lname])

  def update_bins(self,lnames):
    pass


class FastMultiBinSequence(NaiveMultiBinSequence):
  def __init__(self, loci, binsets, lddata, get_tags_required):
    NaiveMultiBinSequence.__init__(self, loci, binsets, lddata, get_tags_required)

    import pqueue
    self.pq = pq = pqueue.PQueue()

    for lname in self.lnames:
      pq[lname] = self.priority(lname)

  def peek(self):
    priority,ref_lname = self.pq.peek()
    return ref_lname

  def update_bins(self,lnames):
    for lname in lnames:
      p = self.priority(lname)
      if p is not None:
        self.pq[lname] = p
      else:
        del self.pq[lname]


def MultiBinSequence(loci, binsets, lddata, get_tags_required):
  try:
    return FastMultiBinSequence(loci, binsets, lddata, get_tags_required)
  except ImportError:
    pass

  return NaiveMultiBinSequence(loci, binsets, lddata, get_tags_required)


def build_binsets(loci, ldpairs, includes, exclude, designscores):
  '''
  Build initial data structures:
    binsets: Dictionary that for each locus, stores the set of all other
             loci that meet the rthreshold and the sum of the MAFs
    lddata:  Dictionary of locus pairs to r-squared values
  '''

  binsets = {}
  lddata  = {}

  for pairs in ldpairs:
    for lname1,lname2,r2,dprime in pairs:
      if lname1 not in binsets:
        binsets[lname1] = Bin([lname1], loci[lname1].maf)
      if lname2 not in binsets:
        binsets[lname2] = Bin([lname2], loci[lname2].maf)

      if lname1 != lname2:
        binsets[lname1].add(lname2, loci[lname2].maf)
        lddata[lname1,lname2] = r2,dprime
        binsets[lname2].add(lname1, loci[lname1].maf)

  # Update the bin disposition if the lname is one of the excludes
  for lname in exclude or []:
    if lname in binsets:
      binsets[lname].disposition = Bin.EXCLUDE

  # Update the bin disposition for all undesignable loci, if design scores
  # are provided
  if designscores:
    for lname,bin in binsets.iteritems():
      if designscores[lname] < epsilon:
        bin.disposition = bin.EXCLUDE
        exclude.add(lname)

  # Build include sets and pre-remove obligates from other include bins
  for lname in includes.untyped:
    if lname in binsets:
      bin = binsets[lname]
      bin.disposition = Bin.INCLUDE_UNTYPED

      for lname2 in includes.untyped & bin:
        if lname != lname2:
          binsets[lname].remove(lname2, loci[lname2].maf)

  for lname in includes.typed:
    if lname in binsets:
      bin = binsets[lname]
      bin.disposition = Bin.INCLUDE_TYPED

  return binsets,lddata


def bin_qualifier(bin, binned_loci, options):
  qualifier = ''
  if ((options.targetbins and  bin.binnum > options.targetbins) or \
      (options.targetloci and binned_loci > options.targetloci)) and \
      bin.disposition != 'obligate-exclude':
    qualifier = 'residual'
    bin.disposition = 'residual'
  elif bin.disposition == 'obligate-exclude':
    qualifier = 'excluded'
  elif bin.disposition == 'obligate-typed':
    qualifier = 'typed_bin'
  elif bin.disposition == 'obligate-untyped':
    qualifier = 'untyped_bin'
  return qualifier


def tag_disposition(lname, bin):
  if bin.disposition == 'obligate-untyped':
    if lname == bin.include:
      disposition = 'untyped-tag'
    elif lname in bin.include_typed:
      disposition = 'redundant-tag'
    else:
      disposition = 'alternate-tag'
  elif bin.disposition == 'obligate-typed':
    if lname == bin.include:
      disposition = 'typed-tag'
    elif lname in bin.include_typed:
      disposition = 'redundant-tag'
    else:
      disposition = 'alternate-tag'
  elif bin.disposition == 'obligate-exclude':
    disposition = 'excluded-tag'
  elif len(bin.tags) > 1:
    disposition = 'candidate-tag'
  elif len(bin) > 1:
    disposition = 'necessary-tag'
  elif bin.maxcovered > 1:
    disposition = 'lonely-tag'
  else:
    disposition = 'singleton-tag'

  if lname in bin.recommended_tags:
    disposition += ',recommended'

  return disposition


def locus_disposition(lname, bin, exclude, qualifier=None):
  if lname in bin.tags:
    disposition = tag_disposition(lname, bin)
  elif exclude is not None and lname in exclude and bin.disposition != 'obligate-exclude':
    disposition = 'exclude'
  else:
    disposition = 'other'

  if qualifier:
    disposition = '%s,%s' % (disposition,qualifier)

  return disposition


def pair_disposition(lname1, lname2, bin, qualifier=None):
  if lname1 == lname2:
    disposition = tag_disposition(lname1, bin)
  else:
    labels = ['other','tag']
    # Don't ask -- you really don't want to know.  Sigh.
    disposition = '%s-%s' % (labels[lname1 in bin.tags], labels[lname2 in bin.tags])

  if qualifier:
    disposition = '%s,%s' % (disposition,qualifier)

  return disposition


def sfloat(n):
  '''Compactly format a float as a string'''
  return ('%.3f' % n).rstrip('0.').lstrip('0') or '0'


def can_tag(bin, reference):
  '''
  Return True if the specified candidate bin can tag the reference bin, False otherwise.

  The following conditions must hold:
    1) If the reference bin disposition is not an exclude, then the candidate bin cannot be either.
    2) The contents of the candidate bin set must be a superset of the reference bin.
  '''
  return (bin.disposition != Bin.EXCLUDE or reference.disposition == Bin.EXCLUDE) \
         and bin.issuperset(reference)


def build_result(lname, largest, bins, lddata, includes, get_tags_required):
  result = BinResult()

  if get_tags_required:
    result.tags_required = get_tags_required(len(largest))
  else:
    result.tags_required = 1

  result.recommended_tags = []
  result.include_typed    = includes.typed & largest
  result.average_maf      = largest.average_maf()
  result.maxcovered       = largest.maxcovered

  if largest.disposition in (Bin.INCLUDE_TYPED,Bin.INCLUDE_UNTYPED):
    result.include = lname
  else:
    result.include = None

  result.tags   = []
  result.others = []
  result.ld     = []

  if largest.disposition == Bin.INCLUDE_UNTYPED:
    result.disposition = 'obligate-untyped'
  elif largest.disposition == Bin.INCLUDE_TYPED:
    result.disposition = 'obligate-typed'
  elif largest.disposition == Bin.EXCLUDE:
    result.disposition = 'obligate-exclude'
  else:
    result.disposition = 'maximal-bin'

  # Process each locus in the bin
  for lname,bin in bins.iteritems():
    # If the current bin is a superset of the reference set, then consider this locus
    # a tag.  The superset is needed to handle the case where the reference locus is
    # an obligate include and may not be the largest bin.
    if can_tag(bin, largest):
      result.tags.append(lname)

      # Update maximum coverage number for all candidate tags, except for
      # obligate include bins
      if largest.disposition not in (Bin.INCLUDE_TYPED,Bin.INCLUDE_UNTYPED):
        result.maxcovered = max(result.maxcovered,bin.maxcovered)
    else:
      result.others.append(lname)

  assert len(result.tags) >= result.tags_required

  for lname in result.tags:
    # Output the tags as self-pairs (r-squared=1,dprime=1)
    result.ld.append( (lname,lname,1.,1.) )

  # For each pair of loci in the bin, yield name, location, and LD info
  for lname1,lname2 in pair_generator(largest):
    if (lname1,lname2) not in lddata:
      lname1,lname2=lname2,lname1

    if (lname1,lname2) in lddata:
      r2,dprime = lddata.pop( (lname1,lname2) )
      result.ld.append( (lname1,lname2,r2,dprime) )

  return result


def binner(loci, binsets, lddata, includes, get_tags_required=None):
  '''
  Greedy tag marker binning algorithm -- similar to the Carlson algorithm.

  The binner utilizes the recomputed binsets and lddata which reflect LD
  values previously filtered against the chosen minimum threshold for bin.

  Given a set of binsets and a sequence of loci, the binner iteratively selects
  the largest bin in the following priority order:
    1) bin dispositions such that include > normal > exclude
    2) bins with the largest number of loci
    3) bins with the largest total MAF

  The binner constructs a BinResult object for each bin with the following
  attributes:
    1) tags: a list of all loci chosen as tags
    2) others: a list of all non-tag loci
    3) average_maf: the average of MAF of all loci in the bin
    4) include: the locus name of the obligate tag or None
    5) disposition: This attribute may take one of three possible values:
           'obligate-untyped' if the reference locus in include set and untyped
           'obligate-typed'   if the reference locus in include set and typed
           'obligate-exclude' if the reference locus in exclude set
           'maximal-bin' otherwise
    6) ld: a list of tuples of pairwise LD data within each bin.  Tags for
           that bin are encoded as records with the locus paired with itself
           and r-squared and dprime of 1:

               (BINNUM,LNAME,LNAME,1,1,DISPOSITION)

               DISPOSITION for tags takes one of the following values:

                 'untyped-tag'      if it is an obligate tag and has not been genotyped
                 'typed-tag'        if it is an obligate tag and is the reference tag of
                                    an obligate-typed bin
                 'redundant-tag'    an obligate tag and has been genotyped, but not the
                                    reference tag of a bin
                 'alternate-tag'    if it is a tag in an obligate-include
                                    bin, but not the obligate tag
                 'excluded-tag'     a tag for a bin that contains all
                                    obligatorily excluded loci
                 'candidate-tag'    a tag for a bin that has more than one
                                    possible non-obligate tag
                 'necessary-tag'    a tag for a bin that has only one tag,
                                    but covers at least one other locus
                 'lonely-tag'       a tag for a bin with no other loci, but
                                    originally covered more loci.  These
                                    additional loci where removed by
                                    previous iterations of the binning
                                    algorithm.  This disposition is
                                    primarily to distinguish these bins from
                                    singletons, which intrinsically are in
                                    insufficient Ld with any other locus.
                 'singleton-tag'    a tag that is not in significant LD with
                                    any other locus, based on the specified
                                    LD criteria

           The pairwise r-squared values within the bin follow in the form:
               (BINNUM,LNAME1,LNAME2,RSQUARED,DPRIME,DISPOSITION).

               DISPOSITION for these pairwise elements takes one of the
               following values:

                 'tag-tag'         for LD between tags within a bin;
                 'other-tag'       for LD between a non-tag and a tag
              &  'tag-other'         or a tag and a non-tag, respectively;
                 'other-other'     for LD between non-tags within a bin.

    7) maxcovered: the maximum number for loci that each candidate tag may
                   have covered.  The actual coverage may be smaller, due to
                   other loci removed by bins that were selected in an
                   earlier iteration of the algorithm.  For obligate include
                   bins, only the obligatory tags is considered, since
                   alternate tags are not considered.

  @type:  binsets: dictionary, key is of string type and value is a Bin object
  @param: binsets: dictionary that for each locus mapped to a Bin which includes
                   all loci satisfying rsquared threshold with this reference locus and its
                   own MAF greater than MAF threshold
  @type      loci: Sequence of (LNAME, LOCATION,MAF,...GENOTYPES...)
  @param     loci: A sequence of loci that may appear in the ldpairs

  @type    lddata: Sequence of (LNAME1,LNAME2,R-SQUARED,DPRIME)
  @param   lddata: A sequence of pairwise LD information that exceed a given
                   threshold.  i.e, they must be pre-filtered by the r-squared
                   criteria before being passed to the binner.
  @rtype:          generator for an ordered sequence of BinResult objects
  @return:         the optimal bins with tags, others, ld information etc. for each
  '''

  bin_sequence = BinSequence(loci, binsets, lddata, get_tags_required)

  # Run while we still have loci to bin
  for lname,largest,bins in bin_sequence:
    yield build_result(lname, largest, bins, lddata, includes, get_tags_required)


def binner_vector(loci, binsets, lddata, includes, get_tags_required=None):
  bin_sequence = MultiBinSequence(loci, binsets, lddata, get_tags_required)

  for lname,largest,bins in bin_sequence:
    results = []
    for pop_largest,pop_bins,pop_lddata in izip(largest,bins,lddata):
      if pop_largest:
        result = build_result(lname, pop_largest, pop_bins, pop_lddata, includes, get_tags_required)
        results.append(result)
      else:
        results.append(None)

    yield lname,results


def generate_ldpairs_vector(options, include, subset, ldsubset):
  labels = get_populations(options.multipopulation)
  pops = len(labels)
  regions = len(options.genotypes) // pops

  if len(options.genotypes) % pops != 0:
    raise TagZillaError('ERROR: The number of input files must be a multiple of the number of populations')

  for i in xrange(regions):
    ldpairs = []
    multi_options = []
    locusmap = []

    for filename in options.genotypes[i*pops:(i+1)*pops]:
      lmap    = {}
      regions = generate_ldpairs_from_file(filename, lmap, include, subset, ldsubset, options)
      pairs   = chain(*regions)

      ldpairs.append(pairs)
      locusmap.append(lmap)

    yield ldpairs,locusmap


def do_tagging_vector(ldpairs, includes, exclude, designscores, options):
  labels = get_populations(options.multipopulation)
  pops = len(labels)

  # If we require a total ordering, then build binsets from all ldpairs
  if options.targetbins or options.targetloci:
    sys.stderr.write('[%s] Building global binsets\n' % time.asctime())
    multi_ldpairs  = [ [] for p in xrange(pops) ]
    multi_locusmap = [ {} for p in xrange(pops) ]

    for region,lmap in ldpairs:
      for pop_ldpairs,pop_locusmap,pop_lmap,pairs in izip(multi_ldpairs,multi_locusmap,lmap,region):
        pop_ldpairs.append(pairs)
        update_locus_map(pop_locusmap,pop_lmap.itervalues())

    binsets = []
    lddata  = []
    for pop_ldpairs,pop_locusmap in izip(multi_ldpairs,multi_locusmap):
      pop_binsets,pop_lddata = build_binsets(pop_locusmap, pop_ldpairs, includes, exclude, designscores)
      binsets.append(pop_binsets)
      lddata.append(pop_lddata)

    sys.stderr.write('[%s] Choosing global bins\n' % time.asctime())
    bins = binner_vector(multi_locusmap, binsets, lddata, includes, get_tags_required_function(options))
    yield bins,lddata,multi_locusmap

  else:
    # Otherwise, process each sequence of ldpairs independently
    for pairs,locusmap in ldpairs:
      sys.stderr.write('[%s] Building binsets\n' % time.asctime())

      binsets = []
      lddata  = []
      for pop_ldpairs,pop_locusmap in izip(pairs,locusmap):
        pop_binsets,pop_lddata = build_binsets(pop_locusmap, [pop_ldpairs], includes, exclude, designscores)
        binsets.append(pop_binsets)
        lddata.append(pop_lddata)

      sys.stderr.write('[%s] Choosing bins\n' % time.asctime())
      bins = binner_vector(locusmap, binsets, lddata, includes, get_tags_required_function(options))
      yield bins,lddata,locusmap


def load_festa_file(filename, locusmap, subset, rthreshold):
  '''
  Load FESTA formatted file that contain pre-computed LD data for pairs of loci
  '''
  ldfile = autofile(filename)
  header = ldfile.readline()

  for line in ldfile:
    lname1,lname2,ldvalue = re_spaces.split(line.strip())
    ldvalue = float(ldvalue)

    if subset is not None and (lname1 not in subset or lname2 not in subset):
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

    if subset is not None and (lname1 not in subset or lname2 not in subset):
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


def build_design_score(designscores,designdefault=0):
  designscores = designscores or []
  aggscores = defaultdict(lambda: designdefault)
  for design in designscores:
    design = design.split(':')
    dfile = design[0]
    threshold = 0.0
    scale = 1.0
    if len(design) > 1:
      threshold = float(design[1])
    elif len(design) > 2:
      scale = float(design[2])
    scores = read_design_score(dfile)
    for lname,score in scores:
      if score < threshold:
        score = 0.0
      aggscores[lname] = aggscores.get(lname,1.0)*score*scale
  return aggscores


def build_tag_criteria(tagcriteria):
  weights = {}
  tagcriteria = tagcriteria or []
  for c in tagcriteria:
    c = c.lower().split(':')
    method = c[0]
    weight = 2 #default weight
    if len(c) > 1:
      weight = float(c[1])
    weights[method] = weight
  return weights


class TagSelector(object):
  default_weight = 2.0

  def __init__(self, scores, weights):
    self.scores  = scores
    self.weights = weights

  def select_tags(self,bin):
    if not self.weights and not self.scores:
      return

    if not self.weights and bin.disposition == 'obligate-exclude':
      return

    if len(bin.tags) == 1:
      bin.recommended_tags = list(bin.tags)[:1]
      return

    weights = {}
    for method,weight in self.weights.iteritems():
      w = self.build_weights(bin, method, weight)
      for lname, weight in w.iteritems():
        weights[lname] = weights.get(lname,1) * weight

    if bin.disposition == 'obligate-exclude':
      allscores = {}
    else:
      allscores = self.scores

    # Default score: 0 if scores exist, 1 otherwise
    default_score = not allscores

    scores = []
    for tag in bin.tags:
      score  = allscores.get(tag,default_score)
      weight = weights.get(tag,1)
      s = score * weight
      scores.append( (s,tag) )

    scores.sort(reverse=True)

    # Store tags in weight order
    bin.tags = [lname for s, lname in scores]

    # Recommend tags
    bin.recommended_tags = bin.tags[:bin.tags_required]
    if bin.include is not None and bin.include not in bin.recommended_tags:
      bin.recommended_tags = [bin.include]+bin.recommended_tags[:bin.tags_required-1]

  def build_weights(self, bin, method, weight):
    if not method:
      return {}

    # Lexically-scoped weighting functions
    def maxsnp():
      w[lname1] = min(w.get(lname1,1),r2)
    def avgsnp():
      w[lname1] = w.get(lname1,0) + r2
    def maxtag():
      if lname2 not in bin.tags:
        maxsnp()
    def avgtag():
      if lname2 not in bin.tags:
        avgsnp()

    func = locals().get(method.lower(),None)

    if not callable(func):
      raise RuntimeError('Invalid tag information criterion specified')

    w = {}
    for lname1,lname2,r2,dprime in bin.ld:
      if lname1==lname2:
        continue
      for lname1,lname2 in [(lname1,lname2),(lname2,lname1)]:
        if lname1 in bin.tags:
          func()

    if not w:
      return {}

    maxval = max(w.itervalues())

    weights = {}
    for tag in bin.tags:
      if abs(w[tag] - maxval) > 1e-10:
        weights[tag] = 1./weight

    return weights


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


def load_genotypes(filename, options):
  loci = load_genostream(filename,format=options.informat,genorepr=options.ingenorepr,
                                  genome=options.loci,phenome=options.pedigree,
                                  transform=options).as_ldat()

  nonfounders = set(f.name for f in loci.phenome.phenos.itervalues() if f.nonfounder())
  loci = loci.transformed(excludesamples=nonfounders,repack=True)

  genome = loci.genome
  for lname,genos in loci:
    loc = genome.get_locus(lname)
    yield Locus(lname,loc.chromosome,loc.location,genos)


def filter_loci(loci, include, subset, options):
  if getattr(options,'obmaf',None) is None:
    options.obmaf = options.maf

  if options.maf or options.obmaf:
    loci = filter_loci_by_maf(loci, options.maf, options.obmaf, include)

  if subset is not None:
    loci = filter_loci_by_inclusion(loci, subset)

  if options.range:
    loci = filter_loci_by_range(loci, options.range)

  if options.mincompletion or options.mincompletionrate:
    loci = filter_loci_by_completion(loci, options.mincompletion, options.mincompletionrate)

  if options.hwp:
    loci = filter_loci_by_hwp(loci, options.hwp)

  return loci


def order_loci(loci):
  counter = count().next
  def locus_key(l):
    return l.chromosome,l.location,counter()

  return sorted(loci,key=locus_key)


def filter_loci_ldsubset(loci, ldsubset, maxd):
  '''
  Yield only loci where there exists at least one monitored location within
  maxd distance.

  This implementation performs two linear passes over the list of loci, one
  forward and one in reverse.  As such, a locus may be yielded twice if it
  is within maxd if a monitored location on both the left and the right.
  '''
  if ldsubset is None:
    return loci

  monitor = [ (l.chromosome,l.location) for l in loci if l.name in ldsubset ]

  n = len(loci)
  keep = set()

  # Scan forward though loci, yielding all following loci within maxd of a
  # monitored location
  pos = 0
  for chr,loc in monitor:
    while pos < n and (loci[pos].chromosome,loci[pos].location) < (chr,loc):
      pos += 1
    while pos < n and loci[pos].chromosome==chr and loci[pos].location-loc <= maxd:
      keep.add(loci[pos])
      pos += 1

  # Scan backward though loci, yielding all prior loci within maxd if a
  # monitored location
  pos = n-1
  for chr,loc in reversed(monitor):
    while pos >= 0 and (loci[pos].chromosome,loci[pos].location) > (chr,loc):
      pos -= 1
    while pos >= 0 and loci[pos].chromosome==chr and loc-loci[pos].location <= maxd:
      keep.add(loci[pos])
      pos -= 1

  return order_loci(keep)


def update_locus_map(locusmap, loci):
  addloci = dict( (locus.name,locus) for locus in loci )
  overlap = set(addloci).intersection(locusmap)
  if overlap:
    raise TagZillaError('ERROR: Genotype files may not contain overlapping loci')
  locusmap.update(addloci)


def generate_ldpairs(options, locusmap, include, subset, ldsubset):
  for filename in options.genotypes:
    regions = generate_ldpairs_from_file(filename, locusmap, include, subset, ldsubset, options)

    for ldpairs in regions:
      yield ldpairs

    # Clear locusmap to save memory
    # FIXME: May need to be disabled due to potentially wonky semantics
    locusmap.clear()


def generate_ldpairs_from_file(filename, locusmap, include, subset, ldsubset, options):
  sys.stderr.write('[%s] Processing input file %s\n' % (time.asctime(),filename))

  if options.informat == 'festa':
    return load_festa_file(filename, locusmap, subset, options.r)

  elif options.informat == 'hapmapld':
    return load_hapmapld_file(filename, locusmap, subset, options.maxdist*1000, options.r, options.d)

  else: # generate from genotype file
    loci = load_genotypes(filename,options)
    loci = filter_loci(loci, include, subset, options)
    loci = order_loci(loci)
    loci = filter_loci_ldsubset(loci, ldsubset, options.maxdist*1000)

    # Locusmap must contain only post-filtered loci
    update_locus_map(locusmap, loci)
    return scan_ldpairs(loci, options.maxdist*1000, options.r, options.d)


def get_populations(option):
  if not option:
    return ['']

  try:
    n = int(option)
    labels = [ str(i) for i in xrange(1,n+1) ]

  except ValueError:
    labels = [ l.strip() for l in option.split(',') if l.strip() ]

  return labels


def get_tags_required_function(options):
  if options.locipertag:
    return lambda n: min(int(n//options.locipertag)+1,n)
  elif options.loglocipertag:
    l = log(options.loglocipertag)
    return lambda n: int(ceil(log(n+1)/l))
  else:
    return None


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', nargs='+', help='Input genotype file(s)')

  inputgroup = parser.add_argument_group('Input options')

  geno_options(inputgroup,input=True,filter=True)

  inputgroup.add_argument('-e', '--excludetag', dest='exclude', metavar='FILE', default='',
                          help='File containing loci that are excluded from being a tag')
  inputgroup.add_argument('-i', '--includeuntyped', dest='include_untyped', metavar='FILE', default='',
                          help='File containing loci that are obligatorily tags and untyped (may not cover another obligate locus)')
  inputgroup.add_argument('-I', '--includetyped', dest='include_typed', metavar='FILE', default='',
                          help='File containing loci that are obligatorily tags but have been typed (may cover another typed locus)')

  inputgroup.add_argument('-s', '--subset', metavar='FILE', default='',
                          help='File containing loci to be used in analysis')
  inputgroup.add_argument('-S', '--ldsubset', metavar='FILE', default='',
                          help='File containing loci within the region these loci LD will be analyzed (see -d/--maxdist)')
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

  outputgroup.add_argument('-b', '--summary', dest='sumfile', metavar='FILE', default='-',
                          help="Output summary tables FILE (default='-' for standard out)")
  outputgroup.add_argument('-B', '--bininfo', metavar='FILE',
                          help='Output summary information about each bin to FILE')
  outputgroup.add_argument('-H', '--histomax', metavar='N', type=int, default=10,
                          help='Largest bin size output in summary histogram output (default=10)')
  outputgroup.add_argument('-k', '--skip', default=0, action='count',
                          help='Skip output of untagged or excluded loci')
  outputgroup.add_argument('-o', '--output', metavar='FILE', default=None,
                          help="Output tabular LD information for bins to FILE ('-' for standard out)")
  outputgroup.add_argument('-O', '--locusinfo', metavar='FILE',
                          help='Output locus information to FILE')
  outputgroup.add_argument('-u', '--saveldpairs', metavar='FILE',
                          help='Output pairwise LD estimates to FILE')
  outputgroup.add_argument('-x', '--extra', action='count',
                          help='Output inter-bin LD statistics')

  genoldgroup = parser.add_argument_group('Genotype and LD estimation options')

  genoldgroup.add_argument('-a', '--minmaf', dest='maf', metavar='FREQ', type=float, default=0.05,
                          help='Minimum minor allele frequency (MAF) (default=0.05)')
  genoldgroup.add_argument('-A', '--minobmaf', dest='obmaf', metavar='FREQ', type=float, default=None,
                          help='Minimum minor allele frequency (MAF) for obligate tags (defaults to -a/--minmaf)')
  genoldgroup.add_argument('-c', '--mincompletion', metavar='N', default=0, type=int,
                          help='Drop loci with less than N valid genotypes (default=0)')
  genoldgroup.add_argument(      '--mincompletionrate', metavar='N', default=0, type=float,
                          help='Drop loci with completion rate less than N (0-1) (default=0)')
  genoldgroup.add_argument('-m', '--maxdist', metavar='D', type=int, default=200,
                          help='Maximum inter-marker distance in kb for LD comparison (default=200)')
  genoldgroup.add_argument('-P', '--hwp', metavar='p', default=None, type=float,
                          help='Filter out loci that fail to meet a minimum significance level (pvalue) for a '
                               'test Hardy-Weinberg proportion (no default)')

  bingroup = parser.add_argument_group('Binning options')

  bingroup.add_argument('-M', '--multipopulation', metavar='N or P1,P2,...',
                          help='Multipopulation tagging where every N input files represent a group of populations. '
                               'May be specified as an integer N or a comma separated list of population labels.')
  bingroup.add_argument('-d', '--dthreshold', dest='d', metavar='DPRIME', type=float, default=0.,
                          help='Minimum d-prime threshold to output (default=0)')
  bingroup.add_argument('-r', '--rthreshold', dest='r', metavar='N', type=float, default=0.8,
                          help='Minimum r-squared threshold to output (default=0.8)')
  bingroup.add_argument('-t', '--targetbins', metavar='N', type=int, default=0,
                          help='Stop when N bins have been selected (default=0 for unlimited)')
  bingroup.add_argument('-T', '--targetloci', metavar='N', type=int, default=0,
                          help='Stop when N loci have been tagged (default=0 for unlimited)')
  bingroup.add_argument('-C', '--tagcriteria', type=str, metavar='crit', action='append',
                          help='Use the specified criteria to choose the optimal tag for each bin')
  bingroup.add_argument('-z', '--locipertag', metavar='N', type=int, default=None,
                          help='Ensure that bins contain more than one tag per N loci.  Bins with an insufficient number of tags will be reduced.')
  bingroup.add_argument('-Z', '--loglocipertag', metavar='B', type=float, default=None,
                          help='Ensure that bins contains more than one tag per log_B(loci).  Bins with an insufficient number of tags will be reduced.')
  bingroup.add_argument('--skipbinning', action='count',
                          help='Skip binning step.  Typically used in conjunction with -u/--saveldpairs')

  return parser


def do_tagging(ldpairs, locusmap, includes, exclude, designscores, options):
  # If we require a total ordering, then build binsets from all ldpairs
  if options.targetbins or options.targetloci:
    sys.stderr.write('[%s] Building global binsets\n' % time.asctime())
    binsets,lddata = build_binsets(locusmap, ldpairs, includes, exclude, designscores)
    sys.stderr.write('[%s] Choosing global bins\n' % time.asctime())
    bins = binner(locusmap, binsets, lddata, includes, get_tags_required_function(options))
    yield bins,lddata
  else:
    # Otherwise, process each sequence of ldpairs independently
    for i,pairs in enumerate(ldpairs):
      sys.stderr.write('[%s] Generating LD and binsets for region %d\n' % (time.asctime(),i+1))
      binsets,lddata = build_binsets(locusmap, [pairs], includes, exclude, designscores)
      sys.stderr.write('[%s] Choosing bins for region %d\n' % (time.asctime(),i+1))
      bins = binner(locusmap, binsets, lddata, includes, get_tags_required_function(options))
      yield bins,lddata


def build_output(options, exclude):
  pairinfofile = None
  if options.output:
    pairinfofile = table_writer(options.output,hyphen=sys.stdout)
    pairinfo = PairwiseBinOutput(pairinfofile, exclude)
  else:
    pairinfo = NullPairwiseBinOutput()

  locusinfofile = None
  if options.locusinfo:
    locusinfofile = table_writer(options.locusinfo,hyphen=sys.stdout)
    locusinfo = LocusOutput(locusinfofile, exclude)
  else:
    locusinfo = NullLocusOutput()

  infofile = None
  if options.bininfo:
    infofile = autofile(hyphen(options.bininfo,sys.stdout),'w')

  if options.bininfo or options.sumfile:
    bininfo = BinInfo(infofile, options.histomax+1)
  else:
    bininfo = NullBinInfo()

  sumfile = autofile(hyphen(options.sumfile,sys.stdout),'w')

  if [pairinfofile,locusinfofile,infofile,sumfile].count(sys.stdout) > 1:
    raise TagZillaError('ERROR: More than one output file directed to standard out.')

  return pairinfo,locusinfo,bininfo,sumfile


class Includes(object):
  def __init__(self, typed, untyped):
    if typed is None:
      typed = set()
    if untyped is None:
      untyped = set()

    self.typed   = typed - untyped
    self.untyped = untyped

  def __contains__(self, other):
    return other in self.typed or other in self.untyped

  def __iter__(self):
    return chain(self.typed,self.untyped)

  def __len__(self):
    return len(self.typed)+len(self.untyped)


def tagzilla_single(options):
  subset          = None
  include_untyped = None
  include_typed   = None
  exclude         = None
  ldsubset        = None

  if options.subset:
    subset = set(list_reader(options.subset))

  if options.ldsubset:
    ldsubset = set(list_reader(options.ldsubset))

  if options.include_untyped:
    include_untyped = set(list_reader(options.include_untyped))

  if options.include_typed:
    include_typed = set(list_reader(options.include_typed))

  if options.exclude:
    exclude = set(list_reader(options.exclude))
  elif options.designscores:
    exclude = set()

  includes     = Includes(include_typed, include_untyped)
  designscores = build_design_score(options.designscores,options.designdefault)
  tagcriteria  = build_tag_criteria(options.tagcriteria)
  tagselector  = TagSelector(designscores, tagcriteria)

  locusmap = {}
  ldpairs = generate_ldpairs(options, locusmap, includes, subset, ldsubset)

  if options.saveldpairs:
    ldpairs = save_ldpairs(options.saveldpairs, ldpairs)

  if options.skipbinning:
    # Fast trick to save all ld results, but not store them
    for pairs in ldpairs:
      list(dropwhile(lambda x: True, pairs))
    return

  results = do_tagging(ldpairs, locusmap, includes, exclude, designscores, options)

  pairinfo,locusinfo,bininfo,sumfile = build_output(options, exclude)

  binned_loci = 0
  binnum = 0
  tags = set()

  try:
    population = get_populations(options.multipopulation)[0]
  except IndexError:
    pass

  for bins,lddata in results:
    for bin in bins:
      binnum += 1
      bin.binnum = binnum

      qualifier = bin_qualifier(bin, binned_loci, options)

      # Update binned loci after determining qualifier
      binned_loci += len(bin)

      tags.update(bin.tags)

      tagselector.select_tags(bin)

      bininfo.emit_bin(bin, locusmap, exclude, population)
      pairinfo.emit_bin(bin, qualifier, population, options)
      locusinfo.emit_bin(bin, locusmap, qualifier, population)

    # Process remaining items in lddata and output residual ld information
    # (i.e. the inter-bin pairwise)
    if options.extra:
      pairinfo.emit_extra(lddata, tags, population)

  # Emit useful bin summary table
  bininfo.emit_summary(sumfile, population)


def tag_intersection(results):
  ires = (r for r in results if r is not None)
  tags = set()

  tags.update(ires.next())

  for r in ires:
    tags.intersection_update(r.tags)

  return tags


def subset_tags(result, tags, recommended=None):
  diff = set(result.tags) - tags
  result.tags = tags
  result.others.extend(diff)
  if recommended:
    result.recommended_tags = list(recommended.intersection(tags))


def tagzilla_multi(options):
  subset          = None
  ldsubset        = None
  include_untyped = None
  include_typed   = None
  exclude         = None

  if options.subset:
    subset = set(list_reader(options.subset))

  if options.ldsubset:
    ldsubset = set(list_reader(options.ldsubset))

  if options.include_untyped:
    include_untyped = set(list_reader(options.include_untyped))

  if options.include_typed:
    include_typed = set(list_reader(options.include_typed))

  if options.exclude:
    exclude = set(list_reader(options.exclude))
  elif options.designscores:
    exclude = set()

  includes     = Includes(include_typed, include_untyped)
  designscores = build_design_score(options.designscores,options.designdefault)
  tagcriteria  = build_tag_criteria(options.tagcriteria)
  tagselector  = TagSelector(designscores, tagcriteria)

  pairinfo,locusinfo,bininfo,sumfile = build_output(options, exclude)

  ldpairs = generate_ldpairs_vector(options, includes, subset, ldsubset)
  results = do_tagging_vector(ldpairs, includes, exclude, designscores, options)

  labels = get_populations(options.multipopulation)

  binnum = 0
  binned_loci = {}
  poptags = {}
  popdtags = {}

  for resultset,lddata,locusmap in results:
    try:
      tags,resultset = zip(*resultset)
    except ValueError:
      tags,resultset = [],[]

    # Refine tags -- once enabled, also add stags to subset_tags below
    if 0:
      results = zip(*resultset)
      mtags   = set(merge_bins(results))
      stags   = set(shrink_tags(mtags, results))
      print len(tags),len(mtags),len(stags)

    for res in resultset:
      tags = tag_intersection(res)
      # FIXME: This should eventually set the intersection as the bin.tags
      #        and exclude all other candidate tags
      recommended = [iter(tags).next()]

      binnum += 1
      disposition = None
      for population,bin,lmap in izip(labels,res,locusmap):
        if bin is not None:
          bin.binnum = binnum

          # FIXME: This should eventually set the intersection as the bin.tags
          #        and exclude all other candidate tags.
          subset_tags(bin, tags)

          qualifier = bin_qualifier(bin, binned_loci.get(population,0), options)

          # Update binned loci after determining qualifier
          binned_loci[population] = binned_loci.get(population,0) + len(bin)

          poptags.setdefault(population, set()).update(bin.tags)
          disposition = bin.disposition

          # FIXME: The recommended tag must be selected for this method to
          #        ensure across-population coverage.
          #        tagselector.select_tags(bin) must be extended to pick the
          #        recommended among several parallel bins.
          bin.recommended_tags = recommended

          bininfo.emit_bin(bin, lmap, exclude, population)
          pairinfo.emit_bin(bin, qualifier, population, options)
          locusinfo.emit_bin(bin, lmap, qualifier, population)

      popdtags[disposition] = popdtags.get(disposition,0) + 1

    # Process remaining items in lddata and output residual ld information
    # (i.e. the inter-bin pairwise)
    if options.extra:
      for ldd,population in izip(lddata,labels):
        pairinfo.emit_extra(ldd, poptags[population], population)

  # Emit useful bin summary table
  for population in labels:
    bininfo.emit_summary(sumfile, population)

  bininfo.emit_multipop_summary(sumfile, popdtags)


def main():
  parser  = option_parser()
  options = parser.parse_args()

  pops = len(get_populations(options.multipopulation))

  if pops > 1:
    return tagzilla_multi(options)
  else:
    return tagzilla_single(options)


if __name__ == '__main__':
  main()
