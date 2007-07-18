from hwp import hwp_biallelic

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
    # For mafs[locus.name in include], the index evaluatates to:
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
      raise TagZillaError,'ERROR: Invalid genomic range: %s' % range

    ranges.append( (start,stop) )

  if range_all in ranges:
    ranges = [range_all]

  for locus in loci:
    for start,stop in ranges:
      if start <= locus.location < stop:
        yield locus
        break


def completion(locus):
  return len(locus) - locus.count('  '),len(locus)


def filter_loci_by_completion(loci, mincompletion, mincompletionrate):
  '''Generator that filters loci by a minimum completion rate'''

  for locus in loci:
    m,n = completion(locus.genos)

    rate = 0
    if n:
      rate = float(m)/n

    if m >= mincompletion and rate >= mincompletionrate:
      yield locus
