# -*- coding: utf-8 -*-
'''
File:          genomerge.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   operator          import itemgetter
from   itertools         import izip,starmap,dropwhile
from   collections       import defaultdict

from   glu.lib.utils     import tally
from   glu.lib.fileutils import autofile,hyphen


UNAMBIGUOUS,CONCORDANT,CONSENSUS,DISCORDANT,MISSING = range(5)
CONCORDANCE_HEADER = ['UNAMBIGUOUS_COUNT','CONCORDANT_COUNT','CONSENSUS_COUNT','DISCORDANT_COUNT',
                      'DISCONCORDANCE_RATE','MISSING_COUNT','CONSENSUS_MISSING_RATE']
SAMPLE_CONCORDANCE_HEADER = ['SAMPLE_ID']+ CONCORDANCE_HEADER
LOCUS_CONCORDANCE_HEADER  = ['LOCUS_ID'] + CONCORDANCE_HEADER


class NonUniqueGenotypeError(ValueError): pass


class Unique(object):
  '''
  Require unique genotypes, ie exactly zero or one genotypes to be
  specified, when forming a consensus genotype.  An error is generated when
  two or more genotypes are specified.  This merge method should be used
  whenever input data are required to be unique, as this will be enforced
  rigorously.

  >>> unique = Unique()
  >>> unique(None)
  (4, (None, None))
  >>> unique(['AA'])
  (0, 'AA')
  >>> unique([None])
  (4, (None, None))
  >>> unique([None,None])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(['AA','AA'])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(['AA','AB'])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(['AA',None])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  '''
  def __init__(self):
    '''
    '''

  def __call__(self,genos):
    '''
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''
    if not genos:
      return MISSING,(None, None)
    elif len(genos) == 1:
      geno = genos[0]
      if not geno:
        return MISSING,(None,None)
      else:
        return UNAMBIGUOUS,geno
    else:
      raise NonUniqueGenotypeError, 'Non-unique genotype found'


class Vote(object):
  '''
  Determine a consensus genotype and a concordance class based on how that
  genotype was chosen.  A voting process is introduced to pick a consensus
  genotype if two or more non-missing genotypes exist.  Only the genotype
  with better percentage than the specified threshold will be returned,
  otherwise a missing genotype is returned.

  Classes:
   0) unambiguous:  exactly one genotype and non-missing
   1) concordant:   two or more concordant non-missing genotypes
   2) consensus:    consensus of two or more non-missing genotypes determined by voting
   3) discordant:   two or more non-missing genotypes that do not meet voting threshold
   4) missing:      zero or more missing genotypes only

  >>> vote = Vote(1.0)
  >>> vote(None)
  (4, (None, None))
  >>> vote(['AA'])
  (0, 'AA')
  >>> vote([1])
  (0, 1)
  >>> vote([1,1])
  (1, 1)
  >>> vote([None,None,None])
  (4, (None, None))
  >>> vote(['AA',None,None,None])
  (1, 'AA')
  >>> vote(['AA','AA',None,'AA','AA','AA'])
  (1, 'AA')
  >>> vote(['AA','AA',None,'AA','AB','AA'])
  (3, (None, None))
  >>> vote = Vote(0.8)
  >>> vote(['AA','AA',None,'AA','AB','AA'])
  (2, 'AA')
  >>> vote = Vote(0.4)
  >>> vote(['AT','AT','AA','TT','AC',None,None])
  (2, 'AT')
  >>> vote(['AA','AA','AA','BB','BB','BB'])
  (3, (None, None))
  >>> vote(['AA','AA','BB','BB','AB'])
  (3, (None, None))
  '''
  def __init__(self, threshold=1.0):
    '''
    @param   threshold: cut-off value to be voted as a consensus
    @type    threshold: float
    '''
    self.threshold   = threshold

  def __call__(self,genos):
    '''
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''
    # Fast path
    if not genos:
      return MISSING,(None, None)
    elif len(genos)==1:
      geno = genos[0]
      if not geno:
        return MISSING,(None,None)
      else:
        return UNAMBIGUOUS,geno

    # Slow path
    genocounts = tally(g for g in genos if g)

    if not genocounts:
      return MISSING,(None, None)
    elif len(genocounts) == 1:
      return CONCORDANT,iter(genocounts).next()

    counts = sorted(genocounts.iteritems(), key=itemgetter(1))

    winner   = counts[-1][1]
    runnerup = counts[-2][1]

    total  = sum(genocounts.itervalues())
    votes  = float(winner)/total

    # If the majority is sufficient and the winner has more votes then the
    # runner-up let it win.
    if votes >= self.threshold and winner > runnerup:
      return CONSENSUS,counts[-1][0]
    else:
      return DISCORDANT,(None, None)


class Ordered(object):
  '''
  Determine a consensus genotype and a concordance class based on how that
  genotype was chosen.  When ambiguous, the first observed genotype will be
  chosen, provided it meets a minimum concordance threshold (default=49.9%).

  A voting process is introduced to pick a consensus genotype if two or more
  non-missing genotypes exist.  Only the first non-missing genotype with
  better frequency higher the specified threshold will be returned,
  otherwise None.

  Classes:
  0) unambiguous:  exactly one genotype and non-missing
  1) concordant:   two or more concordant non-missing genotypes
  2) consensus:    consensus of two or more non-missing genotypes determined by voting
  3) discordant:   two or more non-missing genotypes that do not meet voting threshold
  4) missing:      zero or more missing genotypes only

  >>> ordered = Ordered(0)
  >>> ordered(None)
  (4, (None, None))
  >>> ordered(['AA'])
  (0, 'AA')
  >>> ordered([1])
  (0, 1)
  >>> ordered([1,1])
  (1, 1)
  >>> ordered([None,None,None])
  (4, (None, None))
  >>> ordered(['AA',None,None,None])
  (1, 'AA')
  >>> ordered(['AA','AA',None,'AA','AA','AA'])
  (1, 'AA')
  >>> ordered(['AA','AA',None,'AA','AB','AA'])
  (2, 'AA')
  >>> ordered(['AA','AA',None,'AA','AB','AA'])
  (2, 'AA')
  >>> ordered(['AT','AT','AA','TT','AC',None,None])
  (2, 'AT')
  >>> ordered(['AA','AA','BB','BB','BB'])
  (2, 'AA')
  >>> ordered(['AA','AA','BB','BB','AB'])
  (2, 'AA')
  >>> ordered = Ordered()
  >>> ordered(None)
  (4, (None, None))
  >>> ordered(['AA'])
  (0, 'AA')
  >>> ordered([1])
  (0, 1)
  >>> ordered([1,1])
  (1, 1)
  >>> ordered([None,None,None])
  (4, (None, None))
  >>> ordered(['AA',None,None,None])
  (1, 'AA')
  >>> ordered(['AA','AA',None,'AA','AA','AA'])
  (1, 'AA')
  >>> ordered(['AA','AA',None,'AA','AB','AA'])
  (2, 'AA')
  >>> ordered(['AA','AA',None,'AA','AB','AA'])
  (2, 'AA')
  >>> ordered(['AT','AT','AA','TT','AC',None,None])
  (3, (None, None))
  >>> ordered(['AA','AA','AA','BB','BB','BB'])
  (2, 'AA')
  >>> ordered(['AA','AA','BB','BB','AB'])
  (3, (None, None))
  '''

  def __init__(self, threshold=0.4999999):
    '''
    @param   threshold: cut-off value to be voted as a consensus
    @type    threshold: float
    '''
    self.threshold   = threshold

  def __call__(self,genos):
    '''
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''
    # Fast path
    if not genos:
      return MISSING,(None, None)
    elif len(genos)==1:
      geno = genos[0]
      if not geno:
        return MISSING,(None,None)
      else:
        return UNAMBIGUOUS,geno

    # Slow path
    genocounts = tally(g for g in genos if g)

    if not genocounts:
      return MISSING,(None, None)
    elif len(genocounts) == 1:
      return CONCORDANT,iter(genocounts).next()

    total = sum(genocounts.itervalues())
    geno  = dropwhile(lambda g: not g, genos).next()

    votes  = float(genocounts[geno])/total

    if votes >= self.threshold:
      return CONSENSUS,geno
    else:
      return DISCORDANT,(None, None)


class Merger(object):
  '''
  Object to form consensus/merged genotypes and track statistics on genotypes that
  have been processed.

  Two statistics objects are maintained, samplestats and locustats.  They
  are dictionaries from sample and locus, respectively, to a three element
  list contianing the following genotype counts:

  0) unambiguous:  exactly one genotype and non-missing
  1) concordant:   two or more concordant non-missing gentypes
  2) consensus:    consensus of two or more non-missing genotypes determined by voting
  3) discordant:   two or more non-missing genotypes that do not meet voting threshold
  4) missing:      zero or more missing genotypes only
  '''
  def __init__(self, votefunc):
    '''
    Create a new Merger object
    '''
    self.votefunc = votefunc
    self.clear()

  def clear(self):
    '''
    Clear concordance statistics
    '''
    zeros = lambda: [0,0,0,0,0]
    self.samplestats = defaultdict(zeros)
    self.locusstats  = defaultdict(zeros)

  def __call__(self,sample,locus,genos):
    '''
    Merge a group of genotypes into a consensus for a given sample and locus
    and build merge statistics

    @param sample: sample id
    @type  sample: str
    @param  locus: locus id
    @type   locus: str
    @param  genos: sequence of genotypes
    @type   genos: sequence of objects
    @return      : consensus genotype
    @rtype       : object

    >>> samples  = ['s1','s1','s2','s2','s3','s3']
    >>> loci     = ['l1','l2','l3','l1','l2','l3']
    >>> genosets = [['AA','AA',None,'AA','AB','AA'],None,['AA','AB','AB'],[None,None],[None],['AA']]
    >>> merger = Merger(Vote(threshold = 0.8))
    >>> list(starmap(merger,izip(samples,loci,genosets)))
    ['AA', (None, None), (None, None), (None, None), (None, None), 'AA']
    >>> sorted(merger.samplestats.iteritems())
    [('s1', [0, 0, 1, 0, 1]), ('s2', [0, 0, 0, 1, 1]), ('s3', [1, 0, 0, 0, 1])]
    >>> sorted(merger.locusstats.iteritems())
    [('l1', [0, 0, 1, 0, 1]), ('l2', [0, 0, 0, 0, 2]), ('l3', [1, 0, 0, 1, 0])]
    >>> merger = Merger(Vote(threshold = 0.8))
    >>> list(starmap(mergefunc_transpose_adapter(merger),izip(loci,samples,genosets)))
    ['AA', (None, None), (None, None), (None, None), (None, None), 'AA']
    >>> sorted(merger.samplestats.iteritems())
    [('s1', [0, 0, 1, 0, 1]), ('s2', [0, 0, 0, 1, 1]), ('s3', [1, 0, 0, 0, 1])]
    >>> sorted(merger.locusstats.iteritems())
    [('l1', [0, 0, 1, 0, 1]), ('l2', [0, 0, 0, 0, 2]), ('l3', [1, 0, 0, 1, 0])]
    '''
    concordclass,geno = self.votefunc(genos)

    self.samplestats[sample][concordclass] += 1
    self.locusstats[locus][concordclass]   += 1

    return geno


def UniqueMerger(threshold=None):
  return Merger(Unique())


def VoteMerger(threshold=1.0):
  return Merger(Vote(threshold=threshold))


def OrderedMerger(threshold=0.4999999):
  return Merger(Ordered(threshold=threshold))


def get_genomerger(reprname):
  '''
  Retrieve the supported genotype merge algorithm. Otherwise raises an ValueError exception.

  @param reprname: genotype merger name
  @type  reprname: str
  @return        : genotype merge object
  '''
  parts = reprname.split(':')

  if not parts:
    raise ValueError, 'Missing merge algorithm specified'

  name = parts[0].lower()

  if name == 'unique':
    merger = UniqueMerger
  elif name == 'vote':
    merger = VoteMerger
  elif name == 'ordered':
    merger = OrderedMerger
  else:
    raise ValueError, 'Unknown merge algorithm specified'

  if len(parts) == 2:
    threshold=float(parts[1])
    return merger(threshold=threshold)
  else:
    return merger()


def output_merge_statistics(mergefunc,samplefile=None,locusfile=None):
  '''
  Retrieve merge statistics from the object and write to files by sample/locus

  >>> import StringIO
  >>> samples  = ['s1','s1','s1','s1','s2','s2','s2','s2','s3','s3','s3','s3','s4','s4','s4','s4']
  >>> loci     = ['l1','l2','l3','l4','l1','l2','l3','l4','l1','l2','l3','l4','l1','l2','l3','l4']
  >>> genosets = [['AA','AA',None,'AA','AB','AA'],['AA','AB','AB'],['AA','AA','AA'],[None,None],
  ...            [None],['AA'],['AA',None,None,None],['AT','AT','AA','TT','AC',None,None],
  ...            ['AB'],['AA'],['AA','AA','AA','AA','AB'],['AB','AB'],
  ...            ['AA','AB','AB'],['AA','AA'],['AA',None,None,None],['AA','AA','AA']]
  >>> merger = Merger(Vote(threshold = 0.8))
  >>> list(starmap(merger, izip(samples,loci,genosets)))
  ['AA', (None, None), 'AA', (None, None), (None, None), 'AA', 'AA', (None, None), 'AB', 'AA', 'AA', 'AB', (None, None), 'AA', 'AA', 'AA']
  >>> sampleout = StringIO.StringIO()
  >>> locusout  = StringIO.StringIO()
  >>> output_merge_statistics(merger,sampleout,locusout)
  >>> print sampleout.getvalue() # doctest: +NORMALIZE_WHITESPACE
  SAMPLE_ID   UNAMBIGUOUS_COUNT       CONCORDANT_COUNT        CONSENSUS_COUNT DISCORDANT_COUNT        DISCONCORDANCE_RATE       MISSING_COUNT    CONSENSUS_MISSING_RATE
  s3  2       1       1       0       0.000000        0       0.000000
  s4  0       3       0       1       0.250000        0       0.250000
  s1  0       1       1       1       0.333333        1       0.500000
  s2  1       1       0       1       0.500000        1       0.500000
  <BLANKLINE>
  >>> print locusout.getvalue() # doctest: +NORMALIZE_WHITESPACE
  LOCUS_ID    UNAMBIGUOUS_COUNT       CONCORDANT_COUNT        CONSENSUS_COUNT DISCORDANT_COUNT        DISCONCORDANCE_RATE       MISSING_COUNT    CONSENSUS_MISSING_RATE
  l3  0       3       1       0       0.000000        0       0.000000
  l4  0       2       0       1       0.333333        1       0.500000
  l2  2       1       0       1       0.500000        0       0.250000
  l1  1       0       1       1       0.500000        1       0.500000
  <BLANKLINE>
  '''
  if samplefile is not None:
    f = csv.writer(autofile(samplefile,'w'), dialect='tsv')
    f.writerow(SAMPLE_CONCORDANCE_HEADER)
    f.writerows(sorted(build_concordance_output(mergefunc.samplestats),key=itemgetter(5)))

  if locusfile is not None:
    fl = csv.writer(autofile(locusfile,'w'), dialect='tsv')
    fl.writerow(LOCUS_CONCORDANCE_HEADER)
    fl.writerows(sorted(build_concordance_output(mergefunc.locusstats),key=itemgetter(5)))


def build_concordance_output(stats):
  '''
  Calculate concordance rate and format all the results for output
  '''
  for key,(unambiguous,concordant,consensus,discordant,missing) in stats.iteritems():
    total = unambiguous+concordant+discordant+consensus+missing
    concord_total = concordant+consensus
    # discordance rate for non-missing genotypes only
    discord_rate = 0
    if concord_total or discordant:
      discord_rate  = '%.6f' % (float(discordant)/(concord_total+discordant))
    # calculate the percentage of missing genotypes after mergerd
    missing_rate  = '%.6f' % (float(missing+discordant)/total)
    yield key,unambiguous,concordant,consensus,discordant,discord_rate,missing,missing_rate


def mergefunc_transpose_adapter(mergefunc):
  '''
  Return a new mergefunc that transposes the calling convention for samples and loci
  '''
  def _mergefunc_transpose_adapter(locus,sample,genos):
    return mergefunc(sample,locus,genos)
  return _mergefunc_transpose_adapter


def _test():
  import doctest
  return doctest.testmod()


if __name__=='__main__':
  _test()