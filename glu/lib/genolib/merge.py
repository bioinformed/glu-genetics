# -*- coding: utf-8 -*-

__abstract__  = 'genotype merge algorithms and infrastructure'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   operator                  import itemgetter
from   itertools                 import izip,starmap,dropwhile,repeat,imap
from   collections               import defaultdict

from   glu.lib.utils             import tally
from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib.genoarray import Genotype,model_from_alleles


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

  >>> model = model_from_alleles('AB')
  >>> missing = model[None,None]
  >>> A,B='A','B'
  >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

  >>> unique = Unique()
  >>> unique(model,[])
  (4, (None, None))
  >>> unique(model,[AA])
  (0, ('A', 'A'))
  >>> unique(model,[missing])
  (4, (None, None))
  >>> unique(model,[missing,missing])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(model,[AA,AA,AA])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(model,[AA,AB])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  >>> unique(model,[AA,missing])
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found
  '''
  def __init__(self):
    '''
    '''

  def __call__(self,model,genos):
    '''
    @param       model: new internal representation of genotypes
    @type        model: UnphasedMarkerRepresentation or similar object
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''
    if not genos:
      return MISSING,model[None,None]
    elif len(genos) == 1:
      geno = genos[0]
      if not geno:
        return MISSING,geno
      else:
        return UNAMBIGUOUS,geno
    else:
      raise NonUniqueGenotypeError('Non-unique genotype found')

try:
  from _genoarray import merge_unanimous as unanimous

  def test_unanimous():
    '''
    >>> model = model_from_alleles('AB')
    >>> missing = model[None,None]
    >>> A,B='A','B'
    >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

    >>> unanimous(model,missing)
    (4, (None, None))
    >>> unanimous(model,None)
    (4, (None, None))
    >>> unanimous(model,[])
    (4, (None, None))
    >>> unanimous(model,AA)
    (0, ('A', 'A'))
    >>> unanimous(model,[AA])
    (0, ('A', 'A'))
    >>> unanimous(model,[missing,missing,missing])
    (4, (None, None))
    >>> unanimous(model,[AA,missing,missing,missing])
    (0, ('A', 'A'))
    >>> unanimous(model,[AA,AA,missing,AA,AA,AA])
    (1, ('A', 'A'))
    >>> unanimous(model,[AA,AA,missing,AA,AB,AA])
    (3, (None, None))
    '''

  def test_genotype_merger():
    '''
    >>> model = model_from_alleles('AB')
    >>> missing = model[None,None]
    >>> A,B='A','B'
    >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

    >>> merger = GenotypeMerger(unanimous)
    >>> merger.merge_geno('s1','l1',model,None)
    (None, None)
    >>> merger.merge_geno('s1','l1',model,[])
    (None, None)
    >>> merger.merge_geno('s1','l1',model,missing)
    (None, None)
    >>> merger.merge_geno('s1','l1',model,AA)
    ('A', 'A')
    >>> merger.merge_geno('s1','l1',model,[missing])
    (None, None)
    >>> merger.merge_geno('s1','l1',model,[AA])
    ('A', 'A')
    >>> merger.merge_geno('s1','l1',model,[AA,missing])
    ('A', 'A')
    >>> merger.merge_geno('s1','l1',model,[AA,missing,AA])
    ('A', 'A')
    >>> merger.merge_geno('s1','l1',model,[AA,AA])
    ('A', 'A')
    >>> merger.merge_geno('s1','l1',model,[AA,AA,AB])
    (None, None)

    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [3, 2, 0, 1, 4])]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [3, 2, 0, 1, 4])]

    >>> samples = ['s1','s2','s3','s4','s5','s6']
    >>> locus   = 'l1'
    >>> row     = [[AA,AA,missing,AA,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

    >>> merger = GenotypeMerger(unanimous)
    >>> merger.merge_locus(samples, locus, model, row)
    [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]

    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [1, 1, 0, 1, 3])]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [0, 1, 0, 0, 0]), ('s2', [0, 0, 0, 0, 1]), ('s3', [0, 0, 0, 1, 0]), ('s4', [0, 0, 0, 0, 1]), ('s5', [0, 0, 0, 0, 1]), ('s6', [1, 0, 0, 0, 0])]

    >>> sample = 's1'
    >>> loci   = ['l1','l2','l3','l4','l5','l6']
    >>> row    = [[AA,AA,missing,AA,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

    >>> merger = GenotypeMerger(unanimous)
    >>> merger.merge_sample(sample, loci, [model]*6, row)
    [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]

    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [0, 1, 0, 0, 0]), ('l2', [0, 0, 0, 0, 1]), ('l3', [0, 0, 0, 1, 0]), ('l4', [0, 0, 0, 0, 1]), ('l5', [0, 0, 0, 0, 1]), ('l6', [1, 0, 0, 0, 0])]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [1, 1, 0, 1, 3])]
    '''

except ImportError:
  def unanimous(model,genos):
    '''
    Determine a consensus genotype and a concordance class based on how that
    genotype was chosen.  This method requires unanimous selection of
    non-missing genotypes.

    Classes:
     0) unambiguous:  exactly one genotype and non-missing
     1) concordant:   two or more concordant non-missing genotypes
     2) consensus:    consensus of two or more non-missing genotypes determined by voting
     3) discordant:   two or more non-missing genotypes that do not meet voting threshold
     4) missing:      zero or more missing genotypes only

    @param       model: new internal representation of genotypes
    @type        model: UnphasedMarkerRepresentation or similar object
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str

    >>> model = model_from_alleles('AB')
    >>> missing = model[None,None]
    >>> A,B='A','B'
    >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

    >>> unanimous(model,missing)
    (4, (None, None))
    >>> unanimous(model,[AA])
    (0, ('A', 'A'))
    >>> unanimous(model,[missing,missing,missing])
    (4, (None, None))
    >>> unanimous(model,[AA,missing,missing,missing])
    (0, ('A', 'A'))
    >>> unanimous(model,[AA,AA,missing,AA,AA,AA])
    (1, ('A', 'A'))
    >>> unanimous(model,[AA,AA,missing,AA,AB,AA])
    (3, (None, None))
    '''
    # Fast path
    if not genos:
      return MISSING,model[None,None]
    elif isinstance(genos, Genotype):
      return UNAMBIGUOUS,genos
    elif len(genos)==1:
      geno = genos[0]
      if not geno:
        return MISSING,geno
      else:
        return UNAMBIGUOUS,geno

    # Slow path
    genocounts = tally(g for g in genos if g)

    if not genocounts:
      return MISSING,model[None,None]
    elif len(genocounts) == 1:
      geno,n = genocounts.popitem()
      if n == 1:
        return UNAMBIGUOUS,geno
      else:
        return CONCORDANT,geno
    else:
      return DISCORDANT,model[None,None]


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

  >>> model = model_from_alleles('AB')
  >>> missing = model[None,None]
  >>> A,B='A','B'
  >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

  >>> vote = Vote(1.0)
  >>> vote(model,missing)
  (4, (None, None))
  >>> vote(model,[AA])
  (0, ('A', 'A'))
  >>> vote(model,[missing,missing,missing])
  (4, (None, None))
  >>> vote(model,[AA,missing,missing,missing])
  (0, ('A', 'A'))
  >>> vote(model,[AA,AA,missing,AA,AA,AA])
  (1, ('A', 'A'))
  >>> vote(model,[AA,AA,missing,AA,AB,AA])
  (3, (None, None))
  >>> vote = Vote(0.8)
  >>> vote(model,[AA,AA,missing,AA,AB,AA])
  (2, ('A', 'A'))
  >>> vote = Vote(0.4)
  >>> vote(model,[AB,AB,AA,BB,missing,missing])
  (2, ('A', 'B'))
  >>> vote(model,[AA,AA,AA,BB,BB,BB])
  (3, (None, None))
  >>> vote(model,[AA,AA,BB,BB,AB])
  (3, (None, None))
  '''
  def __init__(self, threshold=1.0):
    '''
    @param   threshold: cut-off value to be voted as a consensus. Default is 1.0
    @type    threshold: float
    '''
    self.threshold = threshold

  def __call__(self,model,genos):
    '''
    @param       model: new internal representation of genotypes
    @type        model: UnphasedMarkerRepresentation or similar object
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''
    # Fast path
    if not genos:
      return MISSING,model[None,None]
    elif isinstance(genos, Genotype):
      return UNAMBIGUOUS,genos
    elif len(genos)==1:
      geno = genos[0]
      if not geno:
        return MISSING,geno
      else:
        return UNAMBIGUOUS,geno

    # Slow path
    genocounts = tally(g for g in genos if g)

    if not genocounts:
      return MISSING,model[None,None]
    elif len(genocounts) == 1:
      geno,n = genocounts.popitem()
      if n == 1:
        return UNAMBIGUOUS,geno
      else:
        return CONCORDANT,geno

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
      return DISCORDANT,model[None,None]


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

  >>> model = model_from_alleles('AB')
  >>> missing = model[None,None]
  >>> A,B='A','B'
  >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

  >>> ordered = Ordered(0)
  >>> ordered(model,missing)
  (4, (None, None))
  >>> ordered(model,[AA])
  (0, ('A', 'A'))
  >>> ordered(model,[missing,missing,missing])
  (4, (None, None))
  >>> ordered(model,[AA,missing,missing,missing])
  (0, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AA,AA])
  (1, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AB,AA])
  (2, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AB,AA])
  (2, ('A', 'A'))
  >>> ordered(model,[AB,AB,AA,BB,AB,missing,missing])
  (2, ('A', 'B'))
  >>> ordered(model,[AA,AA,BB,BB,BB])
  (2, ('A', 'A'))
  >>> ordered(model,[AA,AA,BB,BB,AB])
  (2, ('A', 'A'))

  >>> ordered = Ordered()
  >>> ordered(model,missing)
  (4, (None, None))
  >>> ordered(model,[AA])
  (0, ('A', 'A'))
  >>> ordered(model,[missing,missing,missing])
  (4, (None, None))
  >>> ordered(model,[AA,missing,missing,missing])
  (0, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AA,AA])
  (1, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AB,AA])
  (2, ('A', 'A'))
  >>> ordered(model,[AA,AA,missing,AA,AB,AA])
  (2, ('A', 'A'))
  >>> ordered(model,[AB,AB,AA,BB,AA,missing,missing])
  (3, (None, None))
  >>> ordered(model,[AA,AA,AA,BB,BB,BB])
  (2, ('A', 'A'))
  >>> ordered(model,[AA,AA,BB,BB,AB])
  (3, (None, None))
  '''
  def __init__(self, threshold=0.4999999):
    '''
    @param   threshold: cut-off value to be voted as a consensus. Default is 0.4999999
    @type    threshold: float
    '''
    self.threshold   = threshold

  def __call__(self,model,genos):
    '''
    @param       model: new internal representation of genotypes
    @type        model: UnphasedMarkerRepresentation or similar object
    @param       genos: genotypes
    @type        genos: sequence
    @return           : concordance class, the consensus genotype
    @rtype            : int,str
    '''

    # Fast path
    if not genos:
      return MISSING,model[None,None]
    elif isinstance(genos, Genotype):
      return UNAMBIGUOUS,genos
    elif len(genos)==1:
      geno = genos[0]
      if not geno:
        return MISSING,geno
      else:
        return UNAMBIGUOUS,geno

    # Slow path
    genocounts = tally(g for g in genos if g)

    if not genocounts:
      return MISSING,model[None,None]
    elif len(genocounts) == 1:
      geno,n = genocounts.popitem()
      if n == 1:
        return UNAMBIGUOUS,geno
      else:
        return CONCORDANT,geno

    total = sum(genocounts.itervalues())
    geno  = dropwhile(lambda g: not g, genos).next()

    votes  = float(genocounts[geno])/total

    if votes >= self.threshold:
      return CONSENSUS,geno
    else:
      return DISCORDANT,model[None,None]


try:
  from _genoarray import GenotypeMerger

  def test_merger():
    '''
    >>> model = model_from_alleles('AB')
    >>> missing = model[None,None]
    >>> A,B='A','B'
    >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

    >>> samples  = ['s1','s1','s2','s2','s3','s3']
    >>> loci     = ['l1','l2','l3','l1','l2','l3']
    >>> genosets = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

    >>> merger = GenotypeMerger(Vote(threshold = 0.8))
    >>> list(starmap(merger.merge_geno,izip(samples,loci,repeat(model),genosets)))
    [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [0, 0, 1, 0, 1]), ('s2', [0, 0, 0, 1, 1]), ('s3', [1, 0, 0, 0, 1])]
    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [0, 0, 1, 0, 1]), ('l2', [0, 0, 0, 0, 2]), ('l3', [1, 0, 0, 1, 0])]

    >>> samples = ['s1','s2','s3','s4','s5','s6']
    >>> locus   = 'l1'
    >>> row     = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

    >>> merger = GenotypeMerger(Vote(threshold = 0.8))
    >>> merger.merge_locus(samples, locus, model, row)
    [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [0, 0, 1, 0, 0]), ('s2', [0, 0, 0, 0, 1]), ('s3', [0, 0, 0, 1, 0]), ('s4', [0, 0, 0, 0, 1]), ('s5', [0, 0, 0, 0, 1]), ('s6', [1, 0, 0, 0, 0])]
    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [1, 0, 1, 1, 3])]

    >>> sample = 's1'
    >>> loci   = ['l1','l2','l3','l4','l5','l6']
    >>> row    = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

    >>> merger = GenotypeMerger(Vote(threshold = 0.8))
    >>> merger.merge_sample(sample, loci, [model]*6, row)
    [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
    >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
    [('s1', [1, 0, 1, 1, 3])]
    >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
    [('l1', [0, 0, 1, 0, 0]), ('l2', [0, 0, 0, 0, 1]), ('l3', [0, 0, 0, 1, 0]), ('l4', [0, 0, 0, 0, 1]), ('l5', [0, 0, 0, 0, 1]), ('l6', [1, 0, 0, 0, 0])]
    '''

except ImportError:
  class GenotypeMerger(object):
    '''
    Object to form consensus/merged genotypes and track statistics on genotypes that
    have been processed.

    Two statistics objects are maintained, samplestats and locustats.  They
    are dictionaries from sample and locus, respectively, to a five element
    list contianing the following genotype counts:

    0) unambiguous:  exactly one genotype and non-missing
    1) concordant:   two or more concordant non-missing gentypes
    2) consensus:    consensus of two or more non-missing genotypes determined by voting
    3) discordant:   two or more non-missing genotypes that do not meet voting threshold
    4) missing:      zero or more missing genotypes only
    '''
    def __init__(self, votefunc):
      '''
      Create a new GenotypeMerger object
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

    def merge_geno(self,sample,locus,model,genos):
      '''
      Merge a group of genotypes into a consensus for a given sample and locus
      and build merge statistics

      @param sample: sample identifier
      @type  sample: str
      @param  locus: locus identifier
      @type   locus: str
      @param  model: new internal representation of genotypes
      @type   model: UnphasedMarkerRepresentation or similar object
      @param  genos: genotypes
      @type   genos: sequence
      @return      : consensus genotype
      @rtype       : object

      >>> model = model_from_alleles('AB')
      >>> missing = model[None,None]
      >>> A,B='A','B'
      >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

      >>> samples  = ['s1','s1','s2','s2','s3','s3']
      >>> loci     = ['l1','l2','l3','l1','l2','l3']
      >>> genosets = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

      >>> merger = GenotypeMerger(Vote(threshold = 0.8))
      >>> list(starmap(merger.merge_geno,izip(samples,loci,repeat(model),genosets)))
      [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
      >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
      [('s1', [0, 0, 1, 0, 1]), ('s2', [0, 0, 0, 1, 1]), ('s3', [1, 0, 0, 0, 1])]
      >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
      [('l1', [0, 0, 1, 0, 1]), ('l2', [0, 0, 0, 0, 2]), ('l3', [1, 0, 0, 1, 0])]
      '''
      concordclass,geno = self.votefunc(model,genos)

      self.samplestats[sample][concordclass] += 1
      self.locusstats[locus][concordclass]   += 1

      return geno

    def merge_locus(self, samples, locus, model, row):
      '''
      Merge a sequence of group of genotypes into a consensus for a list of
      samples at a given locus and build merge statistics

      @param samples: sequence of sample identifiers
      @type  samples: list of str
      @param   locus: locus identifier
      @type    locus: str
      @param   model: genotype model
      @type    model: UnphasedMarkerRepresentation or similar object
      @param     row: genotypes to merge
      @type      row: sequence of genotype sequences
      @return       : consensus genotypes
      @rtype        : list of genotypes

      >>> model = model_from_alleles('AB')
      >>> missing = model[None,None]
      >>> A,B='A','B'
      >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

      >>> samples = ['s1','s2','s3','s4','s5','s6']
      >>> locus   = 'l1'
      >>> row     = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

      >>> merger = GenotypeMerger(Vote(threshold = 0.8))
      >>> merger.merge_locus(samples, locus, model, row)
      [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
      >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
      [('s1', [0, 0, 1, 0, 0]), ('s2', [0, 0, 0, 0, 1]), ('s3', [0, 0, 0, 1, 0]), ('s4', [0, 0, 0, 0, 1]), ('s5', [0, 0, 0, 0, 1]), ('s6', [1, 0, 0, 0, 0])]
      >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
      [('l1', [1, 0, 1, 1, 3])]
      '''
      status,new_row = zip(*imap(self.votefunc,repeat(model),row))

      lstat = self.locusstats[locus]
      sstat = self.samplestats
      for sample,s in izip(samples,status):
        lstat[s] += 1
        sstat[sample][s] += 1

      return list(new_row)

    def merge_sample(self, sample, loci, models, row):
      '''
      Merge a sequence of group of genotypes into a consensus for a list of
      loci for a given sample and build merge statistics

      @param samples: sequence of sample identifiers
      @type  samples: list of str
      @param   locus: locus identifier
      @type    locus: str
      @param   model: genotype model
      @type    model: UnphasedMarkerRepresentation or similar object
      @param     row: genotypes to merge
      @type      row: sequence of genotype sequences
      @return       : consensus genotypes
      @rtype        : list of genotypes

      >>> model = model_from_alleles('AB')
      >>> missing = model[None,None]
      >>> A,B='A','B'
      >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

      >>> sample = 's1'
      >>> loci   = ['l1','l2','l3','l4','l5','l6']
      >>> row    = [[AA,AA,missing,AA,AB,AA],missing,[AA,AB,AB],[missing,missing],[missing],[AA]]

      >>> merger = GenotypeMerger(Vote(threshold = 0.8))
      >>> merger.merge_sample(sample, loci, [model]*6, row)
      [('A', 'A'), (None, None), (None, None), (None, None), (None, None), ('A', 'A')]
      >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
      [('s1', [1, 0, 1, 1, 3])]
      >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
      [('l1', [0, 0, 1, 0, 0]), ('l2', [0, 0, 0, 0, 1]), ('l3', [0, 0, 0, 1, 0]), ('l4', [0, 0, 0, 0, 1]), ('l5', [0, 0, 0, 0, 1]), ('l6', [1, 0, 0, 0, 0])]
      '''
      status,new_row = zip(*imap(self.votefunc,models,row))

      lstat = self.locusstats
      sstat = self.samplestats[sample]
      for locus,s in izip(loci,status):
        lstat[locus][s] += 1
        sstat[s] += 1

      return list(new_row)


def UniqueMerger(threshold=None,trackstats=True):
  return GenotypeMerger(Unique(),trackstats=trackstats)


def UnanimousMerger(threshold=None,trackstats=True):
  return GenotypeMerger(unanimous,trackstats=trackstats)


def VoteMerger(threshold=1.0,trackstats=True):
  if threshold>=1.0:
    return UnanimousMerger(trackstats=trackstats)
  return GenotypeMerger(Vote(threshold=threshold),trackstats=trackstats)


def OrderedMerger(threshold=0.4999999,trackstats=True):
  return GenotypeMerger(Ordered(threshold=threshold),trackstats=trackstats)


def get_genomerger(mergername,trackstats=True):
  '''
  Retrieve the supported genotype merge algorithm. Otherwise raises an ValueError exception.

  @param mergername: genotype merger name
  @type  mergername: str
  @return          : genotype merge object
  @rtype           : merger object
  '''
  parts = mergername.split(':')

  if not parts:
    raise ValueError('Missing merge algorithm specified')

  name = parts[0].lower()

  if name == 'unanimous':
    merger = UnanimousMerger
  elif name == 'unique':
    merger = UniqueMerger
  elif name == 'vote':
    merger = VoteMerger
  elif name == 'ordered':
    merger = OrderedMerger
  elif name == 'none':
    return None
  else:
    raise ValueError('Unknown merge algorithm specified')

  if len(parts) == 2:
    threshold=float(parts[1])
    return merger(threshold=threshold,trackstats=trackstats)
  else:
    return merger(trackstats=trackstats)


def output_merge_statistics(mergefunc,samplefile=None,locusfile=None):
  '''
  Retrieve merge statistics from the object and write to files by sample/locus

  @param  mergefunc: function to merge multiple genotypes into a consensus genotype.
  @type   mergefunc: callable
  @param samplefile: output file name or file object for merge stats by sample
  @type  samplefile: str or file object
  @param  locusfile: output file name or file object for merge stats by locus
  @type   locusfile: str or file object

  >>> model = model_from_alleles('AB')
  >>> missing = model[None,None]
  >>> A,B='A','B'
  >>> AA,AB,BB = model[A,A],model[A,B],model[B,B]

  >>> import StringIO
  >>> samples  = ['s1','s1','s1','s1','s2','s2','s2','s2','s3','s3','s3','s3','s4','s4','s4','s4']
  >>> loci     = ['l1','l2','l3','l4','l1','l2','l3','l4','l1','l2','l3','l4','l1','l2','l3','l4']
  >>> genosets = [[AA,AA,missing,AA,AB,AA],[AA,AB,AB],[AA,AA,AA],[missing,missing],
  ...            [missing],[AA],[AA,missing,missing,missing],[AB,AB,AA,'BB',AB,missing,missing],
  ...            [AB],[AA],[AA,AA,AA,AA,AB],[AB,AB],
  ...            [AA,AB,AB],[AA,AA],[AA,missing,missing,missing],[AA,AA,AA]]

  >>> merger = GenotypeMerger(Vote(threshold = 0.8))
  >>> list(starmap(merger.merge_geno, izip(samples,loci,repeat(model),genosets))) # doctest: +NORMALIZE_WHITESPACE
  [('A', 'A'), (None, None), ('A', 'A'), (None, None), (None, None), ('A', 'A'),
   ('A', 'A'), (None, None), ('A', 'B'), ('A', 'A'), ('A', 'A'), ('A', 'B'),
   (None, None), ('A', 'A'), ('A', 'A'), ('A', 'A')]

  >>> sampleout = StringIO.StringIO()
  >>> locusout  = StringIO.StringIO()
  >>> output_merge_statistics(merger,sampleout,locusout)
  >>> print sampleout.getvalue() # doctest: +NORMALIZE_WHITESPACE
  SAMPLE_ID   UNAMBIGUOUS_COUNT       CONCORDANT_COUNT        CONSENSUS_COUNT DISCORDANT_COUNT        DISCONCORDANCE_RATE     MISSING_COUNT     CONSENSUS_MISSING_RATE
  s3  2       1       1       0       0.000000        0       0.000000
  s1  0       1       1       1       0.333333        1       0.500000
  s4  1       2       0       1       0.333333        0       0.250000
  s2  2       0       0       1       1.000000        1       0.500000
  >>> print locusout.getvalue() # doctest: +NORMALIZE_WHITESPACE
  LOCUS_ID    UNAMBIGUOUS_COUNT       CONCORDANT_COUNT        CONSENSUS_COUNT DISCORDANT_COUNT        DISCONCORDANCE_RATE     MISSING_COUNT     CONSENSUS_MISSING_RATE
  l3  2       1       1       0       0.000000        0       0.000000
  l4  0       2       0       1       0.333333        1       0.500000
  l1  1       0       1       1       0.500000        1       0.500000
  l2  2       1       0       1       0.500000        0       0.250000
  '''
  if samplefile is not None:
    f = table_writer(samplefile,dialect='tsv')
    f.writerow(SAMPLE_CONCORDANCE_HEADER)
    f.writerows(sorted(build_concordance_output(mergefunc.samplestats),key=itemgetter(5,0)))

  if locusfile is not None:
    fl = table_writer(locusfile,dialect='tsv')
    fl.writerow(LOCUS_CONCORDANCE_HEADER)
    fl.writerows(sorted(build_concordance_output(mergefunc.locusstats),key=itemgetter(5,0)))


def build_concordance_output(stats):
  '''
  Calculate concordance rate and format all the results for output

  @param  stats: a list of merge stats for each sample/locus on the genotypes that was processed
  @type   stats: dict from str -> list
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


def _test():
  import doctest
  return doctest.testmod()


if __name__=='__main__':
  _test()
