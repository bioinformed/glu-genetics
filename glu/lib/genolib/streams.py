# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GLU genotype data management objects'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import optparse

from   types             import NoneType
from   operator          import itemgetter, getitem
from   collections       import defaultdict
from   itertools         import izip,ifilter,imap,chain,groupby,repeat

from   glu.lib.utils     import as_set,tally,izip_exact,gcdisabled,is_str
from   glu.lib.imerge    import imerge
from   glu.lib.xtab      import xtab,rowsby

from   glu.lib.genolib.reprs     import snp
from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.phenos    import Phenome,merge_phenome_list
from   glu.lib.genolib.transform import GenoTransform, prove_unique_transform
from   glu.lib.genolib.merge     import UniqueMerger, VoteMerger
from   glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,Genotype,   \
                                        GenotypeLookupError, GenotypeRepresentationError, \
                                        model_from_alleles,pick,pick_columns,place_list

# Debugging flag
DEBUG=False


class GenotypeStream(object):
  __slots__ = []
  def __init__(self, stream):
    self.__private_stream_do_not_touch = stream

  def __iter__(self):
    '''
    Returns the embedded genotriple stream and marks it as used (ie
    unavailable for further operations) if not already used.  Otherwise,
    raises a RuntimeError exception.

    @return: genotriple stream
    @rtype : sequence of sample, locus, and genotype
    '''
    return iter(self.use_stream())

  def used_stream(self):
    '''
    Returns True if the genotriple stream has been used as an iterable and is
    thus no longer available.  Otherwise, returns False.

    @return: availability of the genotriple stream
    @rtype : bool
    '''
    return self.__private_stream_do_not_touch is None

  def use_stream(self):
    '''
    Returns the embedded genotriple stream and marks it as used (ie
    unavailable for further operations) if not already used.  Otherwise,
    raises a RuntimeError exception.

    @return: genotriple stream
    @rtype :  sequence of sample, locus, and genotype
    '''
    if self.used_stream():
      raise RuntimeError('Genotriple stream already used')

    if not self.materialized:
      self.__private_stream_do_not_touch,genos = None,self.__private_stream_do_not_touch
    else:
      genos = self.__private_stream_do_not_touch

    return genos


class GenotripleStream(GenotypeStream):
  '''
  A stream of genotriples with optional metadata
  '''
  format = 'genotriple'

  def __init__(self, triples, samples=None, loci=None, genome=None, phenome=None,
                              order=None, unique=False, materialized=False):
    '''
    Create a new GenotripleStream object

    Optional metadata on the stream can be supplied, including the samples,
    loci, ordering, and uniqueness of each genotype.  These metadata are
    used to optimize many operations and must be accurate or else incorrect
    results are virtually guaranteed.  When in doubt, do not specify them,
    as each algorithm can compensate, although this may require full
    materialization of the data.

    GenotripleStream objects do support data sources that are materialized
    via the materialized flag.  These behave identically to a non-
    materialized stream except that it is not marked as being used after
    many operations that would normally consume a non-materialized stream.
    Conversely, a materialized GenotripleStream only supports streaming
    operatings, with no additional random-access features of a true
    materialized class.

    @param      triples: sequence of genotriples(str,str,genotype representation)
    @type       triples: sequence
    @param      samples: optional set of samples refered to by the triples
    @type       samples: sequence, set, or None
    @param         loci: optional set of loci refered to by the triples
    @type          loci: sequence, set, or None
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    @param        order: sort order, 'sample' or 'locus', or None, Default is None
    @type         order: str or None

    @type        unique: bool
    @param materialized: flag indicating if genos is a materialized
                         (multiply iterable) data type
    @type  materialized: bool
    '''
    assert genome is not None
    assert isinstance(genome,Genome)
    assert phenome is not None
    assert isinstance(phenome,Phenome)

    if order not in (None,'locus','sample'):
      raise ValueError('invalid GenotripleStream order specified')

    GenotypeStream.__init__(self, triples)

    if not isinstance(samples,(NoneType,set)):
      samples = set(samples)

    if not isinstance(loci,(NoneType,set)):
      loci = set(loci)

    self.samples      = samples
    self.loci         = loci
    self.order        = order
    self.genome       = genome
    self.phenome      = phenome
    self.unique       = bool(unique)
    self.materialized = materialized or isinstance(triples, (list,tuple))

  def _model_pairs(self):
    get_model = self.genome.get_model
    return ( (locus,get_model(locus)) for locus in self.genome.loci or [] )

  model_pairs = property(_model_pairs)

  @staticmethod
  def from_streams(genos, mergefunc=None, order=None):
    '''
    Combine multiple genostreams into one genotriple stream

    @param     genos: genostreams
    @type      genos: list
    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @param     order: sort order, 'sample' or 'locus', or None, Default is None
    @type      order: str or None
    @return         : combined genotriple stream
    @rtype          : sequence of sample, locus, and genotype

    >>> trip1 = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...          ('s2','l2', ('A', 'T')),('s3','l1', ('T', 'T'))]
    >>> trip1 = GenotripleStream.from_tuples(trip1)
    >>> trip2 = [('s2','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> trip2 = GenotripleStream.from_tuples(trip2)
    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)

    >>> streams = [trip1,genos,trip2]
    >>> combined = GenotripleStream.from_streams(streams, mergefunc=VoteMerger())
    >>> for row in combined:
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', (None, None))
    ('s2', 'l2', (None, None))
    ('s3', 'l1', ('T', 'T'))
    ('s3', 'l2', (None, None))
    '''
    if not genos:
      raise ValueError('empty stream list')

    # FIXME: Consider merging before converting to triples if inputs are all
    #        sdat/ldat
    triples = [ g.as_genotriples() for g in genos ]

    if len(triples) == 1:
      triples = triples[0]

    else:
      if mergefunc is not None and order is None:
        order = 'sample'

      # FIXME: Incore sort is deadly
      if order is not None:
        triples = [ t.sorted(order) for t in triples ]

      if order is None and mergefunc is None:
        triples = combine_unsorted_genotriple_list(triples)

      else:
        triples = combine_sorted_genotriple_list(triples)

    if mergefunc is not None:
      triples = triples.merged(mergefunc)

    return triples

  @staticmethod
  def from_tuples(triples, samples=None, loci=None, order=None, unique=False,
                           genome=None, phenome=None):
    '''
    Alternate constructor that builds a new GenotripleStream object from a
    sequence of triples with genotypes in tuple

    @param      triples: sequence of genotriples(str,str,genotype representation)
    @type       triples: sequence
    @param      samples: optional set of samples refered to by the triples. Default is None
    @type       samples: sequence, set, or None
    @param         loci: optional set of loci refered to by the triples. Default is None
    @type          loci: sequence, set, or None
    @param        order: sort order, 'sample' or 'locus', or None, Default is None
    @type         order: str or None
    @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
    @type        unique: bool
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples)
    >>> for row in triples:
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', ('G', 'T'))
    ('s2', 'l2', ('T', 'T'))
    ('s3', 'l1', ('G', 'G'))
    ('s3', 'l2', ('A', 'A'))
    '''
    if genome is None:
      genome = Genome()

    if phenome is None:
      phenome = Phenome()

    triples = encode_genotriples_from_tuples(triples, genome)
    return GenotripleStream(triples, samples=samples, loci=loci, order=order, genome=genome, phenome=phenome,
                                     unique=unique)

  @staticmethod
  def from_strings(triples, genorepr, samples=None, loci=None, order=None, unique=False,
                            genome=None, phenome=None):
    '''
    Alternate constructor that builds a new GenotripleStream object from a
    sequence of triples with genotypes in a string format

    @param      triples: sequence of genotriples(str,str,genotype representation)
    @type       triples: sequence
    @param     genorepr: internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    @param      samples: optional set of samples refered to by the triples. Default is None
    @type       samples: sequence, set, or None
    @param         loci: optional set of loci refered to by the triples. Default is None
    @type          loci: sequence, set, or None
    @param        order: sort order, 'sample' or 'locus', or None, Default is None
    @type         order: str or None
    @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
    @type        unique: bool
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance

    >>> triples = [('l1','s1','AA'),('l1','s1','  '),('l1','s2','AB'),('l2','s1','AA'),
    ...            ('l2','s1','AA'),('l3','s1','BB'),('l3','s1','BB'),('l3','s1','AB')]
    >>> triples = GenotripleStream.from_strings(triples,snp)
    >>> for row in triples:
    ...   print row
    ('l1', 's1', ('A', 'A'))
    ('l1', 's1', (None, None))
    ('l1', 's2', ('A', 'B'))
    ('l2', 's1', ('A', 'A'))
    ('l2', 's1', ('A', 'A'))
    ('l3', 's1', ('B', 'B'))
    ('l3', 's1', ('B', 'B'))
    ('l3', 's1', ('A', 'B'))
    '''
    if genome is None:
      genome = Genome()

    if phenome is None:
      phenome = Phenome()

    triples = encode_genotriples_from_strings(triples, genorepr, genome)
    return GenotripleStream(triples, samples=samples, loci=loci, order=order, genome=genome, phenome=phenome,
                                     unique=unique)

  def to_tuples(self):
    '''
    Iterate over genotriples with genotypes coded as tuples

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples).to_tuples()
    >>> for row in triples:
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', ('G', 'T'))
    ('s2', 'l2', ('T', 'T'))
    ('s3', 'l1', ('G', 'G'))
    ('s3', 'l2', ('A', 'A'))
    '''
    for sample,locus,geno in self:
      yield sample,locus,geno.alleles()

  def to_strings(self,genorepr):
    '''
    Iterate over genotriples with genotypes coded as strings

    @param genorepr: internal representation of genotypes
    @type  genorepr: UnphasedMarkerRepresentation or similar object

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples).to_strings(snp)
    >>> for row in triples:
    ...   print row
    ('s1', 'l1', 'GG')
    ('s1', 'l2', 'AA')
    ('s2', 'l1', 'GT')
    ('s2', 'l2', 'TT')
    ('s3', 'l1', 'GG')
    ('s3', 'l2', 'AA')
    '''
    repr = genorepr.to_string
    for sample,locus,geno in self:
      yield sample,locus,repr(geno)

  def clone(self, triples, **kwargs):
    '''
    Alternative constructor that builds a new GenotripleStream object
    with attributes based on self, but updated with the specified keyword
    arguments.

    @param      triples: sequence of genotriples(str,str,genotype representation)
    @type       triples: sequence
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param        order: sort order, either 'samples', 'locus'
    @type         order: str
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param materialized: flag indicating if this stream is materialized and
                         allows iteration multiple times
    @type  materialized: bool

    For example:

    stream.clone(triples, materialized=False)

    is equivalent to

    GenotripleStream(triples, samples=stream.samples, loci=stream.loci,
                              order=stream.order, genome=stream.genome,
                              phenome=stream.phenome, unique=stream.unique,
                              materialized=False)
    '''
    kwargs.setdefault('samples',      self.samples)
    kwargs.setdefault('loci',         self.loci)
    kwargs.setdefault('order',        self.order)
    kwargs.setdefault('genome',       self.genome)
    kwargs.setdefault('phenome',      self.phenome)
    kwargs.setdefault('unique',       self.unique)
    kwargs.setdefault('materialized', self.materialized)

    return GenotripleStream(triples, **kwargs)

  def materialize(self):
    '''
    Returns a materialized genotriple stream.

    GenotripleStream objects do support data sources that are materialized
    via the materialized flag.  These behave identically to a non-
    materialized stream except that it is not marked as being used after
    many operations that would normally consume a non-materialized stream.
    Conversely, a materialized GenomatrixStream only supports streaming
    operatings, with no additional random-access features of a true
    materialized class.
    '''
    if self.materialized:
      return self

    genos   = list(self.use_stream())

    # FIXME: These extraction operators should work in parallel to avoid thrashing
    samples = set(imap(itemgetter(0),genos))
    loci    = set(imap(itemgetter(1),genos))

    return self.clone(genos, loci=loci, samples=samples, materialized=True)

  def transformed(self, transform=None, mergefunc=None, **kwargs):
    '''
    Apply filtering and renaming transformations to a genotriple stream.
    Transformations may be specified by key-word arguments or a
    transformation object, though not both.  Inclusions and exclusions are
    always performed before renaming operations.  If a merge function is
    specified when performing renaming transformations, the resulting
    genotriple stream will be guaranteed to contain unique rows and columns
    (see merged function).

    @param       transform: transformation object (optional)
    @type        transform: GenoTransform object
    @param       mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type        mergefunc: callable
    @param include_samples: filter samples such that they must appear in the set (optional)
    @type  include_samples: set
    @param exclude_samples: filter samples such that they must not appear in the set (optional)
    @type  exclude_samples: set
    @param    include_loci: filter loci such that they must appear in the set (optional)
    @type     include_loci: set
    @param    exclude_loci: filter loci such that they must not appear in the set (optional)
    @type     exclude_loci: set
    @param  rename_samples: rename any samples that appear in the supplied dictionary to the
                            associated value (optional)
    @type   rename_samples: dict from str -> str
    @param     rename_loci: rename any loci that appear in the supplied dictionary to the
                            associated value (optional)
    @type      rename_loci: dict from str -> str
    @param   order_samples: reorder samples such based on the order of the supplied list (optional)
    @type    order_samples: list
    @param      order_loci: reorder loci such based on the order of the supplied list (optional)
    @type       order_loci: list
    @param  filter_missing: filter missing genotypes from the stream
    @type   filter_missing: bool
    @return               : transformed genotriple stream
    @rtype                : GenotripleStream

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples)
    >>> for row in triples.transformed(include_loci=['l1'],exclude_samples=['s3']):
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s2', 'l1', ('G', 'T'))
    '''
    if transform and kwargs:
      raise ValueError('Ambiguous transformation specification')

    if transform:
      if isinstance(transform, optparse.OptionContainer):
        transform = GenoTransform.from_options(options)
      elif isinstance(transform, dict):
        transform = GenoTransform.from_kwargs(transform)
      elif not isinstance(transform, GenoTransform):
        raise ValueError('Invalid genotype transformation specification')
    elif kwargs:
      transform = GenoTransform.from_kwargs(**kwargs)

    if transform is None:
      return self

    triples = self

    # Filter missing genotyes
    if transform.filter_missing_genotypes:
      triples = filter_genotriples_missing(triples)

    # Sample and locus includes
    if transform.samples.include is not None or transform.loci.include is not None:
      triples = filter_genotriples(triples,transform.samples.include,transform.loci.include)

    # Sample and locus excludes
    if transform.samples.exclude is not None or transform.loci.exclude is not None:
      triples = filter_genotriples(triples,transform.samples.exclude,transform.loci.exclude,exclude=True)

    # Determine if resulting triples will be unique (before renaming samples and loci)
    triples.unique = prove_unique_transform(transform=transform,loci=triples.loci,samples=triples.samples,unique=self.unique)

    # Sample and locus renaming
    if transform.samples.rename is not None or transform.loci.rename is not None:
      triples = rename_genotriples(triples,transform.samples.rename,transform.loci.rename)

    # FIXME: recoding and renaming can be combined
    if transform.rename_alleles is not None:
      triples = rename_genotriples_alleles(triples, transform.rename_alleles)
    if transform.recode_models is not None:
      triples = recode_genotriples(triples, transform.recode_models)

    # FIXME: I'm not sure what I meant here, but it doesn't look quite right
    if transform.samples.order is not None and transform.loci.order is not None:
      if transform.samples.order is not None:
        order = 'samples'
      else:
        order = 'loci'
      triples = triples.sorted(order,sampleorder=transform.samples.order,locusorder=transform.loci.order)

    # If the results are not unique and a merge function is specified, perform a merge
    if mergefunc is not None:
      triples = triples.merged(mergefunc)

    return triples

  def sorted(self, order='sample', locusorder=None, sampleorder=None):
    '''
    Returns a new GenotripleStream ordered by sample or locus, where samples
    and loci are placed in the specified order or sorted lexicographically.

    @param       order: sort order, either 'samples', 'locus'.  Default is 'sample'
    @type        order: str
    @param  locusorder: ordered list of loci (optional)
    @type   locusorder: sequence
    @param sampleorder: ordered list of samples (optional)
    @type  sampleorder: sequence
    @return           : sorted genotriple stream
    @rtype            : GenotripleStream

    >>> triples = [('s3','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s3','l2', ('T', 'T')),
    ...            ('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples)
    >>> for row in triples.sorted():
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', ('G', 'T'))
    ('s2', 'l2', ('A', 'A'))
    ('s3', 'l1', ('G', 'G'))
    ('s3', 'l2', ('T', 'T'))
    '''
    if self.order is not None and self.order == order:
      return self

    return sort_genotriples(self,order=order,locusorder=locusorder,sampleorder=sampleorder)

  def merged(self, mergefunc=None, order='sample'):
    '''
    Returns a new sorted genotriple stream with all genotypes for the same
    sample and locus merged into a consensus using the supplied merge
    function.  The optional order argument is used to determine a preferred
    sorting order for the triples, if the stream is not already sorted.
    Thus, one must not expect that the resulting stream be in the suggested
    order.

    @param  mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type   mergefunc: callable
    @param      order: sort order, either 'samples', 'locus'.  Default is 'sample'
    @type       order: str
    @return          : sorted and merged genotriple stream
    @rtype           : GenotripleStream

    >>> triples = [('s1','l1', ('T', 'T')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('A', 'A')),
    ...            ('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples)
    >>> for row in triples.merged(VoteMerger(),order='locus'):
    ...   print row
    ('s1', 'l1', (None, None))
    ('s2', 'l1', ('G', 'T'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l2', ('A', 'A'))
    '''
    if self.unique:
      return self

    genos = self

    if genos.order not in ('sample','locus'):
      genos = genos.sorted(order)

    return merge_sorted_genotriples(genos,mergefunc)

  def as_genotriples(self):
    '''
    Return the current genotriple stream as a genotriple stream.  i.e., returns self

    @return: current genotriple stream
    @rtype : GenotripleStream
    '''
    return self

  def as_ldat(self, mergefunc=None):
    '''
    Return the current genotriple data as a GenomatrixStream by locus.
    Ordered triple streams can be transformed without full materialization,
    though unordered streams do require materialization.  A merge function is
    required for non-unique streams.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> triples = [('s1','l1',('G', 'G')),('s1','l2',('A', 'A')),
    ...            ('s2','l1',('G', 'T')),('s2','l2',('T', 'T')),
    ...            ('s3','l1',('G', 'G')),('s3','l2',('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples).materialize()

    >>> merge = VoteMerger()
    >>> ldat = triples.as_ldat(merge)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = triples.clone(triples,loci=['l1','l2'],samples=None).as_ldat(merge)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = triples.clone(triples,samples=['s1','s2','s3'],loci=None).as_ldat(merge)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
    '''
    if mergefunc is None:
      mergefunc = UniqueMerger()

    return build_genomatrixstream_from_genotriples(self, 'ldat', mergefunc=mergefunc)

  def as_sdat(self, mergefunc=None):
    '''
    Return the current genotriple data as a GenomatrixStream by sample.
    Ordered triple streams can be transformed without full materialization,
    though unordered streams do require materialization.  A merge function is
    required for non-unique streams.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples).materialize()

    >>> merge = VoteMerger()
    >>> sdat = triples.as_sdat(merge)
    >>> sdat.loci
    ('l1', 'l2')
    >>> sdat.samples
    ('s1', 's2', 's3')
    >>> for row in sdat:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])

    >>> sdat = triples.clone(triples,loci=['l1','l2'],samples=None).as_sdat(merge)
    >>> sdat.loci
    ('l1', 'l2')
    >>> sdat.samples
    ('s1', 's2', 's3')
    >>> for row in sdat:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])

    >>> sdat = triples.clone(triples,samples=['s1','s2','s3'],loci=None).as_sdat(merge)
    >>> sdat.loci
    ('l1', 'l2')
    >>> sdat.samples
    ('s1', 's2', 's3')
    >>> for row in sdat:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])
    '''
    if mergefunc is None:
      mergefunc = UniqueMerger()

    return build_genomatrixstream_from_genotriples(self, 'sdat', mergefunc=mergefunc)


class GenomatrixStream(GenotypeStream):
  '''
  A stream of genomatrix by sample or locus with optional metadata
  '''
  def __init__(self, genos, format, samples=None, loci=None, models=None,
                     genome=None, phenome=None, unique=True, materialized=False,
                     packed=False):
    '''
    Create a new GenomatrixStream object

    Metadata on the stream can be supplied, including the samples, loci,
    ordering, and uniqueness of each genotype. Samples and loci are a list
    of identifiers that correspond to the precise order of the rows or
    columns.  ldat formated matrix streams must have samples specified, sdat
    formated matrix streams must have loci specified.  Otherwise, the
    optional metadata are used to optimize many operations and must be
    accurate or else incorrect results are virtually guaranteed.  When in
    doubt, do not specify them, as each algorithm can compensate, although
    this may require full materialization of the data.

    GenomatrixStream objects support genotypes that are materialized, like
    lists of rows, via the materialized flag.  These behave identically to a
    non- materialized stream except that it is not marked as being used
    after many operations that would normally consume a non-materialized
    stream.  Conversely, a materialized GenomatrixStream only supports
    streaming operatings, with no additional random-access features of a
    true materialized class.

    @param        genos: genomatrix stream
    @type         genos: sequence
    @param       format: format of input genomatrix, either 'ldat' or 'sdat'
    @type        format: str
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param       models: new internal representation of genotypes
    @type        models: UnphasedMarkerRepresentation or similar object
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param materialized: flag indicating if this stream is materialized and
                         allows iteration multiple times
    @type  materialized: bool
    @param       packed: flag indicating if genotypes are packed into a
                         compressed array format
    @type        packed: bool
    '''
    assert models is not None
    assert format=='ldat' or len(models) == len(loci)

    assert genome is not None
    assert isinstance(genome,Genome)
    assert phenome is not None
    assert isinstance(phenome,Phenome)

    if phenome is None:
      phenome = Phenome()

    if format not in ('sdat','ldat'):
      raise ValueError("Invalid genomatrix format '%s'.  Must be either sdat or ldat" % format)

    GenotypeStream.__init__(self, genos)

    if loci and not isinstance(loci,tuple):
      loci = tuple(loci)

    if samples and not isinstance(samples,tuple):
      samples = tuple(samples)

    if format=='ldat':
      if samples is None:
        raise ValueError('GenomatrixStream in ldat format requires sample metadata')
    elif format=='sdat':
      if loci is None:
        raise ValueError('GenomatrixStream in sdat format requires locus metadata')

    materialized = materialized or isinstance(genos, (list,tuple))

    # If materialized, we can supplement our metadata and perform additional
    # sanity checking
    if materialized:
      rowlabels = tuple(imap(itemgetter(0),genos))

      if format == 'sdat':
        if samples is not None and samples!=rowlabels:
          raise ValueError('GenomatrixStream row labels do not match samples')
        samples = rowlabels

      elif format == 'ldat':
        if loci is not None and loci!=rowlabels:
          raise ValueError('GenomatrixStream row labels do not match loci')
        loci = rowlabels

    self.format       = format
    self.samples      = samples
    self.loci         = loci
    self.models       = models
    self.genome       = genome
    self.phenome      = phenome
    self.unique       = unique
    self.materialized = materialized
    self.packed       = packed

  if DEBUG:
    def __iter__(self):
      '''
      Returns the embedded genotriple stream and marks it as used (ie
      unavailable for further operations) if not already used.  Otherwise,
      raises a RuntimeError exception.

      @return: genotriple stream
      @rtype :  sequence of sample, locus, and genotype
      '''
      if self.format=='ldat':
        def _check(rows):
          for (locus,row),model in izip_exact(rows,self.models):
            assert not self.packed or isinstance(row,GenotypeArray)
            assert not self.packed or model is row.descriptor.models[0]
            assert all(g.model is model for g in row)
            assert model is self.genome.loci[locus].model
            assert all(isinstance(g,Genotype) for g in row)
            yield locus,row
        return _check(self.use_stream())
      else:
        assert all(self.genome.loci[locus].model in (model,None) for locus,model in izip_exact(self.loci,self.models))
        def _check(rows):
          for sample,row in rows:
            assert all(isinstance(g,Genotype) for g in row)
            assert len(self.models) == len(row)
            assert all(g.model is model for model,g in izip_exact(self.models,row))
            yield sample,row
        return _check(self.use_stream())

      return iter(self.use_stream())

  @staticmethod
  def from_streams(genos, format, mergefunc=None):
    '''
    Combine multiple genostreams into one genomatrix stream

    @param     genos: genostreams
    @type      genos: list
    @param    format: format of input genomatrix, either 'ldat' or 'sdat'
    @type     format: str
    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @return         : combined genotriple stream
    @rtype          : sequence of sample, locus, and genotype

    >>> trip1 = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...          ('s2','l2', ('T', 'T')),('s3','l1', ('T', 'T'))]
    >>> trip1 = GenotripleStream.from_tuples(trip1).materialize()
    >>> trip2 = [('s2','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> trip2 = GenotripleStream.from_tuples(trip2).materialize()

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples).materialize()

    >>> streams = [trip1,genos,trip2]
    >>> merger=VoteMerger()
    >>> combined = GenomatrixStream.from_streams(streams, 'ldat', mergefunc=merger)
    >>> combined.columns
    ('s1', 's2', 's3')
    >>> for row in combined:
    ...   print row
    ('l1', [('G', 'G'), (None, None), ('T', 'T')])
    ('l2', [('A', 'A'), ('T', 'T'), (None, None)])

    >>> streams = [genos,genos]
    >>> combined = GenomatrixStream.from_streams(streams, 'sdat', mergefunc=merger)
    >>> combined.columns
    ('l1', 'l2')
    >>> for row in combined:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('T', 'T'), ('A', 'T')])
    '''
    if format not in ('sdat','ldat'):
      raise ValueError("Invalid genomatrix format '%s'.  Must be either sdat or ldat" % format)

    formats  = set(g.format for g in genos)
    informat = list(formats)[0] if len(formats)==1 else None
    headers  = [ g.columns for g in genos ] if informat in ('ldat','sdat') else None

    # Single input is trivial -- just merge
    if len(genos) == 1:
      genos = genos[0].transformed(mergefunc=mergefunc)

    # Inputs are all matricies, without knowledge of uniqueness of rows.
    elif formats <= set(['ldat','sdat']):
      genos = [ (g.as_ldat() if format=='ldat' else g.as_sdat()) for g in genos ]
      genos = merge_genomatrixstream_list(genos, mergefunc)

    # Very general strategy of converting all inputs to genotriples, sorting
    # by the appopriate clustering, and transforming into ldat or sdat.
    # Unfortunately, very slow.
    else:
      order = 'locus' if format=='ldat' else 'sample'
      genos = GenotripleStream.from_streams(genos, order=order)

    if format == 'ldat':
      genos = genos.as_ldat(mergefunc)
    else:
      genos = genos.as_sdat(mergefunc)

    return genos

  def _get_rows(self):
    return self.loci if self.format == 'ldat' else self.samples

  def _set_rows(self, rows):
    if self.format == 'ldat':
      self.loci    = rows
    else:
      self.samples = rows

  def _get_columns(self):
    return self.loci if self.format == 'sdat' else self.samples

  def _set_columns(self, columns):
    if self.format == 'sdat':
      self.loci    = columns
    else:
      self.samples = columns

  @staticmethod
  def from_tuples(genos, format, samples=None, loci=None, unique=True,
                         genome=None, phenome=None):
    '''
    Alternate constructor that builds a new GenomatrixStream object from a
    genotype matrix stream with genotypes in a string format.

    @param        genos: genomatrix stream
    @type         genos: sequence
    @param       format: format of input genomatrix, either 'ldat' or 'sdat'
    @type        format: str
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
    >>> for row in genos:
    ...  print row
    ('l1', [('G', 'G'), ('G', 'T'), ('T', 'T')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
    '''
    if format=='ldat':
      columns = samples
    elif format=='sdat':
      columns = loci
    else:
      raise ValueError('Invalid genotype matrix format')

    if phenome is None:
      phenome = Phenome()

    columns,models,genome,genos = encode_genomatrixstream_from_tuples(columns,genos,format,
                                   genome=genome,unique=unique)
    return GenomatrixStream(genos, format, samples=samples, loci=loci, models=models,
                                   genome=genome, phenome=phenome, unique=unique,
                                   packed=True, materialized=False)

  @staticmethod
  def from_strings(genos, format, genorepr, samples=None, loci=None, unique=True,
               packed=False, genome=None, phenome=None):
    '''
    Alternate constructor that builds a new GenomatrixStream object from a
    genotype matrix stream with genotypes in a string format.

    @param        genos: genomatrix stream
    @type         genos: sequence
    @param       format: format of input genomatrix, either 'ldat' or 'sdat'
    @type        format: str
    @param     genorepr: internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param       packed: flag indicating if genotypes are packed into a
                         compressed array format
    @type        packed: bool
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    '''
    if format=='ldat':
      columns = samples
    elif format=='sdat':
      columns = loci
    else:
      raise ValueError('Invalid genotype matrix format')

    if phenome is None:
      phenome = Phenome()

    columns,models,genome,genos = encode_genomatrixstream_from_strings(columns,genos,format,genorepr,
                                   genome=genome,unique=unique)
    return GenomatrixStream(genos, format, samples=samples, loci=loci, models=models,
                                   genome=genome, phenome=phenome, unique=unique,
                                   packed=True, materialized=False)

  def to_tuples(self):
    '''
    Iterate over genotypes coded as tuples

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples).to_tuples()
    >>> for row in genos:
    ...  print row
    ('l1', [('G', 'G'), ('G', 'T'), ('T', 'T')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
    '''
    for label,row in self:
      yield label,[ g.alleles() for g in row ]

  def to_strings(self,genorepr):
    '''
    Iterate over genotypes coded as strings

    @param     genorepr: internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples).to_strings(snp)
    >>> for row in genos:
    ...  print row
    ('l1', ['GG', 'GT', 'TT'])
    ('l2', ['AA', 'TT', 'AT'])
    '''
    repr = genorepr.to_strings
    for label,row in self:
      yield label,repr(row)

  def clone(self, genos, **kwargs):
    '''
    Alternative constructor that builds a new GenomatrixStream object
    with a new data stream and attributes based on self, but updated with
    specified keyword arguments.

    @param        genos: genomatrix stream
    @type         genos: sequence
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param       models: new internal representation of genotypes
    @type        models: UnphasedMarkerRepresentation or similar object
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param       genome: genome descriptor
    @type        genome: Genome instance
    @param      phenome: phenome descriptor
    @type       phenome: Phenome instance
    @param materialized: flag indicating if this stream is materialized and
                         allows iteration multiple times
    @type  materialized: bool
    @param       packed: flag indicating if genotypes are packed into a
                         compressed array format
    @type        packed: bool

    For example:

      stream.clone(genos, materialized=False)

    is equivalent to

      GenomatrixStream(genos, format=stream.format, samples=stream.samples,
                              loci=stream.loci, model=stream.models,
                              genome=stream.genome, phenome=stream.phenome,
                              unique=stream.unique, materialized=stream.materialized,
                              packed=stream.packed)
    '''
    kwargs.setdefault('format',       self.format)
    kwargs.setdefault('samples',      self.samples)
    kwargs.setdefault('loci',         self.loci)
    kwargs.setdefault('models',       self.models)
    kwargs.setdefault('genome',       self.genome)
    kwargs.setdefault('phenome',      self.phenome)
    kwargs.setdefault('unique',       self.unique)
    kwargs.setdefault('materialized', self.materialized)
    kwargs.setdefault('packed',       self.packed)

    return GenomatrixStream(genos, **kwargs)

  rows    = property(_get_rows,   _set_rows)
  columns = property(_get_columns,_set_columns)

  def _model_pairs(self):
    return izip(self.loci,self.models)

  model_pairs = property(_model_pairs)

  def __len__(self):
    '''
    Return the number of rows of a GenotypematrixStream, if known.
    Otherwise a TypeError is raised.

    @return: number of rows
    @rtype : int

    >>> samples = ('s1','s2','s3')
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
    >>> len(genos)
    Traceback (most recent call last):
       ...
    TypeError: Unknown GenomatrixStream length

    >>> genos = genos.materialize()
    >>> len(genos)
    2
    >>> len(genos.transposed())
    3
    '''
    if self.materialized:
      # Does not invalidate stream since we are materialized
      assert len(self.rows) == len(self.use_stream())
      return len(self.use_stream())
    elif self.rows is not None:
      return len(self.rows)

    raise TypeError('Unknown GenomatrixStream length')

  def __getitem__(self, index):
    '''
    Return one or more rows of a materialized GenotypematrixStream.  If
    index is an integer, the corresponding row label and genotype data are
    returned using Python index semantics.  Similarly if a slice is
    specified, a sequence of row labels and genotype data are returned.

    @param index: requested index or slice
    @type  index: integer or slice object
    @return     : integer index: a tuple of row label and genotype array
                  slice: a list of row label and genotype array tuples

    >>> samples = ('s1','s2','s3')
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)

    >>> genos[1]
    Traceback (most recent call last):
       ...
    IndexError: Random access required a materialized GenomatrixStream

    >>> genos = genos.materialize()
    >>> genos[1]
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
    >>> genos[-1]
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
    >>> genos[:1]
    [('l1', [('G', 'G'), ('G', 'T'), ('T', 'T')])]
    >>> genos[1:]
    [('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])]
    '''
    if not self.materialized:
      raise IndexError('Random access required a materialized GenomatrixStream')

    # Does not invalidate stream since we are materialized
    return self.use_stream()[index]

  def materialize(self):
    '''
    Returns a materialized genomatrix stream.

    GenomatrixStream objects do support data sources that are materialized
    via the materialized flag.  These behave identically to a non-
    materialized stream except that it is not marked as being used after
    many operations that would normally consume a non-materialized stream.
    Conversely, a materialized GenomatrixStream only supports streaming
    operatings, with no additional random-access features of a true
    materialized class.
    '''
    if self.materialized:
      return self

    genos = self.transformed(repack=True)
    return genos.clone(list(genos.use_stream()), materialized=True)

  def transformed(self, transform=None, mergefunc=None, **kwargs):
    '''
    Apply filtering and renaming transformations to a genomatrix stream.
    Transformations may be specified by key-word arguments or a
    transformation object, though not both.  Inclusions and exclusions are
    always performed before renaming operations.  If a mergefunction is
    specified when performing renaming transformations, the resulting
    genomatrix stream will be guaranteed to contain unique rows and columns
    (see merged function).

    @param           genos: genomatrix stream
    @type            genos: sequence
    @param       transform: transformation object (optional)
    @type        transform: GenoTransform object
    @param       mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type        mergefunc: callable
    @param include_samples: filter samples such that they must appear in the set (optional)
    @type  include_samples: set
    @param exclude_samples: filter samples such that they must not appear in the set (optional)
    @type  exclude_samples: set
    @param    include_loci: filter loci such that they must appear in the set (optional)
    @type     include_loci: set
    @param    exclude_loci: filter loci such that they must not appear in the set (optional)
    @type     exclude_loci: set
    @param  rename_samples: rename any samples that appear in the supplied dictionary to the
                            associated value (optional)
    @type   rename_samples: dict from str -> str
    @param     rename_loci: rename any loci that appear in the supplied dictionary to the
                            associated value (optional)
    @type      rename_loci: dict from str -> str
    @param   order_samples: reorder samples such based on the order of the supplied list (optional)
    @type    order_samples: list
    @param      order_loci: reorder loci such based on the order of the supplied list (optional)
    @type       order_loci: list
    @param  filter_missing: filter missing genotypes from the stream
    @type   filter_missing: bool
    @return               : possibly materialized genotype matrix
    @rtype                : generator or list

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ ('G', 'G'),   ('G', 'T'),   ('T', 'T') ]),
    ...         ('l2', [ ('A', 'A'),   ('T', 'T'),   ('A', 'T') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
    >>> genos = genos.transformed(include_loci=['l1'],exclude_samples=['s3'])
    >>> genos.samples
    ('s1', 's2')
    >>> for row in genos:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T')])
    '''
    if transform and kwargs:
      raise ValueError('Ambiguous transformation specification')

    if transform:
      if isinstance(transform, optparse.OptionContainer):
        transform = GenoTransform.from_options(options)
      elif isinstance(transform, dict):
        transform = GenoTransform.from_kwargs(transform)
      elif not isinstance(transform, GenoTransform):
        raise ValueError('Invalid genotype transformation specification')
    elif kwargs:
      transform = GenoTransform.from_kwargs(**kwargs)

    if transform is None:
      return self

    if self.format == 'ldat':
      rowtransform = transform.loci
      coltransform = transform.samples
    else:
      rowtransform = transform.samples
      coltransform = transform.loci

    # Apply 6 stage transformation process:
    #
    # 1) Row include and exclude
    #    o update row metadata, since it would be lost otherwise
    # 2) Column include and exclude
    #    o column metadata is readily avaialble via the stream head, so no
    #      need to maintain it
    # 3) Update unique-status of the resulting matrix
    # 4) Rename rows and columns
    #    o again, update row medatadata
    # 5) Adjust representation
    # 6) Order and repack

    genos = self

    # Apply row includes and excludes
    if rowtransform.exclude:
      genos = filter_genomatrixstream_by_row(genos,rowtransform.exclude,exclude=True)
    if rowtransform.include is not None:
      genos = filter_genomatrixstream_by_row(genos,rowtransform.include)

    # Apply column includes and excludes
    if coltransform.exclude:
      genos = filter_genomatrixstream_by_column(genos,coltransform.exclude,exclude=True)
    if coltransform.include is not None:
      genos = filter_genomatrixstream_by_column(genos,coltransform.include)

    genos.unique = prove_unique_transform(transform=transform,loci=genos.loci,samples=genos.samples,unique=genos.unique)

    # Apply renamings
    if rowtransform.rename:
      genos = rename_genomatrixstream_row(genos,rowtransform.rename)
    if coltransform.rename:
      genos = rename_genomatrixstream_column(genos,coltransform.rename)

    # Filter rows and columns with all missing data
    if transform.filter_missing_genotypes:
      genos = filter_genomatrixstream_missing(genos)

    if transform.rename_alleles:
      genos = rename_genomatrixstream_alleles(genos,transform.rename_alleles)

    if transform.recode_models is not None:
      genos = recode_genomatrixstream(genos, transform.recode_models)

    # Ordering by [] and None are distinct cases: the first will order in lexicographical order
    # as a side-effect, while the latter should not alter order.
    if transform.samples.order is not None or transform.loci.order is not None:
      genos = genos.sorted(locusorder=transform.loci.order, sampleorder=transform.samples.order)

    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if transform.repack and not genos.packed:
      genos = pack_genomatrixstream(genos)

    return genos

  def sorted(self, locusorder=None, sampleorder=None):
    '''
    Returns a new GenomatrixStream ordered samples and/or loci.  Any sample
    or locus not appearing in the ordering list will be returned at the end,
    preserving the input ordering.

    @param  locusorder: ordered list of loci (optional)
    @type   locusorder: sequence
    @param sampleorder: ordered list of samples (optional)
    @type  sampleorder: sequence
    @return           : sorted genomatrix stream
    @rtype            : GenomatrixStream

    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ ('G', 'G'),   ('A', 'A') ]),
    ...         ('s2', [ ('G', 'T'),   ('T', 'T') ]),
    ...         ('s3', [ ('G', 'G'),   ('A', 'A') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)

    >>> genos = genos.sorted(locusorder=['l2'],sampleorder=['s2','s1','s3']).materialize()
    >>> genos.loci
    ('l2', 'l1')
    >>> for row in genos:
    ...   print row
    ('s2', [('T', 'T'), ('G', 'T')])
    ('s1', [('A', 'A'), ('G', 'G')])
    ('s3', [('A', 'A'), ('G', 'G')])

    >>> genos = genos.transposed()
    >>> genos = genos.sorted(locusorder=['l2'],sampleorder=['s2','s1','s3'])
    >>> genos.samples
    ('s2', 's1', 's3')
    >>> for row in genos:
    ...   print row
    ('l2', [('T', 'T'), ('A', 'A'), ('A', 'A')])
    ('l1', [('G', 'T'), ('G', 'G'), ('G', 'G')])
    '''
    if self.format == 'sdat':
      columnorder = locusorder
      roworder    = sampleorder
    else:
      columnorder = sampleorder
      roworder    = locusorder

    genos = self
    if roworder is not None:
      genos = reorder_genomatrixstream_rows(genos,roworder)

    if columnorder is not None:
      genos = reorder_genomatrixstream_columns(genos,columnorder)

    return genos

  def merged(self, mergefunc):
    '''
    Merge genotypes for rows and columns with the same labels

    @param       mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type        mergefunc: callable
    '''
    return merge_genomatrixstream(self, mergefunc)

  def unique_checked(self):
    '''
    Return a genomatrix stream that is verified to have unique row and
    column labels

    Non-unique columns:

    >>> genos = GenomatrixStream([],'sdat',loci=['L1','L2','L3','L1'],models=[snp]*4,
    ...                                    genome=Genome(),phenome=Phenome())
    >>> genos.unique_checked()
    Traceback (most recent call last):
         ...
    NonUniqueError: Non-unique loci: L1:2

    Non-unique rows:

    >>> loci=('L1','L2')
    >>> rows=[('R1',['AA','AC']),
    ...       ('R1',['AA','AC'])]
    >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
    >>> genos.unique_checked().materialize()
    Traceback (most recent call last):
         ...
    NonUniqueError: Non-unique row name: R1

    Known unique rows and columns:

    >>> loci=('L1','L2')
    >>> samples=('R1', 'R2')
    >>> rows=[('R1',['AA','AC']),
    ...       ('R2',['AA','AC'])]
    >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci,samples=samples)
    >>> ugenos = genos.unique_checked()
    >>> genos is ugenos
    True

    Known columns, unknown but unique rows:

    >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
    >>> genos = genos.unique_checked()
    >>> for sample,row in genos:
    ...   print sample,row
    R1 [('A', 'A'), ('A', 'C')]
    R2 [('A', 'A'), ('A', 'C')]
    '''
    return unique_check_genomatrixstream(self)

  def transposed(self):
    '''
    Return the transpose of this genomatrix stream; the same genotypes but
    with the rows and columns swapped.  This is also equivalent to toggling
    between ldat and sdat formats.  Take care using this method as it
    requires full materialization of the data.

    @return: transposed genomatrix stream
    @rtype : genomatrix stream

    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ ('G', 'G'),   ('A', 'A') ]),
    ...         ('s2', [ ('G', 'T'),   ('T', 'T') ]),
    ...         ('s3', [ ('G', 'G'),   ('A', 'A') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
    >>> genos.format
    'sdat'
    >>> genos.columns
    ('l1', 'l2')
    >>> genos.rows
    >>> new_genos = genos.transposed()
    >>> new_genos.format
    'ldat'
    >>> new_genos.rows
    ('l1', 'l2')
    >>> new_genos.columns
    ('s1', 's2', 's3')
    >>> for row in new_genos:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
    '''
    genos = self.materialize()
    rows,tgenos = transpose_generator(genos.columns, genos.use_stream())

    if genos.format == 'ldat':
      return genos.clone(tgenos, format='sdat', loci=tuple(rows),    packed=False, materialized=False)
    else:
      return genos.clone(tgenos, format='ldat', samples=tuple(rows), packed=False, materialized=False)


  def as_genotriples(self):
    '''
    Return the current genomatrix data as a GenotripleStream.
    Ordered triple streams can be transformed without full materialization,
    though unordered streams do require materialization.  A merge function is
    required for non-unique streams.

    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ ('G', 'G'),   ('A', 'A') ]),
    ...         ('s2', [ ('G', 'T'),   ('T', 'T') ]),
    ...         ('s3', [ ('G', 'G'),   ('A', 'A') ])]
    >>> triples = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_genotriples()
    >>> sorted(triples.loci)
    ['l1', 'l2']
    >>> triples.order
    >>> for row in triples:
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', ('G', 'T'))
    ('s2', 'l2', ('T', 'T'))
    ('s3', 'l1', ('G', 'G'))
    ('s3', 'l2', ('A', 'A'))
    '''
    return build_genotriples_from_genomatrix(self)

  def as_ldat(self, mergefunc=None):
    '''
    Return a genomatrix stream in ldat format.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @return         : ldat genomatrix stream
    @rtype          : GenomatrixStream

    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ ('G', 'G'),   ('A', 'A') ]),
    ...         ('s2', [ ('G', 'T'),   ('T', 'T') ]),
    ...         ('s3', [ ('G', 'G'),   ('A', 'A') ])]
    >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).as_ldat()
    >>> genos.samples
    ('s1', 's2', 's3')
    >>> for row in genos:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
    '''
    genos = self
    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if genos.format != 'ldat':
      genos = genos.transposed()

    return genos

  def as_sdat(self, mergefunc=None):
    '''
    Return a genomatrix stream in sdat format.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
    @type  mergefunc: callable
    @return         : sdat genomatrix stream
    @rtype          : GenomatrixStream

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')]),
    ...         ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])]
    >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples).as_sdat()
    >>> genos.loci
    ('l1', 'l2')
    >>> for row in genos:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])
    '''
    genos = self
    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if genos.format != 'sdat':
      genos = genos.transposed()

    return genos


#######################################################################################


class NonUniqueError(ValueError): pass


def unique_check_genomatrixstream(genos):
  '''
  Check that all row and column labels of a genomatrix are unique.  Raises
  a NonUniqueError if they are not.

  @param rows: genotype matrix data with the first row
               being the column meta-data
  @type rows: sequence

  Non-unique columns:

  >>> genos = GenomatrixStream([],'sdat',loci=['L1','L2','L3','L1'],models=[snp]*4,
  ...                                    genome=Genome(),phenome=Phenome())
  >>> unique_check_genomatrixstream(genos)
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique loci: L1:2

  Non-unique rows:

  >>> loci=('L1','L2')
  >>> rows=[('R1',['AA','AC']),
  ...       ('R1',['AA','AC'])]
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> genos = unique_check_genomatrixstream(genos)
  >>> list(genos)
  Traceback (most recent call last):
       ...
  NonUniqueError: Non-unique row name: R1

  Known unique rows and columns:

  >>> loci=('L1','L2')
  >>> samples=('R1', 'R2')
  >>> rows=[('R1',['AA','AC']),
  ...       ('R2',['AA','AC'])]
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci,samples=samples)
  >>> ugenos = unique_check_genomatrixstream(genos)
  >>> genos is ugenos
  True

  Known columns, unknown but unique rows:

  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> genos = unique_check_genomatrixstream(genos)
  >>> for sample,row in genos:
  ...   print sample,row
  R1 [('A', 'A'), ('A', 'C')]
  R2 [('A', 'A'), ('A', 'C')]
  '''
  assert genos.columns is not None

  if genos.loci is not None:
    dup_loci = [ (k,n) for k,n in tally(genos.loci).iteritems() if n>1 ]
    if dup_loci:
      msg = ','.join( '%s:%d' % kv for kv in dup_loci )
      raise NonUniqueError('Non-unique loci: %s' % msg)

  if genos.samples is not None:
    dup_samples = [ (k,n) for k,n in tally(genos.samples).iteritems() if n>1 ]
    if dup_samples:
      msg = ','.join( '%s:%d' % kv for kv in dup_samples )
      raise NonUniqueError('Non-unique samples: %s' % msg)

  # FASTPATH: Unique samples and loci
  if None not in (genos.samples,genos.loci):
    genos.unique = True
    return genos

  # SLOWPATH: Check rows as they stream past
  def _check():
    drows = set()
    for label,row in genos:
      if label in drows:
        raise NonUniqueError('Non-unique row name: %s' % label)

      drows.add(label)

      yield label,row

  return genos.clone(_check(),materialized=False,unique=True)


#######################################################################################


def _encoding_error(locus,item,model,warn=False):
  '''
  Handle genotype encoding error by either producing an informative exception
  or a warning message.
  '''
  if is_str(item):
    item = 'allele %s' % item
  else:
    item = 'genotype %s' % (','.join(item))

  msg = 'Locus model %s cannot accommodate %s (max_alleles=%d,alleles=%s)' \
                      % (locus,item,model.max_alleles,','.join(model.alleles[1:]))

  if warn:
    sys.stderr.write('[WARNING] %s\n' % msg)
  else:
    raise GenotypeRepresentationError(msg)


def _sample_encoding_error(loci,models,genos,warn=False):
  '''
  Handle sample genotype encoding errors by either producing an informative
  exception or a warning message.
  '''
  for i in xrange(len(genos)):
    model = models[i]
    geno  = genos[i]
    try:
      model[geno]
    except GenotypeLookupError:
      _encoding_error(loci[i],geno,model,warn)
      genos[i] = model[None,None]


def recode_genomatrixstream(genos, genome, warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded with representations
  defined by the supplied genome object.  Locus metadata other than models
  are merged and discrepencies raise errors, If genotype models change, then
  all genotypes are recoded to use the same representation provided the
  models are compatible.

  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : new genomatrixstream with encoding identical to the
                       supplied genome
  @rtype             : GenomatrixStream

  >>> defmodel = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome = Genome()
  >>> genome.set_locus('l1',model=defmodel)
  >>> genome.set_locus('l2',model=defmodel)
  >>> genome.set_locus('l3',model=defmodel)
  >>> genome.set_locus('l4',model=defmodel)

  Test ldat "unknown" models remapped to a default model

  >>> samples = ('s1', 's2', 's3')
  >>> rows = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...         ('l2', [(None, None), (None, None),  (None, None)]),
  ...         ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...         ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples)
  >>> len(genos.models)
  0
  >>> genos = recode_genomatrixstream(genos,genome).materialize()
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> all(model is defmodel for model in genos.models)
  True
  >>> for (label,row),model in izip_exact(genos,genos.models):
  ...   print label,row,model is defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test sdat known models

  >>> samples = ('s1', 's2', 's3')
  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).as_sdat()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,genome).materialize()
  >>> genos.loci
  ('l1', 'l2', 'l3', 'l4')
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> all(model is defmodel for model in genos.models)
  True
  >>> for label,row in genos:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  Test ldat fastpath for no recoding

  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).materialize()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,Genome())
  >>> for (label,row),model in izip_exact(genos,genos.models):
  ...   print label,row,model is not defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test ldat fastpath for known models

  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).materialize()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,genome)
  >>> for (label,row),model in izip_exact(genos,genos.models):
  ...   print label,row,model is defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test sequential recoding to ensure that the same model is used throughout

  >>> samples =          ('s1',         's2',        's3')
  >>> rows1 = [('l1',[ ('G', 'G'),   ('G', 'T'),  ('T', 'T')]),
  ...          ('l2',[ ('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> rows2 = [('l1',[(None, None),  ('T', 'T'),  ('G', 'G')]),
  ...          ('l3',[ ('A', 'A'),  (None, None), ('A', 'T')])]
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples)
  >>> genos2 = GenomatrixStream.from_tuples(rows2,'ldat',samples=samples)
  >>> genome = Genome()
  >>> genos1 = recode_genomatrixstream(genos1, genome).materialize()
  >>> sorted(genome.loci)
  ['l1', 'l2']
  >>> for locus,row in genos1:
  ...   print locus,row
  l1 [('G', 'G'), ('G', 'T'), ('T', 'T')]
  l2 [('A', 'A'), ('T', 'T'), ('A', 'T')]

  >>> genos2 = recode_genomatrixstream(genos2, genome).materialize()
  >>> for locus,row in genos2:
  ...   print locus,row
  l1 [(None, None), ('T', 'T'), ('G', 'G')]
  l3 [('A', 'A'), (None, None), ('A', 'T')]

  >>> sorted(genome.loci)
  ['l1', 'l2', 'l3']
  >>> for locus,model in genos1.model_pairs:
  ...   assert genome.get_model(locus) is model
  >>> for locus,model in genos2.model_pairs:
  ...   assert genome.get_model(locus) is model
  '''
  # Fastpath for null recoding
  if genos.genome is genome:
    return genos

  # Slowpath
  models = []

  # All loci and models are known
  if genos.loci is not None and len(genos.models) == len(genos.loci):
    recode = False
    for i,locus in enumerate(genos.loci):
      # Get the new model or fix the old model
      old_model = genos.models[i]
      old_locus = genos.genome.loci[locus]
      assert old_locus.model is old_model or None in (old_model,old_locus.model)

      if locus not in genome.loci and old_locus.model is not None:
        loc = genome.loci[locus] = old_locus
      else:
        genome.merge_locus(locus, None, old_locus.fixed,    old_locus.chromosome,
                                        old_locus.location, old_locus.strand, warn)

        loc = genome.get_locus(locus)

        if loc.model is None:
          loc.model = old_model

      model = loc.model

      # Check to see if a full recoding or update is necessary
      if model is not old_model:
        recode = True
        try:
          for g in old_model.genotypes[1:]:
            model.add_genotype(g)
        except GenotypeRepresentationError:
          _encoding_error(locus,g,model,warn)

      models.append(model)

    # FASTPATH: No models change, so return with the updated genome
    if not recode:
      assert genos.models == models
      return genos.clone(genos.use_stream(), genome=genome, materialized=False)

  # MEDIUM PATH: Ldat, recoding needed, known models
  if models and genos.format=='ldat':
    def _recode_genomatrixstream():
      n = len(genos.samples)

      descrcache = {}
      packed = genos.packed

      for (locus,row),old_model,model in izip(genos,genos.models,models):
        # Cache the descriptor for this model, since we're likely to see it again
        if packed:
          assert old_model is row.descriptor.models[0]
          descr = descrcache[old_model] = row.descriptor

        # Get or build the new descriptor
        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        # If the model changed, recode by adding all genotypes and packing
        if old_model is not model:
          loc = genome.get_locus(locus)
          if not loc.fixed:
            # Unpack to speed updates and repacking
            row = row[:]
            try:
              for g in set(row):
                model.add_genotype(g)
            except GenotypeRepresentationError:
              _encoding_error(locus,g,model,warn)
              row = None

          row = GenotypeArray(descr,row)

        # Otherwise, only repack if necessary
        elif not packed:
          row = GenotypeArray(descr,row)

        yield locus,row

  # SLOWPATH: Ldat without model information
  elif genos.format=='ldat':
    def _recode_genomatrixstream():
      n = len(genos.samples)

      descrcache = {}
      packed = genos.packed

      for i,(locus,row) in enumerate(genos):
        # Get the new model or fix the old model
        old_model = genos.models[i]
        old_locus = genos.genome.loci[locus]
        assert old_locus.model is old_model or None in (old_model,old_locus.model)

        genome.merge_locus(locus, fixed=old_locus.fixed, chromosome=old_locus.chromosome,
                                  location=old_locus.location, strand=old_locus.strand, warn=warn)

        loc = genome.get_locus(locus)
        if loc.model is None:
          loc.model = old_model

        model = loc.model

        # Cache the descriptor for this model, since we're likely to see it again
        if packed:
          assert old_model is row.descriptor.models[0]
          descr = descrcache[old_model] = row.descriptor

        # Get or build the new descriptor
        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        # If the model changed, recode by adding all genotypes and packing
        if old_model is not model:
          if not loc.fixed:
            # Unpack to speed updates and repacking
            row = row[:]
            try:
              for g in set(row):
                model.add_genotype(g)
            except GenotypeRepresentationError:
              _encoding_error(locus,g,model,warn)
              row = None

          row = GenotypeArray(descr,row)

        # Otherwise, only repack if necessary
        elif not packed:
          row = GenotypeArray(descr,row)

        models.append(model)
        yield locus,row

  # sdat format
  elif genos.format=='sdat':
    assert genome is not genos.genome
    assert genos.loci is not None and len(genos.loci) == len(genos.models) == len(models)

    # Find all models that must be updated
    updates = []
    for i,(locus,old_model,model) in enumerate(izip(genos.loci,genos.models,models)):
      if model is not old_model:
        loc = genome.get_locus(locus)
        if not loc.fixed:
          updates.append( (i,model.add_genotype) )

    # FASTERPATH: If all models are fixed, recoding is straightforward
    if not updates:
      def _recode_genomatrixstream():
        descr = GenotypeArrayDescriptor(models)
        for sample,row in genos:
          try:
            # Recode and yield new row
            row = GenotypeArray(descr,row)
          except GenotypeLookupError:
            _sample_encoding_error(genos.loci,models,row,warn)

          yield sample,row

    # SLOWPATH: Otherwise, recode by adding genotypes from all changed
    # models and packing.
    else:
      def _recode_genomatrixstream():
        descr = GenotypeArrayDescriptor(models)
        for sample,row in genos:
          # Unpack row to speed access -- both updates and GenotypeArray will
          # need to unpack
          row = row[:]

          # Try to yield updated genotype array, hoping that all alleles are represented
          try:
            yield sample,GenotypeArray(descr,row)

          except GenotypeRepresentationError:
            # Update all changed models to ensure they contain the needed alleles
            for i,add in updates:
              try:
                add(row[i])
              except GenotypeRepresentationError:
                _encoding_error(genos.loci[i],row[i],models[i],warn)
                row[i] = models[i][None,None]

            # Recode and yield new row
            yield sample,GenotypeArray(descr,row)

  else:
    raise ValueError('Uknown format')

  return genos.clone(_recode_genomatrixstream(),models=models,genome=genome,packed=True,materialized=False)


def sdat_model_lookahead_from_strings(loci,genos,genome,genorepr,min_unknown=10,max_lookahead=50,warn=False):
  '''
  Lookahead in an sdat genotype stream to determine alleles and fixed
  models.  This is a major optimization that can usually avoid the majority
  of the overhead associated with creating len(loci) default models
  '''
  # Do not attempt to lookahead if a fixed model is used for unknown loci
  if genome.default_fixed:
    return genos

  # Do not attempt to lookahead if only a small number of loci are unknown
  unknown = sum(1 for locus in loci if genome.get_locus(locus).model is None)
  if unknown < min_unknown:
    return genos

  # Track the alleles seen
  alleles_seen = [ [] for i in range(len(loci)) ]

  # Start looking ahead up to max_lookahead rows
  lookahead_rows = []
  for sample,row in genos:
    lookahead_rows.append( (sample,row) )

    row = genorepr.from_strings(row)

    for g,seen in izip(row,alleles_seen):
      seen.append(g[0])
      seen.append(g[1])

    if len(lookahead_rows) == max_lookahead:
      break

  # If there are rows in the lookahead buffer
  if not lookahead_rows:
    return genos

  max_alleles = genome.max_alleles
  modelcache = {}

  try:
    # Review the alleles seen at each locus
    for locus,seen in izip(loci,alleles_seen):
      loc = genome.get_locus(locus)
      model = loc.model

      # If the model is unknown, then check the alleles
      if model is None:
        seen = set(seen)
        seen.discard(None)

        # Create or reuse a fixed model if all alleles have been seen
        # FIXME: add support for hemizygote models
        if len(seen) == max_alleles:
          seen  = tuple(sorted(seen))
          model = modelcache.get(seen)
          if model is None:
            model = modelcache[seen] = model_from_alleles(seen, max_alleles=max_alleles)
          loc.model = model
          loc.fixed = True
          continue

        # Otherwise create an empty default model
        model = genome.get_model(locus)

      # Populate the observed alleles when a fixed model cannot be used
      for allele in seen:
        model.add_allele(allele)

  except GenotypeRepresentationError:
    _encoding_error(locus,allele,model,warn)

  return chain(lookahead_rows,genos)


def sdat_model_lookahead_from_tuples(loci,genos,genome,min_unknown=10,max_lookahead=25,warn=False):
  '''
  Lookahead in an sdat genotype stream to determine alleles and fixed
  models.  This is a major optimization that can usually avoid the majority
  of the overhead associated with creating len(loci) default models
  '''
  # Do not attempt to lookahead if a fixed model is used for unknown loci
  if genome.default_fixed:
    return genos

  # Do not attempt to lookahead if only a small number of loci are unknown
  unknown = sum(1 for locus in loci if genome.get_locus(locus).model is None)
  if unknown < min_unknown:
    return genos

  # Track the alleles seen
  alleles_seen = [ [] for i in range(len(loci)) ]

  # Start looking ahead up to max_lookahead rows
  lookahead_rows = []
  for sample,row in genos:
    lookahead_rows.append( (sample,row) )

    for g,seen in izip(row,alleles_seen):
      seen.append(g[0])
      seen.append(g[1])

    if len(lookahead_rows) == max_lookahead:
      break

  # If there are rows in the lookahead buffer
  if not lookahead_rows:
    return genos

  max_alleles = genome.max_alleles
  modelcache = {}

  try:
    # Review the alleles seen at each locus
    for locus,seen in izip(loci,alleles_seen):
      loc = genome.get_locus(locus)
      model = loc.model

      # If the model is unknown, then check the alleles
      if model is None:
        seen = set(seen)
        seen.discard(None)

        # Create or reuse a fixed model if all alleles have been seen
        # FIXME: add support for hemizygote models
        if len(seen) == max_alleles:
          seen  = tuple(sorted(seen))
          model = modelcache.get(seen)
          if model is None:
            model = modelcache[seen] = model_from_alleles(seen, max_alleles=max_alleles)
          loc.model = model
          loc.fixed = True
          continue

        # Otherwise create an empty default model
        model = genome.get_model(locus)

      # Populate the observed alleles when a fixed model cannot be used
      for allele in seen:
        model.add_allele(allele)

  except GenotypeRepresentationError:
    _encoding_error(locus,allele,model,warn)

  return chain(lookahead_rows,genos)


def encode_genomatrixstream_from_tuples(columns, genos, format, genome=None,
                                                 unique=False, warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       format: format of input genomatrix, either 'ldat' or 'sdat'
  @type        format: str
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
  @type        unique: bool
  @return            : tuple of columns and a genomatrix generator in packed format
  @rtype             : 2-tuple of list of str and genomatrix generator

  >>> defmodel = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome  = Genome(default_model=defmodel, default_fixed=True)

  With genome:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
  >>> samples
  ('s1', 's2', 's3')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True
  l2 [(None, None), (None, None), (None, None)] True
  l3 [('A', 'A'), (None, None), (None, None)] True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')]),
  ...          ('s2', [(None, None), (None, None), (None, None), (None, None)]),
  ...          ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])]
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  No genome:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat')
  >>> samples
  ('s1', 's2', 's3')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True
  l2 [(None, None), (None, None), (None, None)] True
  l3 [('A', 'A'), (None, None), (None, None)] True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')]),
  ...          ('s2', [(None, None), (None, None), (None, None), (None, None)]),
  ...          ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])]
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  See if we can provide a subtle cache bug when models are cached too
  aggressively for non-unique loci:

  >>> samples = ('s1',)
  >>> genos = [('l1', [ ('A','A') ]),
  ...          ('l2', [ ('A','A') ]),
  ...          ('l1', [ ('A','T') ]),
  ...          ('l2', [ ('A','G') ])]
  >>> genome = Genome()
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
  >>> new_rows = list(new_rows)
  >>> samples
  ('s1',)
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A')] True
  l2 [('A', 'A')] True
  l1 [('A', 'T')] True
  l2 [('A', 'G')] True
  '''
  if genome is None:
    genome = Genome()

  models = []

  if format=='ldat':
    def _encode():
      n = len(columns)
      m = genome.max_alleles+1

      # FIXME: add support for hemizygote models
      modelcache = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                               if locus.model is not None and len(locus.model.alleles)==m )
      descrcache = {}

      def _genokey(genos):
        gset = set(a for g in set(genos) for a in g)
        gset.discard(None)
        return tuple(sorted(gset))

      for locus,row in genos:
        key = None

        loc = genome.get_locus(locus)

        if loc.model is None:
          key = _genokey(row)
          cachable = genome.default_model is None and len(key) == genome.max_alleles

          if cachable:
            # Aggressively reuse models with fully compatible alleles
            loc.model = modelcache.get(key)

          if loc.model is None:
            genome.get_locus_model(locus)

          if cachable:
            modelcache[key] = loc.model

        model = loc.model

        assert model is not None

        if not loc.fixed:
          key = key or _genokey(row)
          try:
            for a in key:
              model.add_allele(a)
          except GenotypeRepresentationError:
            _encoding_error(locus,a,model,warn)
            row = None

        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        assert descr.models[0] is model
        models.append(model)
        yield locus,GenotypeArray(descr,row)

  elif format=='sdat':

    genos = sdat_model_lookahead_from_tuples(columns,genos,genome)

    updates = []

    for i,locus in enumerate(columns):
      loc = genome.get_locus_model(locus)
      models.append(loc.model)
      if not loc.fixed:
        updates.append( (i,loc.model.add_genotype) )

    def _encode():
      descr = GenotypeArrayDescriptor(models)

      if not updates:
        for sample,row in genos:
          yield sample,GenotypeArray(descr,row)
      else:
        for sample,row in genos:
          for i,add in updates:
            try:
              add(row[i])
            except GenotypeRepresentationError:
              _encoding_error(columns[i],row[i],models[i],warn)
              row[i] = models[i][None,None]

          yield sample,GenotypeArray(descr,row)
  else:
    raise ValueError('Uknown format')

  return columns,models,genome,_encode()


def encode_genomatrixstream_from_strings(columns,genos,format,genorepr,genome=None,
                                                 unique=False,warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       format: format of input genomatrix, either 'ldat' or 'sdat'
  @type        format: str
  @param     genorepr: internal representation of genotypes for the input/output
  @type      genorepr: UnphasedMarkerRepresentation or similar object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
  @type        unique: bool
  @return            : tuple of columns and a genomatrix generator in packed format
  @rtype             : 2-tuple of list of str and genomatrix generator

  >>> from reprs import snp
  >>> defmodel  = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome = Genome(default_model=defmodel, default_fixed=True)

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', ['AA', '  ', 'AA', 'GT']),
  ...          ('s2', ['  ', '',   '  ', '  ']),
  ...          ('s3', ['GG', '  ', '  ', 'TT'])]
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp,genome)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in new_rows:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', ['AA', '  ', 'AA', 'GT']),
  ...          ('s2', ['  ', '',   '  ', '  ']),
  ...          ('s3', ['GG', '  ', '  ', 'TT'])]
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in new_rows:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])

  See if we can provide a subtle cache bug when models are cached too
  aggressively for non-unique loci:

  >>> samples = ('s1',)
  >>> genos = [('l1', [ ('AA') ]),
  ...          ('l2', [ ('AA') ]),
  ...          ('l1', [ ('AT') ]),
  ...          ('l2', [ ('AG') ])]
  >>> genome = Genome()
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
  >>> new_rows = list(new_rows)
  >>> samples
  ('s1',)
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A')] True
  l2 [('A', 'A')] True
  l1 [('A', 'T')] True
  l2 [('A', 'G')] True
  '''
  if genome is None:
    genome = Genome()

  models = []

  if format=='ldat':
    def _encode():
      n = len(columns)
      m = genome.max_alleles+1

      # FIXME: add support for hemizygote models
      modelcache   = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                               if locus.model is not None and len(locus.model.alleles)==m )
      descrcache   = {}
      strcache     = {}
      to_string    = genorepr.to_string
      from_strings = genorepr.from_strings

      def _genokey(genos):
        gset = set(a for g in set(from_strings(genos)) for a in g)
        gset.discard(None)
        return tuple(sorted(gset))

      for locus,row in genos:
        key = None
        loc = genome.get_locus(locus)

        if loc.model is None:
          key = _genokey(row)
          cachable = genome.default_model is None and len(key) == genome.max_alleles

          if cachable:
            # Aggressively reuse models with fully compatible alleles
            loc.model = modelcache.get(key)

          if loc.model is None:
            genome.get_locus_model(locus)

          if cachable:
            modelcache[key] = loc.model

        model = loc.model

        if not loc.fixed:
          key = key or _genokey(row)
          try:
            for a in key:
              model.add_allele(a)
          except GenotypeRepresentationError:
            _encoding_error(locus,a,model,warn)
            row = ['']*len(row)

        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        cache = strcache.get(model)
        if cache is None:
          cache = strcache[model] = dict( (to_string(g),g) for g in model.genotypes )
          for g in genorepr.missing_geno_strs:
            cache[g] = model[None,None]

        models.append(model)

        try:
          yield locus,GenotypeArray(descr,imap(getitem, repeat(cache), row))
        except KeyError:
          gset = set(row)
          try:
            cache.update( (g,model[r]) for g,r in izip(gset,from_strings(gset)) )
            yield locus,GenotypeArray(descr,imap(getitem, repeat(cache), row))
          except KeyError,g:
            _encoding_error(locus,g,model,warn)
            yield locus,GenotypeArray(descr)

  elif format=='sdat':

    genos = sdat_model_lookahead_from_strings(columns,genos,genome,genorepr)

    n = len(columns)
    updates   = []
    cachemap  = {}
    cachelist = []

    to_string    = genorepr.to_string
    from_strings = genorepr.from_strings

    for i,locus in enumerate(columns):
      loc = genome.get_locus_model(locus)
      model = loc.model
      models.append(model)
      if loc.fixed:
        update = model.get_genotype
      else:
        update = model.add_genotype

      cache = cachemap.get(model)
      if cache is None:
        cache = cachemap[model] = dict( (to_string(g),g) for g in model.genotypes )
        cache.update( (genorepr.to_string_from_alleles((g[1],g[0])),g) for g in model.genotypes )
        for g in genorepr.missing_geno_strs:
          cache[g] = model[None,None]

      cachelist.append(cache)
      updates.append( (locus,update,cache) )

    def _encode():
      repr   = genorepr.from_string
      descr  = GenotypeArrayDescriptor(models)
      errors = (KeyError,ValueError,GenotypeRepresentationError)

      for sample,row in genos:
        try:
          row = GenotypeArray(descr,imap(getitem, cachelist, row) )
        except errors:
          geno_tuples = from_strings(row)

          for (locus,update,cache),gstr,g in izip(updates,row,geno_tuples):
            # Aggressively form homozygote genotypes and cache them.  Thus
            # we only see cache misses when we encounter previously
            # unobserved alleles or when genotype formatting is off.
            g1,g2 = g
            if g1!=g2:
              gh = g1,g1
              try:
                cache[genorepr.to_string_from_alleles(gh)] = update(gh)
              except errors:
                pass

              gh = g2,g2
              try:
                cache[genorepr.to_string_from_alleles(gh)] = update(gh)
              except errors:
                pass

            try:
              cache[gstr] = update(g)
            except errors:
              pass

          try:
            row = GenotypeArray(descr,imap(getitem, cachelist, row) )
          except errors:
            _sample_encoding_error(columns,models,geno_tuples,warn)

        yield sample,row
  else:
    raise ValueError('Uknown format')

  return columns,models,genome,_encode()


def recode_genotriples(triples,genome):
  '''
  Returns a new genotriples with the genotypes encoded to a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation). e.g.
                       ('s1','l1','AA'),...
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> triples = [('s3', 'l1', ('G', 'G')),('s3', 'l2', ('A', 'A')),
  ...            ('s2', 'l3', ('G', 'T')),('s1', 'l1', ('T', 'T')),
  ...            ('s1', 'l1', ('G', 'G')),('s2', 'l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> genome = Genome()
  >>> for row in recode_genotriples(triples,genome):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  >>> sorted(genome.loci)
  ['l1', 'l2', 'l3']
  '''
  def _recode():
    updates = {}
    try:
      for sample,locus,geno in triples:
        ud = updates.get(locus)

        if ud is None:
          old_model = geno.model
          old_locus = triples.genome.loci[locus]
          assert old_model is triples.genome.loci[locus].model
          assert old_locus.model is old_model

          genome.merge_locus(locus, fixed=old_locus.fixed, chromosome=old_locus.chromosome,
                                          location=old_locus.location)

          loc = genome.get_locus(locus)
          if loc.model is None:
            loc.model = old_model

          if old_model is loc.model or loc.fixed:
            ud = loc.model.get_genotype
          else:
            ud = loc.model.add_genotype

          updates[locus] = ud

        yield sample,locus,ud(geno)

    except GenotypeRepresentationError:
      model = triples.genome.get_model(locus)
      _encoding_error(locus,geno,model)

  return triples.clone(_recode(),genome=genome,materialized=False)


def encode_genotriples_from_tuples(triples,genome):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation)
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> triples = [('s3', 'l1', ('G', 'G')),('s3', 'l2', ('A', 'A')),
  ...            ('s2', 'l3', ('G', 'T')),('s1', 'l1', ('T', 'T')),
  ...            ('s1', 'l1', ('G', 'G')),('s2', 'l2', ('A', 'A'))]
  >>> for row in encode_genotriples_from_tuples(triples,Genome()):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  try:
    for sample,locus,geno in triples:
      loc = genome.get_locus_model(locus)

      if loc.fixed:
        geno = loc.model.get_genotype(geno)
      else:
        geno = loc.model.add_genotype(geno)

      yield sample,locus,geno
  except GenotypeRepresentationError:
    _encoding_error(locus,geno,loc.model)


def encode_genotriples_from_strings(triples,genorepr,genome):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation)
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> from reprs import snp
  >>> triples = [('s3', 'l1', 'GG'),('s3', 'l2', 'AA'),
  ...            ('s2', 'l3', 'GT'),('s1', 'l1', 'TT'),
  ...            ('s1', 'l1', 'GG'),('s2', 'l2', 'AA')]
  >>> for row in encode_genotriples_from_strings(triples,snp,Genome()):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  local_repr = genorepr.from_string
  if genome is None:
    genome = Genome()

  updates = {}

  try:
    for sample,locus,geno in triples:
      ud = updates.get(locus)

      if ud is None:
        loc = genome.get_locus_model(locus)
        if loc.fixed:
          ud = lambda g,get=loc.model.get_genotype: get(local_repr(g))
        else:
          ud = lambda g,add=loc.model.add_genotype: add(local_repr(g))

        updates[locus] = ud

      yield sample,locus,ud(geno)

  except GenotypeRepresentationError:
    _encoding_error(locus,geno,loc.model)


#######################################################################################


def sort_genotriples(triples,order,locusorder=None,sampleorder=None,maxincore=1000000):
  '''
  Sort genotriples by the specified order.  When

  order='sample' -> sort by sample,locus
  order='locus'  -> sort by locus,sample

  Samples and loci are ordered by the sequence found in locusorder and
  sampleorder, if specified, and the remaining loci or samples are sorted in
  lexicographical order.

  This current version fully materializes the input triple-stream and
  performs the sort in core memory.  Future versions will honor the incore
  parameter and spill very large streams to disk and perform a multi-phase
  offline merge sort.

  @param     triples: sequence of genotriples(str,str,genotype representation)
  @type      triples: sequence
  @param       order: sort order, either 'samples', 'locus'
  @type        order: str
  @param  locusorder: ordered list of loci (optional)
  @type   locusorder: sequence
  @param sampleorder: ordered list of samples (optional)
  @type  sampleorder: sequence
  @param   maxincore: maximum number of triples to process in core (not currently honored)
  @type    maxincore: int or None
  @return           : samples, loci, sorted genotriples
  @rtype            : tuple of list, list, genotriple sequence

  >>> triples = [('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A')),
  ...            ('s2','l2', ('A', 'T')),('s1','l1', ('T', 'T')),
  ...            ('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples).materialize()

  >>> striples = sort_genotriples(triples,order='sample',sampleorder=['s1','s2','s3'],
  ...                                                     locusorder=['l1','l2','l3'])
  >>> sorted(striples.samples)
  ['s1', 's2', 's3']
  >>> sorted(striples.loci)
  ['l1', 'l2']
  >>> list(striples)
  [('s1', 'l1', ('G', 'G')), ('s1', 'l1', ('T', 'T')), ('s2', 'l2', ('A', 'A')), ('s2', 'l2', ('A', 'T')), ('s3', 'l1', ('G', 'G')), ('s3', 'l2', ('A', 'A'))]

  >>> striples = sort_genotriples(triples,order='sample',locusorder=['l3','l2','l1'])
  >>> sorted(striples.samples)
  ['s1', 's2', 's3']
  >>> sorted(striples.loci)
  ['l1', 'l2']
  >>> list(striples)
  [('s1', 'l1', ('G', 'G')), ('s1', 'l1', ('T', 'T')), ('s2', 'l2', ('A', 'A')), ('s2', 'l2', ('A', 'T')), ('s3', 'l2', ('A', 'A')), ('s3', 'l1', ('G', 'G'))]

  >>> striples = sort_genotriples(triples,order='sample',sampleorder=['s3','s2','s1'])
  >>> sorted(striples.samples)
  ['s1', 's2', 's3']
  >>> sorted(striples.loci)
  ['l1', 'l2']
  >>> list(striples)
  [('s3', 'l1', ('G', 'G')), ('s3', 'l2', ('A', 'A')), ('s2', 'l2', ('A', 'A')), ('s2', 'l2', ('A', 'T')), ('s1', 'l1', ('G', 'G')), ('s1', 'l1', ('T', 'T'))]
  '''
  def makeorderdict(o):
    od = {}
    for x in o or []:
      od.setdefault(x, (len(od),x) )
    return od

  def keyfunctor(i1,orderdict1,i2,orderdict2):
    n1 = len(orderdict1)
    n2 = len(orderdict2)
    def _keyfunc(r):
      r1,r2 = r[i1],r[i2]
      k1 = orderdict1.get(r1, (n1,r1))
      k2 = orderdict2.get(r2, (n2,r2))
      return k1,k2,r[2]
    return _keyfunc

  locusorder  = makeorderdict(locusorder)
  sampleorder = makeorderdict(sampleorder)

  if order == 'sample':
    keyfunc = keyfunctor(0, sampleorder, 1, locusorder)
  elif order == 'locus':
    keyfunc = keyfunctor(1, locusorder,  0, sampleorder)
  else:
    raise ValueError("Unknown ordering specified: '%s'" % order)

  # In memory sort should eventually be augmented by an offline multiphase
  # merge sort for very large sets.
  striples = sorted(triples,key=keyfunc)

  # Extract sample and locus sets since this can be done quickly with materialized lists
  samples = set(imap(itemgetter(0), striples))
  loci    = set(imap(itemgetter(1), striples))

  return triples.clone(striples,samples=samples,loci=loci,materialized=True,order=order)


def combine_unsorted_genotriple_list(triplelist):
  '''
  Combine multiple genotriples into one unsorted triple stream

  @param triplelist: list of genotriples
  @type  triplelist: list
  @return          : combined genotriple stream
  @rtype           : sequence of sample, locus, and genotype

  >>> trip1 = [('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A')),
  ...          ('s2','l2', ('A', 'T')),('s1','l1', ('T', 'T'))]
  >>> trip2 = [('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
  >>> trip1 = GenotripleStream.from_tuples(trip1)
  >>> trip2 = GenotripleStream.from_tuples(trip2)
  >>> triplelist = [trip1,trip2]
  >>> combined = combine_unsorted_genotriple_list(triplelist)
  >>> for row in combined:
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l2', ('A', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  # Make sure that input list is materialized and can be iterated
  if not isinstance(triplelist, list):
    triplelist = list(triplelist)

  if not triplelist:
    raise TypeError('cannot combine a single genotriple stream')
  elif len(triplelist) == 1:
    return triplelist[0]

  # Normalize all triples to the same models and phenome
  genome     = triplelist[0].genome
  triplelist = [ triples.transformed(recode_models=genome) for triples in triplelist ]
  phenome    = merge_phenome_list([triples.phenome for triples in triplelist ])

  # Extract parts of all of the triples
  samples = [ triples.samples  for triples in triplelist ]
  loci    = [ triples.loci     for triples in triplelist ]
  order   = [ triples.order    for triples in triplelist ]
  unique  = [ triples.unique   for triples in triplelist ]

  # If any of the triples have unknown samples or loci, mark the results as unknown
  if None in samples:
    samples = None
  if None in loci:
    loci = None

  # Otherwise, count the samples and loci to determine if they are unique
  # and consolidate the sets
  if samples is not None and loci is not None:
    samples = tally(chain(*samples))
    loci    = tally(chain(*loci))
    unique  = all( v<=1 for v in chain(samples.itervalues(),loci.itervalues()) )
    samples = set(samples)
    loci    = set(loci)
  else:
    unique  = False

  triples = chain(*triplelist)

  # Return a new baby triplestream object
  return GenotripleStream(triples,samples=samples,loci=loci,genome=genome,phenome=phenome,
                                  unique=unique)


def combine_sorted_genotriple_list(triplelist):
  '''
  Combine multiple sorted genotriples into one sorted triple stream

  @param triplelist: list of genotriples
  @type  triplelist: list
  @return          : combined genotriple stream
  @rtype           : sequence of sample, locus, and genotype

  >>> trip1 = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
  ...          ('s2','l2', ('A', 'T')),('s3','l1', ('T', 'T'))]
  >>> trip2 = [('s2','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
  >>> trip1 = GenotripleStream.from_tuples(trip1,order='sample')
  >>> trip2 = GenotripleStream.from_tuples(trip2,order='sample')
  >>> triplelist = [trip1,trip2]
  >>> combined = combine_sorted_genotriple_list(triplelist)
  >>> for row in combined:
  ...   print row
  ('s1', 'l1', ('G', 'G'))
  ('s1', 'l2', ('A', 'A'))
  ('s2', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'T'))
  ('s3', 'l1', ('T', 'T'))
  ('s3', 'l2', ('A', 'A'))
  >>> combined.order
  'sample'
  '''
  # Make sure that input list is materialized and can be iterated
  if not isinstance(triplelist, list):
    triplelist = list(triplelist)

  if not triplelist:
    raise TypeError('cannot combine a single genotriple stream')
  elif len(triplelist) == 1:
    return triplelist[0]

  # Normalize all triples to the same models and phenome
  genome     = triplelist[0].genome
  triplelist = [ triples.transformed(recode_models=genome) for triples in triplelist ]
  phenome    = merge_phenome_list([triples.phenome for triples in triplelist ])

  # Extract parts of all of the triples
  samples = [ triples.samples  for triples in triplelist ]
  loci    = [ triples.loci     for triples in triplelist ]
  order   = [ triples.order    for triples in triplelist ]
  unique  = [ triples.unique   for triples in triplelist ]

  # If any of the triples have unknown samples or loci, mark the results as unknown
  if None in samples:
    samples = None
  if None in loci:
    loci = None

  # Otherwise, count the samples and loci to determine if they are unique
  # and consolidate the sets
  if samples is not None and loci is not None:
    samples = tally(chain(*samples))
    loci    = tally(chain(*loci))
    unique  = all( v<=1 for v in chain(samples.itervalues(),loci.itervalues()) )
  else:
    unique  = False

  # Check that everything is in the same order
  order = set(order)

  if len(order) != 1:
    raise ValueError('Cannot merge triplestreams in disparate orders')

  order = order.pop()

  if order == 'sample':
    key = itemgetter(0,1)
  elif order == 'locus':
    key = itemgetter(1,0)
  else:
    raise ValueError('Invalid order specified')

  # Perform merge
  triples = imerge(triplelist,key=key)

  # Return a new baby triplestream object
  return GenotripleStream(triples,samples=samples,loci=loci,genome=genome,phenome=phenome,
                                  order=order,unique=unique)


#######################################################################################


def merge_sorted_genotriples(triples,mergefunc):
  '''
  Merge genotypes from a genotriple list sorted by both sample and locus
  (although it does not matter which is primary and secondary key).  The
  supplied merge function is used to produce a consensus genotype for each
  sample and locus.

  @param   triples: sequence of sorted genotriples(str,str,genotype representation)
  @type    triples: sequence
  @param mergefunc: function to merge multiple genotypes into a consensus genotype.
  @type  mergefunc: callable
  @return         : sorted, merged genotriple stream
  @rtype          : sequence

  >>> from reprs import snp
  >>> triples = [('l1','s1','AA'),('l1','s1','  '),('l1','s2','AB'),('l2','s1','AA'),
  ...            ('l2','s1','AA'),('l3','s1','BB'),('l3','s1','BB'),('l3','s1','AB')]
  >>> triples = GenotripleStream.from_strings(triples,snp)
  >>> list(merge_sorted_genotriples(triples,VoteMerger()))
  [('l1', 's1', ('A', 'A')), ('l1', 's2', ('A', 'B')), ('l2', 's1', ('A', 'A')), ('l3', 's1', (None, None))]
  '''
  def _merge():
    get2 = itemgetter(2)
    genome = triples.genome
    merge = mergefunc.merge_geno
    for (sample,locus),trips in groupby(triples, itemgetter(0,1)):
      genos = map(get2,trips)
      assert len(genos) < 2 or all(genos[0].model is g.model for g in genos)
      geno  = merge(sample,locus,genome.get_model(locus),genos)
      yield sample,locus,geno

  return triples.clone(_merge(),unique=True,materialized=False)


def merge_genomatrixstream(genos, mergefunc):
  '''
  Merge genotypes with identical row or column labels from a genomatrix
  stream using the specified merge function.  The algorithm used requires a
  full materialization of the genomatrix stream, since row labels must be
  known.  Results are in an unpacked list representation that may need to be
  repacked.

  @param      genos: genomatrix stream
  @type       genos: sequence
  @param  mergefunc: function to merge multiple genotypes into a consensus genotype.
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  >>> loci = ('l1','l2','l3')
  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'G')]),
  ...         ('s2',[(None, None), (None, None), (None, None)]),
  ...         ('s3',[ ('A', 'A'),  (None, None), (None, None)]),
  ...         ('s4',[ ('A', 'T'),  (None, None),  ('T', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False)
  >>> merger= VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('G', 'G')])
  ('s2', [(None, None), (None, None), (None, None)])
  ('s3', [('A', 'A'), (None, None), (None, None)])
  ('s4', [('A', 'T'), (None, None), ('T', 'T')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 1]), ('s2', [0, 0, 0, 0, 3]), ('s3', [1, 0, 0, 0, 2]), ('s4', [2, 0, 0, 0, 1])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 1]), ('l2', [0, 0, 0, 0, 4]), ('l3', [2, 0, 0, 0, 2])]

  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False).as_ldat()
  >>> genos.loci = None
  >>> merger= VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.samples
  ('s1', 's2', 's3', 's4')
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> len(genos.models)
  0
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('A', 'A'), ('A', 'T')])
  ('l2', [(None, None), (None, None), (None, None), (None, None)])
  ('l3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])
  >>> len(genos.models)
  3
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 1]), ('s2', [0, 0, 0, 0, 3]), ('s3', [1, 0, 0, 0, 2]), ('s4', [2, 0, 0, 0, 1])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 1]), ('l2', [0, 0, 0, 0, 4]), ('l3', [2, 0, 0, 0, 2])]

  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'T') ]),
  ...         ('s2',[(None, None), ('A', 'C'),   (None, None)]),
  ...         ('s1',[ ('A', 'A'),  ('A', 'A'),   (None, None)]),
  ...         ('s1',[ ('A', 'T'), (None, None),  ('G', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False)
  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), ('A', 'C'), (None, None)])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [1, 1, 0, 1, 0]), ('s2', [1, 0, 0, 0, 2])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [0, 0, 0, 1, 1]), ('l2', [2, 0, 0, 0, 0]), ('l3', [0, 1, 0, 0, 1])]

  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False).as_ldat()
  >>> merger= VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.samples
  ('s1', 's2')
  >>> for row in genos:
  ...   print row
  ('l1', [(None, None), (None, None)])
  ('l2', [('A', 'A'), ('A', 'C')])
  ('l3', [('G', 'T'), (None, None)])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [1, 1, 0, 1, 0]), ('s2', [1, 0, 0, 0, 2])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [0, 0, 0, 1, 1]), ('l2', [2, 0, 0, 0, 0]), ('l3', [0, 1, 0, 0, 1])]

  >>> loci = ('l1','l2','l1')
  >>> rows = [('s1',[(None, None), (None, None),  ('C', 'T')]),
  ...         ('s2',[(None, None),  ('A', 'G'),   ('T', 'T')]),
  ...         ('s1',[ ('C', 'C'),   ('A', 'A'),  (None, None)]),
  ...         ('s1',[ ('C', 'T'),  (None, None),  ('C', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False)
  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), ('A', 'A')])
  ('s2', [('T', 'T'), ('A', 'G')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [1, 0, 0, 1, 0]), ('s2', [2, 0, 0, 0, 0])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [1, 0, 0, 1, 0]), ('l2', [2, 0, 0, 0, 0])]

  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci,unique=False).as_ldat()
  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.samples
  ('s1', 's2')
  >>> for row in genos:
  ...   print row
  ('l1', [(None, None), ('T', 'T')])
  ('l2', [('A', 'A'), ('A', 'G')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [1, 0, 0, 1, 0]), ('s2', [2, 0, 0, 0, 0])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [1, 0, 0, 1, 0]), ('l2', [2, 0, 0, 0, 0])]
  '''
  assert mergefunc is not None

  # FIXME: Does not update stats
  if genos.unique:
    return genos

  merge_rows = genos.rows    is None or len(genos.rows)    != len(set(genos.rows))
  merge_cols = genos.columns is None or len(genos.columns) != len(set(genos.columns))

  # FIXME: Does not update stats
  if not merge_rows and not merge_cols:
    genos.uniqe = True
    return genos

  new_rows      = {}
  new_columns   = {}
  column_dest   = []
  merge_columns = [],[]

  # We must materialize all streams, so ensure they are packed
  genos = genos.transformed(repack=True)

  # Form mappings from old schemas to new
  column_dest = [ new_columns.setdefault(c,len(new_columns)) for c in genos.columns ]

  # Invert row and column dictionaries to recover insertion orderings
  new_columns   = tuple(imap(itemgetter(0),sorted(new_columns.iteritems(),key=itemgetter(1))))
  merge_columns = (new_columns != genos.columns)

  # Collect row labels and materialize all rows for later merging
  if merge_rows:
    rowdict = defaultdict(list)
    for i,(row_label,row) in enumerate(genos):
      new_rows.setdefault(row_label,len(new_rows))
      rowdict[row_label].append(row)

    new_rows = tuple(imap(itemgetter(0),sorted(new_rows.iteritems(),   key=itemgetter(1))))
    rowdata = ( (label,rowdict.pop(label)) for label in new_rows )

  else:
    rowdata = ( (label,[row]) for label,row in genos.use_stream() )
    new_rows = genos.rows

  if genos.format=='ldat':
    models = []
    samples,loci = new_columns,new_rows
  else:
    models = [ genos.genome.get_model(column) for column in new_columns ]
    samples,loci = new_rows,new_columns

  def _merger():
    n = len(new_columns)

    # Fully general merge over duplicate rows and columns
    for label,rows in rowdata:
      with gcdisabled():
        # All columns are extracted from rows of schema i
        if len(rows) == 1:
          new_row = tuple(rows[0])
        else:
          new_row = pick_columns(rows)

        if merge_columns:
          new_row = place_list([None]*n, new_row, column_dest)

        # Merge genotypes
        if genos.format=='ldat':
          model   = genos.genome.get_model(label)
          new_row = mergefunc.merge_locus(samples, label, model, new_row)
          models.append(model)
        else:
          new_row = mergefunc.merge_sample(label, loci, models, new_row)

      # Yield new row
      yield label,new_row

  return genos.clone(_merger(),samples=samples,loci=loci,models=models,
                               packed=False,materialized=False,unique=True)


def merge_genomatrixstream_list(genos, mergefunc):
  '''
  Take a sequence of genotype matrix streams and merge all genotypes
  identical row or column labels into a single genotype matrix stream using
  the specified merge function.

  All input matricies must meet the following requirements:
    1) share the same orientation (either sdat or ldat)
    2) share the same genotype representation

  Unless all row labels are known and all columns are homogeneous, the
  current implementation requires a full materialization of all input
  genomatrix streams.

  Depending on the code path chosen, results may be in an unpacked list
  representation that may need to be repacked.

  Pending optimizations yet to be implemented:

    * The requirement that columns be homogeneous in order to avoid full
      materialization can be relaxed.  This will result in significant
      performance improvements.

    * the case of disjoint input columns with known and unknown rows can be
      optimized to varying degrees

  @param      genos: genomatrix streams
  @type       genos: sequence
  @param  mergefunc: function to merge multiple genotypes into a consensus genotype.
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  Test slow-path for heterogeneous schema:

  >>> samples1 =        ('s1',       's2',       's3')
  >>> rows1 = [('l1',[('G', 'G'),   ('A', 'A'),  ('A', 'G')]),
  ...          ('l2',[('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples1).materialize()

  >>> samples2 =         ('s1',       's3',        's4')
  >>> rows2 = [('l1',[(None, None), ('G', 'G'),  ('A', 'G')]),
  ...          ('l3',[('A', 'A'),  (None, None), ('A', 'T')])]
  >>> genos2 = GenomatrixStream.from_tuples(rows2,'ldat',samples=samples2).materialize()

  >>> genos = [genos1,genos2,genos1]
  >>> genos = merge_genomatrixstream_list(genos,VoteMerger())
  >>> genos.samples
  ('s1', 's2', 's3', 's4')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), ('A', 'A'), (None, None), ('A', 'G')])
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T'), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None), ('A', 'T')])

  >>> genos = [ g.as_sdat() for g in [genos1,genos2,genos1] ]
  >>> genos = merge_genomatrixstream_list(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A'), ('A', 'A')])
  ('s2', [('A', 'A'), ('T', 'T'), (None, None)])
  ('s3', [(None, None), ('A', 'T'), (None, None)])
  ('s4', [('A', 'G'), (None, None), ('A', 'T')])

  Test fast-path for homogeneous schema:

  >>> samples =          ('s1',         's2',        's3')
  >>> rows1 = [('l1',[ ('G', 'G'),   ('G', 'T'),  ('T', 'T')]),
  ...          ('l2',[ ('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> rows2 = [('l1',[(None, None),  ('T', 'T'),  ('G', 'G')]),
  ...          ('l3',[ ('A', 'A'),  (None, None), ('A', 'T')])]
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples).materialize()
  >>> genos2 = GenomatrixStream.from_tuples(rows2,'ldat',samples=samples).materialize()

  >>> genos = [genos1,genos2,genos1]
  >>> genos = merge_genomatrixstream_list(genos,VoteMerger())
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), (None, None), (None, None)])
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
  ('l3', [('A', 'A'), (None, None), ('A', 'T')])

  >>> genos = [ g.as_sdat() for g in [genos1,genos2,genos1] ]
  >>> genos = merge_genomatrixstream_list(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A'), ('A', 'A')])
  ('s2', [(None, None), ('T', 'T'), (None, None)])
  ('s3', [(None, None), ('A', 'T'), ('A', 'T')])

  Test really fast-path for homogeneous schema and disjoint rows:

  >>> samples =          ('s1',         's2',        's3')
  >>> rows1 = [('l1',[ ('G', 'G'),   ('G', 'T'),  ('T', 'T')]),
  ...          ('l2',[ ('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> rows2 = [('l3',[(None, None),  ('T', 'T'),  ('G', 'G')]),
  ...          ('l4',[ ('A', 'A'),  (None, None), ('A', 'T')])]
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples).materialize()
  >>> genos2 = GenomatrixStream.from_tuples(rows2,'ldat',samples=samples).materialize()

  >>> genos = [genos1,genos2]
  >>> genos = merge_genomatrixstream_list(genos,VoteMerger())
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), ('G', 'T'), ('T', 'T')])
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
  ('l3', [(None, None), ('T', 'T'), ('G', 'G')])
  ('l4', [('A', 'A'), (None, None), ('A', 'T')])
  '''
  assert mergefunc is not None

  if not genos:
    raise ValueError('Invalid empty genotype list')

  # Fastpath for nothing to merge
  if len(genos) == 1:
    return genos[0].merged(mergefunc)

  # Fastpath for identical schema: Concatenate all genotype rows and merge
  # using the single-matrix merge function.  Also optimizes the single input
  # case.
  format = genos[0].format
  if not all(g.format==format for g in genos):
    raise ValueError('Input genotypes must all be in same format')

  genome  = genos[0].genome
  genos   = [ g.transformed(recode_models=genome) for g in genos ]
  phenome = merge_phenome_list([g.phenome for g in genos ])

  columns = [ g.columns for g in genos ]
  if all(columns[0]==c for c in columns):
    # Detect if rows constitute a known disjoint partition
    rowset  = set()
    rowlist = []
    unique  = True
    for g in genos:
      rows = set(g.rows) if g.rows is not None else None
      if rows is None or len(rows)!=len(g.rows) or (rowset&rows):
        rowlist = None
        unique  = False
        break
      rowset |= rows
      rowlist.extend(g.rows)

    # Pass-through to merge_genomatrixsteam
    if format=='sdat':
      genos = genos[0].clone(chain(*genos),genome=genome,samples=rowlist,unique=unique,
                             materialized=False)
    else:
      models = []
      def _combine(genos):
        for g in genos:
          for labelrow,model in izip_exact(g,g.models):
            models.append(model)
            yield labelrow

      genos = genos[0].clone(_combine(genos),models=models,genome=genome,loci=rowlist,
                             unique=unique,materialized=False)

    # Merge if not unique
    genos = genos.merged(mergefunc)

    return genos

  # Slow path to handle heterogeneous columns
  new_rows      = {}
  new_columns   = {}
  merge_columns = defaultdict(lambda: ([],[]))
  merge_rows    = defaultdict(lambda: defaultdict(list))

  # We must materialize all streams, so ensure they are packed
  genos = [ g.transformed(repack=True) for g in genos ]

  for i,g in enumerate(genos):
    # Collect columns and form mappings from old schemas to new
    for j,column in enumerate(g.columns):
      k=new_columns.setdefault(column,len(new_columns))
      c=merge_columns[i]
      c[0].append(j)
      c[1].append(k)

    # Collect row labels and materialize all rows for later merging
    for row_label,row in g:
      new_rows.setdefault(row_label,len(new_rows))
      merge_rows[row_label][i].append(row)

  # Invert row and column dictionaries to recover insertion orderings
  new_columns = tuple(imap(itemgetter(0),sorted(new_columns.iteritems(),key=itemgetter(1))))
  new_rows    = tuple(imap(itemgetter(0),sorted(new_rows.iteritems(),   key=itemgetter(1))))
  n = len(new_columns)

  # FIXME: Refactor to clarify the following cases
  #        1) Identical columns
  #        2) Disjoint rows, overlapping columns
  #        3) Disjoint columns, overlapping rows
  #        4) Fully general merge

  if format=='ldat':
    models = []
    samples,loci = new_columns,new_rows
  else:
    models = [ genome.get_model(column) for column in new_columns ]
    samples,loci = new_rows,new_columns

  def _merger():
    # Fully general merge over duplicate rows and columns
    for label in new_rows:
      with gcdisabled():
        # Form null genotype lists at each new column
        # (place_list understands None is a null list)
        new_row = [None]*n

        # Iterate over input rows and schema, find the cooresponding column
        # mappings, and append the relevant genotypes
        for i,rows in merge_rows.pop(label).iteritems():
          j,k  = merge_columns[i]

          # All columns are extracted from rows of schema i
          if len(rows) == 1:
            # place_list is smart enough to handle bare genotypes as well as
            # lists, so use pick when possible
            cols = pick(rows[0], j)
          else:
            # Otherwise use the more general pick_columns
            cols = pick_columns(rows, j)

          # Place columns in the correct destination
          place_list(new_row, cols, k)

        # Merge genotypes
        if format=='ldat':
          model = genome.get_model(label)
          models.append(model)
          new_row = mergefunc.merge_locus(samples, label, model, new_row)

        else:
          new_row = mergefunc.merge_sample(label, loci, models, new_row)

      # Yield new row
      yield label,new_row

  return genos[0].clone(_merger(),samples=samples,loci=loci,models=models,genome=genome,
                                  packed=False,materialized=False,unique=True)


#######################################################################################


def build_genomatrixstream_from_genotriples(triples, format, mergefunc):
  '''
  Build genomatrix from genotriples using either the xtab or the rowsby
  function.  The rowsby function would be chosen over xtab if and only if:

    1. the necessary columns are given, and
    2. triples have been ordered appropriately (specified by the order argument).

  @param   triples: sequence of genotriples(str,str,genotype representation)
  @type    triples: sequence
  @param    format: format genomatrix
  @type     format: string
  @param mergefunc: function to merge multiple genotypes into a consensus genotype
  @type  mergefunc: callable
  @return         : genomatrix formed from the input triples
  @rtype          : genomatrix generator

  >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
  ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
  ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples).sorted('sample')
  >>> merge = UniqueMerger()
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merge)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A')])
  ('s2', [('G', 'T'), ('T', 'T')])
  ('s3', [('G', 'G'), ('A', 'A')])

  >>> triples = [('s1','l1', ('G', 'G')),('s1','l1', ('G', 'T')),
  ...            ('s2','l1', ('G', 'T')),('s2','l1', ('T', 'T')),
  ...            ('s3','l1', ('G', 'G')),('s3','l1', ('G', 'G'))]
  >>> triples = GenotripleStream.from_tuples(triples).sorted('sample')
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merge)
  >>> for row in genos:
  ...   print row
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found

  >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
  ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
  ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples).sorted('sample').materialize()
  >>> merger = VoteMerger()
  >>> triples.samples = None
  >>> triples.loci=['l1','l2']
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merger)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A')])
  ('s2', [('G', 'T'), ('T', 'T')])
  ('s3', [('G', 'G'), ('A', 'A')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merger = VoteMerger()
  >>> triples.loci=None
  >>> triples.samples=['s1','s2','s3']
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merger)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A')])
  ('s2', [('G', 'T'), ('T', 'T')])
  ('s3', [('G', 'G'), ('A', 'A')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merger = VoteMerger()
  >>> genos = build_genomatrixstream_from_genotriples(triples,'ldat',merger)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
  >>> sorted( (l,list(c)) for l,c in merger.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted( (l,list(c)) for l,c in merger.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]
  '''
  genome = triples.genome
  merger = mergefunc.merge_geno

  columns = None
  if format == 'sdat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(0),itemgetter(1),itemgetter(2)
    columns = tuple(sorted(triples.loci)) if triples.loci else None

    def aggfunc(sample,locus,genos):
      return merger(sample,locus,genome.get_model(locus),genos)

  elif format == 'ldat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(1),itemgetter(0),itemgetter(2)
    columns = tuple(sorted(triples.samples)) if triples.samples else None

    def aggfunc(locus,sample,genos):
      return merger(sample,locus,genome.get_model(locus),genos)

  else:
    raise NotImplementedError('triple to %s format conversion is not supported' % format)

  order = triples.order
  if format == 'sdat' and order != 'sample':
    order = False
  elif format == 'ldat' and order != 'locus':
    order = False

  # SLOWPATH: full xtab because of unknown columns or unordered rows
  if columns is None or not order:
    columns,rows,data = xtab(triples, rowkeyfunc, colkeyfunc, valuefunc, aggfunc)
    genos  = tuple(izip(rows,data))
    rows   = tuple(rows)
    if format=='ldat':
      models = [ genome.get_model(row) for row in rows ]
    else:
      models = [ genome.get_model(column) for column in columns ]

  # FASTPATH: rowsby using order and known columns
  else:
    columns,genos = rowsby(triples, columns, rowkeyfunc, colkeyfunc, valuefunc, aggfunc)
    rows = None

    models = []
    if format=='sdat':
      models = [ genome.get_model(locus) for locus in columns ]
    else:
      def _build(genos):
        for locus,row in genos:
          assert locus in genome.loci
          assert isinstance(row[0],Genotype)
          assert row[0].model is genome.loci[locus].model
          models.append( genome.get_model(locus) )
          yield locus,row
      genos = _build(genos)

  columns = tuple(columns)
  if format=='ldat':
    samples,loci = columns,rows
  else:
    samples,loci = rows,columns

  return GenomatrixStream(genos,format,samples=samples,loci=loci,models=models,
                                genome=genome,phenome=triples.phenome,unique=True)


#######################################################################################


def pack_genomatrixstream(genos):
  '''
  Transform a genomatrix into an internal packed representation

  @param      genos: genomatrix stream
  @type       genos: sequence
  @param     genome: genome descriptor
  @type      genome: Genome instance
  @return          : genomatrix with a packed internal format
  @rtype           : genomatrix generator

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[ ('A', 'A'), (None, None), ('G', 'G') ]),
  ...          ('l2',[(None, None),(None, None),(None, None)]),
  ...          ('l3',[ ('A', 'A'), (None, None),(None, None)]),
  ...          ('l4',[ ('G', 'T'), (None, None), ('T', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> genos = pack_genomatrixstream(genos)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])
  '''
  def _pack(genos):
    if genos.format=='sdat':
      descr = GenotypeArrayDescriptor(genos.models)
      for label,row in genos:
        yield label,GenotypeArray(descr,row)
    else:
      n = len(genos.columns)
      descrcache = {}

      for label,row in genos:
        model = row[0].model
        descr = descrcache.get(model)
        if descr is None:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        yield label,GenotypeArray(descr,row)

  if not genos.columns:
    return genos

  return genos.clone(_pack(genos),packed=True)


def rename_genomatrixstream_alleles(genos, rename_alleles, warn=False):
  '''
  Returns a new genomatrix with the alleles renamed

  @param           genos: genomatrix stream
  @type            genos: sequence
  @param  rename_alleles: rename alleles for any loci in the supplied dictionary from old allele name to new allele name
  @type   rename_alleles: dict from old_allele str -> new_allele str

  @return            : genomatrix with a packed internal format
  @rtype             : genomatrix generator

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[('A', 'A'),(None,None),('G', 'G')]),
  ...         ('l2',[(None,None),(None,None),(None,None)]),
  ...         ('l3',[('A', 'A'),(None,None),(None,None)]),
  ...         ('l4',[('G', 'T'),(None,None),('T', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples).materialize()

  >>> complement = {'A':'T','T':'A','C':'G','G':'C',None:None}
  >>> rename = {'l1':complement, 'l4':complement}
  >>> renamed = rename_genomatrixstream_alleles(genos,rename)
  >>> renamed.samples
  ('s1', 's2', 's3')
  >>> for row in renamed:
  ...   print row
  ('l1', [('T', 'T'), (None, None), ('C', 'C')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('A', 'C'), (None, None), ('A', 'A')])

  >>> loci = ('l1','l2','l3','l4')
  >>> rows = [('s1',[ ('A', 'A'),(None,None), ('A', 'A'),('G', 'T')]),
  ...         ('s2',[(None,None),(None,None),(None,None),(None,None)]),
  ...         ('s3',[ ('G', 'G'),(None,None),(None,None),('T', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)

  >>> renamed = rename_genomatrixstream_alleles(genos,rename)
  >>> renamed.loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in renamed:
  ...   print row
  ('s1', [('T', 'T'), (None, None), ('A', 'A'), ('A', 'C')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('C', 'C'), (None, None), (None, None), ('A', 'A')])
  '''
  # FIXME: optimized to only recode models that are remapped
  # FIXME: This need not convert to tuples as an intermediate step
  # FIXME: Strand remapping is not supported
  genome = Genome()

  def _rename():
    if genos.format=='ldat':
      for locus,row in genos:
        old_locus = genos.genome.get_locus(locus)

        genome.merge_locus(locus, chromosome=old_locus.chromosome,
                                  location=old_locus.location,
                                  strand=old_locus.strand, warn=warn)

        if locus in rename_alleles:
          r = rename_alleles[locus]
          row = [ ((r[g[0]],r[g[1]]) if g else g.alleles()) for g in row ]

        yield locus,row

    elif genos.format=='sdat':
      for locus in genos.loci:
        old_locus = genos.genome.get_locus(locus)
        genome.merge_locus(locus, chromosome=old_locus.chromosome,
                                  location=old_locus.location,
                                  strand=old_locus.strand, warn=warn)

      remaps = [ rename_alleles.get(h) for h in genos.loci ]
      for sample,row in genos:
        row = [ ((r[g[0]],r[g[1]]) if g and r else g.alleles()) for g,r in izip_exact(row,remaps) ]
        yield sample,row

    else:
      raise ValueError('Invalid genomatrixstream format')

  return GenomatrixStream.from_tuples(_rename(),genos.format,samples=genos.samples,
                                      loci=genos.loci,unique=genos.unique, genome=genome)


def rename_genotriples_alleles(triples, rename_alleles):
  '''
  Returns a new genotriple stream with the alleles renamed

  @param         triples: sequence of genotriples(str,str,genotype representation)
  @type          triples: sequence
  @param  rename_alleles: rename alleles for any loci in the supplied dictionary from old allele name to new allele name
  @type   rename_alleles: dict from old_allele str -> new_allele str
  @return               : genotriple with the new internal format
  @rtype                : genotriple generator

  >>> triples = [('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A')),
  ...            ('s2','l2', ('A', 'T')),('s1','l1', ('T', 'T')),
  ...            ('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> rename = {'l1':{'A':'T','T':'A','C':'G','G':'C',None:None}}
  >>> for row in rename_genotriples_alleles(triples,rename):
  ...   print row
  ('s3', 'l1', ('C', 'C'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l2', ('A', 'T'))
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l1', ('C', 'C'))
  ('s2', 'l2', ('A', 'A'))
  '''
  # FIXME: optimized to only recode models that are remapped
  # FIXME: Genome is not preserved!!!
  def _rename():
    for sample,locus,geno in triples:
      if locus in rename_alleles:
        remap = rename_alleles[locus]
        geno  = (remap[geno[0]],remap[geno[1]])
      else:
        geno  = geno.alleles()
      yield sample,locus,geno

  return GenotripleStream.from_tuples(_rename(),loci=triples.loci,samples=triples.samples,
                                      order=None,unique=triples.unique)


def filter_genomatrixstream_missing(genos):
  '''
  Filter samples or loci if all genotypes for a locus or sample are missing.
  Will result in full materialization of the dataset when a column contains
  only missing data.  If there is no such column, only the minimum necessary
  portion of the dataset is materialized.

  @param   genos: genomatrix stream
  @type    genos: sequence
  @return       : possibly materialized genotype matrix
  @rtype        : generator or list

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[ ('A','A'), (None,None), ('G','G') ]),
  ...         ('l2',[(None,None),(None,None),(None,None)]),
  ...         ('l3',[ ('A','A'), (None,None),(None,None)]),
  ...         ('l4',[ ('G','T'), (None,None), ('T','T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> genos = filter_genomatrixstream_missing(genos)
  >>> genos.samples
  ('s1', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('G', 'G')])
  ('l3', [('A', 'A'), (None, None)])
  ('l4', [('G', 'T'), ('T', 'T')])

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[(None,None),(None,None),(None,None)]),
  ...         ('l2',[ ('A', 'A'), ('T','T'),  ('A','T')]),
  ...         ('l3',[(None,None),(None,None),(None,None)]),
  ...         ('l4',[ ('G', 'T'),(None,None), ('T','T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> genos = filter_genomatrixstream_missing(genos)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'T')])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])
  '''

  # Helper generator function used to filter for completely empty rows.
  # Since the columns_notseen set is invariant to these rows, they may be
  # silently filtered.  This is also desirable because if all columns are
  # seen, then only this simple filter is required to process the remaining
  # rows and the rest need not be materialized.  This can be a huge
  # performance improvement for many datasets with no missing columns.

  if genos.format=='ldat':
    models = []
    def _filter():
      for (lname,row),model in izip_exact(genos,genos.models):
        if any(row):
          models.append(model)
          yield lname,row

    new_genos = genos.clone(_filter(),loci=None,models=models,materialized=False)
  else:
    def _filter():
      for lname,row in genos:
        if any(row):
          yield lname,row

    new_genos = genos.clone(_filter(),samples=None,materialized=False)

  rows = []

  # Set of column indices not yet observed to have data
  columns_notseen = set(range(len(new_genos.columns)))

  data = iter(new_genos)
  for lname,row in data:
    # Store row in materialized list
    rows.append( (lname,row) )

    # Remove any column indices with non-missing data
    columns_notseen.difference_update( [ i for i in columns_notseen if row[i] ] )

    # Stop materializing if there are no more indices to remove
    if not columns_notseen:
      return new_genos.clone(chain(rows,data),materialized=False)

  # Full materialize was necessary and some columns need to be filtered
  columns_notseen = set(pick(new_genos.columns, columns_notseen))
  new_genos = new_genos.clone(rows)
  return filter_genomatrixstream_by_column(new_genos,columns_notseen,exclude=True)


def filter_genotriples_missing(triples):
  '''
  Filter out the genotriple if its genotype is missing.

  @param triples: sequence of genotriples(str,str,genotype representation)
  @type  triples: sequence
  @return       : iterator with missing genotype in the triples being filtered out.
  @rtype        : iterator

  >>> triples=[('l1','s1',(None,None)),
  ...          ('l2','s2', ('A','T' ))]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> list(filter_genotriples_missing(triples))
  [('l2', 's2', ('A', 'T'))]
  '''
  return triples.clone(ifilter(itemgetter(2), triples))


def build_genotriples_from_genomatrix(genos):
  '''
  Generate genotype triples from the locus major genotype matrix.
  @param   genos: genomatrix stream
  @type    genos: sequence
  @rtype:         generator
  @returns:       a genotype triplet stream

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['  ','CT','  ']),
  ...         ('l3',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)
  >>> triples = build_genotriples_from_genomatrix(genos)
  >>> sorted(triples.samples)
  ['s1', 's2', 's3']
  >>> for s,l,g in triples:
  ...   print s,l,g
  s1 l1 ('A', 'A')
  s2 l1 ('A', 'G')
  s3 l1 ('G', 'G')
  s1 l2 (None, None)
  s2 l2 ('C', 'T')
  s3 l2 (None, None)
  s1 l3 ('A', 'A')
  s2 l3 ('A', 'T')
  s3 l3 ('T', 'T')

  >>> loci =        ('l1','l2','l3')
  >>> rows = [('s1',['AA','AG','GC']),
  ...         ('s2',['AT','GG','CC']),
  ...         ('s3',['  ','AA','  '])]
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> triples = build_genotriples_from_genomatrix(genos)
  >>> sorted(triples.loci)
  ['l1', 'l2', 'l3']
  >>> for s,l,g in triples:
  ...   print s,l,g
  s1 l1 ('A', 'A')
  s1 l2 ('A', 'G')
  s1 l3 ('C', 'G')
  s2 l1 ('A', 'T')
  s2 l2 ('G', 'G')
  s2 l3 ('C', 'C')
  s3 l1 (None, None)
  s3 l2 ('A', 'A')
  s3 l3 (None, None)
  '''
  if genos.format=='ldat':
    order = 'locus'

    def _build():
      samples = genos.samples
      for i,(locus,row) in enumerate(genos):
        for sample,geno in izip_exact(samples, row):
          yield sample,locus,geno
  else:
    order = 'sample'
    def _build():
      loci = genos.loci
      for sample,row in genos:
        for locus,geno in izip_exact(loci, row):
          yield sample,locus,geno

  # Unset the result order if any of the required ordering constraints
  # cannot be verified
  if (genos.samples is None or genos.loci is None
                            or sorted(genos.samples) != list(genos.samples)
                            or sorted(genos.loci)    != list(genos.loci)):
    order = None

  return GenotripleStream(_build(),samples=genos.samples,loci=genos.loci,
                                   genome=genos.genome,phenome=genos.phenome,
                                   order=order, unique=genos.unique)


def _genome_merge_loci(old_genome, old_name, new_genome, new_name, warn):
  '''
  Helper to merge loci
  '''
  old_locus = old_genome.loci[old_name]
  old_model = old_locus.model

  new_genome.merge_locus(new_name, fixed=old_locus.fixed,
                                   chromosome=old_locus.chromosome,
                                   location=old_locus.location,
                                   strand=old_locus.strand, warn=warn)

  new_locus = new_genome.loci[new_name]
  new_model = new_locus.model

  if new_model is None:
    new_locus.model = old_model
  elif old_model is None or new_model is not old_model:
    new_locus.model = None
    return True

  return False


def _genome_rename(old_genome, locusmap, warn=False):
  '''
  Helper to merge loci
  '''
  if not locusmap:
    return old_genome,False

  recode = False

  new_genome = Genome()

  for old_name in old_genome.loci:
    new_name = locusmap.get(old_name,old_name)
    recode |= _genome_merge_loci(old_genome, old_name, new_genome, new_name, warn)

  return new_genome,recode


def _phenome_rename(old_phenome, samplemap, warn=False):
  '''
  Helper to merge loci
  '''
  if not samplemap or old_phenome is None:
    return old_phenome

  new_phenome = Phenome()

  # FIXME: Parents are not renamed correctly
  for old_name in old_phenome.phenos:
    new_name = samplemap.get(old_name,old_name)
    _phenome_merge_individuals(old_phenome, old_name, new_phenome, new_name, warn)

  return new_phenome


def _phenome_merge_individuals(old_phenome, old_name, new_phenome, new_name, warn):
  old_phenos = old_phenome.get_phenos(old_name)

  # FIXME: individual name is muddled
  # FIXME: Parents are not renamed correctly
  new_phenome.merge_phenos(new_name,family     = old_phenos.family,
                                    parent1    = old_phenos.parent1,
                                    parent2    = old_phenos.parent2,
                                    phenoclass = old_phenos.phenoclass,
                                    sex        = old_phenos.sex,
                                    warn       = warn)


def rename_genotriples(triples,samplemap,locusmap,warn=False):
  '''
  Rename the sample and locus for each genotriple according to the samplemap
  and locusmap. If there is not mapping for a particular sample or locus,
  the original name will be used.

  @param   triples: sequence of genotriples(str,str,genotype representation)
  @type    triples: sequence
  @param samplemap: map between the current sample name and a new name
  @type  samplemap: dict
  @param  locusmap: map between the current locus name and a new name
  @type   locusmap: dict
  @return         : renamed sequence of genotriples
  @rtype          : generator

  >>> triples = [('s1','l1','AT'),('s1','l2','AG'),
  ...            ('s2','l1','TT'),('s2','l2','AA')]
  >>> triples = GenotripleStream.from_strings(triples,snp)

  >>> samplemap = dict([('s1','S1'),('s2','S2')])
  >>> locmap    = dict([('l1','L1'),('l2','L2')])
  >>> triples   = rename_genotriples(triples,samplemap,locmap)
  >>> for sample,loc,geno in triples:
  ...   print sample,loc,geno
  S1 L1 ('A', 'T')
  S1 L2 ('A', 'G')
  S2 L1 ('T', 'T')
  S2 L2 ('A', 'A')
  '''
  samplemap = samplemap or {}
  locusmap  = locusmap  or {}

  genome,recode = _genome_rename(triples.genome,   locusmap,  warn)
  phenome       = _phenome_rename(triples.phenome, samplemap, warn)

  def _rename(triples):
    for sample,locus,geno in triples:
      sample = samplemap.get(sample,sample)
      locus  = locusmap.get(locus,locus)
      yield sample,locus,geno

  samples = triples.samples
  if samples and samplemap:
    samples = set(samplemap.get(s,s) for s in samples)

  loci = triples.loci
  if loci and locusmap:
    loci = set(locusmap.get(l,l) for l in loci)

  # FIXME: We can do more to prove uniqueness
  triples = triples.clone(_rename(triples),samples=samples,loci=loci,genome=genome,order=None,materialized=False)

  if recode:
    triples = triples.transformed(recode_models=genome)

  return triples


def filter_genotriples(triples,sampleset,locusset,exclude=False):
  '''
  Filter out the genotriples against the sampelset and locusset.
  Depending on the value of exclude flag, both sets will be treated
  as either inclusion sets(exclude=False) or exclusion sets(exclude=True).

  @param   triples: sequence of genotriples(str,str,genotype representation)
  @type    triples: sequence
  @param sampleset: set of sample names
  @type  sampleset: set
  @param  locusset: set of locus names
  @type   locusset: set
  @param   exclude: flag to exclude rather than include items in sampleset
                    and locusset
  @type    exclude: bool
  @return         : filtered genotriples
  @rtype          : generator

  >>> triples = [('s1','l1','AT'),('s1','l2','AA'),
  ...            ('s2','l1','TT'),('s2','l2','AG')]
  >>> triples = GenotripleStream.from_strings(triples,snp).materialize()

  >>> sset = set(['s1'])
  >>> lset = set(['l1','l2'])
  >>> for s,l,g in filter_genotriples(triples,sset,lset):
  ...   print s,l,g
  s1 l1 ('A', 'T')
  s1 l2 ('A', 'A')

  >>> lset = None
  >>> for s,l,g in filter_genotriples(triples,sset,lset,exclude=True):
  ...   print s,l,g
  s2 l1 ('T', 'T')
  s2 l2 ('A', 'G')
  '''
  if sampleset is not None:
    sampleset = as_set(sampleset)

  if locusset is not None:
    locusset = as_set(locusset)

  def _filter():
    if exclude:
      for sample,locus,geno in triples:
        if sampleset is not None and sample in sampleset:
          continue
        if locusset is not None and locus in locusset:
          continue
        yield sample,locus,geno
    else:
      for sample,locus,geno in triples:
        if sampleset is not None and sample not in sampleset:
          continue
        if locusset is not None and locus not in locusset:
          continue
        yield sample,locus,geno

  loci = triples.loci
  if loci is not None and locusset is not None:
    if exclude:
      loci = set(l for l in loci if l not in locusset)
    else:
      loci = set(l for l in loci if l in locusset)
  elif loci is None and locusset is not None and not exclude:
    loci = locusset

  samples = triples.samples
  if samples is not None and sampleset is not None:
    if exclude:
      samples = set(s for s in samples if s not in sampleset)
    else:
      samples = set(s for s in samples if s in sampleset)
  elif samples is None and sampleset is not None and not exclude:
    samples = sampleset

  return triples.clone(_filter(),loci=loci,samples=samples,order=None)


def rename_genomatrixstream_column(genos,colmap,warn=False):
  '''
  Rename the columns for the genotype matrix data
  according to a name mapping. If the name of the column
  is not in the mapping, its original name will be used.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data
  @type     rows: sequence
  @param  colmap: map between the current column label and a new label
  @type   colmap: dict
  @rtype        : sequence
  @return       : genotype matrix with renamed columns

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> colmap = {'s1':'S1','s2':'S2','s3':'S3'}
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)
  >>> genos = rename_genomatrixstream_column(genos,colmap)
  >>> genos.samples
  ('S1', 'S2', 'S3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('A', 'A'), ('A', 'T'), ('T', 'T')])

  >>> loci = ('l1','l2','l3')
  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'G') ]),
  ...         ('s2',[(None, None),  ('T','T'),   (None, None)]),
  ...         ('s3',[ ('A', 'A'),  (None, None), (None, None)]),
  ...         ('s4',[ ('A', 'T'),  (None, None),  ('T', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci).materialize()

  >>> colmap = {'l2':'l1'}
  >>> genos = rename_genomatrixstream_column(genos,colmap)
  >>> genos.loci
  ('l1', 'l1', 'l3')
  >>> genos.unique
  False
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('G', 'G')])
  ('s2', [(None, None), ('T', 'T'), (None, None)])
  ('s3', [('A', 'A'), (None, None), (None, None)])
  ('s4', [('A', 'T'), (None, None), ('T', 'T')])
  '''
  new_columns = tuple(colmap.get(name,name) for name in genos.columns)

  if new_columns == genos.columns:
    return genos

  unique_columns = (len(genos.columns) == len(set(new_columns)))
  unique = genos.unique and unique_columns

  if genos.format=='ldat':
    phenome = Phenome()

    for old_name in genos.columns:
      new_name = colmap.get(old_name, old_name)
      _phenome_merge_individuals(genos.phenome, old_name, phenome, new_name, warn)

    genos = genos.clone(genos.use_stream(), samples=new_columns, phenome=phenome, unique=unique)
  else:
    genome = Genome()

    recode = False
    for old_name in genos.columns:
      new_name = colmap.get(old_name,old_name)
      recode |= _genome_merge_loci(genos.genome, old_name, genome, new_name, warn)

    genos = genos.clone(genos.use_stream(), loci=new_columns, genome=genome, unique=unique)

    if recode:
      genos = genos.transformed(recode_models=genome)

  return genos


def filter_genomatrixstream_by_column(genos,colset,exclude=False):
  '''
  Filter the genotype matrix data by a column set.  Depending on the
  value of the exclude flag, the column set will be used either for
  inclusion (exclude=False) or exclusion (exclude=True).

  @param   genos: genomatrix stream
  @type    genos: sequence
  @param  colset: set of the column names
  @type   colset: set
  @param exclude: flag to exclude rather than include items in colset
  @type  exclude: bool
  @return       : filtered genotype matrix
  @rtype        : sequence

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> colset = set(['s1','s3'])
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)
  >>> genos = filter_genomatrixstream_by_column(genos,colset)
  >>> genos.samples
  ('s1', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('G', 'G')])
  ('l2', [('A', 'A'), ('T', 'T')])

  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)
  >>> genos = filter_genomatrixstream_by_column(genos,colset,exclude=True)
  >>> genos.samples
  ('s2',)
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'G')])
  ('l2', [('A', 'T')])

  >>> loci =  ('l1','l2')
  >>> rows = [('s1',['AA','AA']),
  ...         ('s2',['AG','AT']),
  ...         ('s3',['GG','TT'])]
  >>> colset = set(['l1','l3'])
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> genos = filter_genomatrixstream_by_column(genos,colset)
  >>> genos.loci
  ('l1',)
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A')])
  ('s2', [('A', 'G')])
  ('s3', [('G', 'G')])

  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci)
  >>> genos = filter_genomatrixstream_by_column(genos,colset,exclude=True)
  >>> genos.loci
  ('l2',)
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A')])
  ('s2', [('A', 'T')])
  ('s3', [('T', 'T')])
  '''
  colset = as_set(colset)

  columns = genos.columns
  if exclude:
    columns = ((name,i) for i,name in enumerate(columns) if name not in colset)
  else:
    columns = ((name,i) for i,name in enumerate(columns) if name     in colset)

  columns,indices = tuple(izip(*columns)) or ((),())

  if columns == genos.columns:
    return genos

  def _filter():
    for label,row in genos:
      yield label,pick(row,indices)

  if genos.format=='ldat':
    new_genos=genos.clone(_filter(),samples=columns,packed=False,materialized=False)
  else:
    assert len(genos.columns) == len(genos.models)
    models = pick(genos.models,indices)
    new_genos=genos.clone(_filter(),loci=columns,models=models,packed=False,materialized=False)

  return new_genos


def reorder_genomatrixstream_columns(genos,labels):
  '''
  Reorder and filter the genotype matrix columns to match the sequence of
  labels provided.  If not all labels appear, the remainder are retained at
  the end of the list in lexicographical order.

  @param   genos: genomatrix stream
  @type    genos: sequence
  @param  labels: desired column labels
  @type   labels: sequence
  @return       : genomatrix generator
  @rtype        : generator

  >>> samples =        ('s1','s2','s3')
  >>> rows    = [('l1',['AA','AG','GG']),
  ...            ('l2',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples).materialize()

  >>> reordered = reorder_genomatrixstream_columns(genos,['s2','s1','s3','s4'])
  >>> reordered.samples
  ('s2', 's1', 's3')
  >>> for row in reordered:
  ...   print row
  ('l1', [('A', 'G'), ('A', 'A'), ('G', 'G')])
  ('l2', [('A', 'T'), ('A', 'A'), ('T', 'T')])

  >>> reordered = reorder_genomatrixstream_columns(genos,['s2','s1'])
  >>> reordered.samples
  ('s2', 's1', 's3')
  >>> for row in reordered:
  ...   print row
  ('l1', [('A', 'G'), ('A', 'A'), ('G', 'G')])
  ('l2', [('A', 'T'), ('A', 'A'), ('T', 'T')])

  >>> reorder_genomatrixstream_columns(genos,['s2','s2'])
  Traceback (most recent call last):
       ...
  ValueError: Duplicated column label: s2
  '''
  columnset = set(genos.columns)
  labelset  = set(labels)
  extras    = sorted(l for l in genos.columns if l not in labelset)
  order     = [ l for l in labels if l in columnset ] + list(extras)
  remap     = dict( (l,i) for i,l in enumerate(genos.columns) )

  try:
    indices = [ remap.pop(l) for l in order ]
  except KeyError:
    raise ValueError('Duplicated column label: %s' % l)

  new_columns = pick(genos.columns,indices)
  rows = iter(genos)
  def _reorder():
    for label,row in rows:
      yield label,pick(row, indices)

  models = genos.models
  if genos.format=='sdat' and models is not None:
    models = pick(models,indices)

  if genos.format=='ldat':
    new_genos = genos.clone(_reorder(),samples=new_columns,materialized=False,packed=False)
  else:
    models = genos.models
    if models is not None:
      models = pick(models,indices)
    new_genos = genos.clone(_reorder(),loci=new_columns,models=models,materialized=False,packed=False)

  return new_genos


def reorder_genomatrixstream_rows(genos, labels):
  '''
  Reorder and filter the genotype matrix rows to match the sequence of
  labels provided.  If not all labels appear, the remainder are retained at
  the end of the list in in lexicographical order.

  @param   genos: genomatrix stream
  @type    genos: sequence
  @param  labels: desired row labels
  @type   labels: sequence
  @return       : genomatrix generator
  @rtype        : generator

  >>> loci =        ('l1','l2','l3')
  >>> rows = [('s1',['AA','AG','CT']),
  ...         ('s2',['TT','GG','CC']),
  ...         ('s3',['AA','GG','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'sdat',snp,loci=loci).materialize()

  >>> reordered = reorder_genomatrixstream_rows(genos,['s2','s1','s3','s4'])
  >>> reordered.loci
  ('l1', 'l2', 'l3')
  >>> for row in reordered:
  ...   print row
  ('s2', [('T', 'T'), ('G', 'G'), ('C', 'C')])
  ('s1', [('A', 'A'), ('A', 'G'), ('C', 'T')])
  ('s3', [('A', 'A'), ('G', 'G'), ('T', 'T')])

  >>> reordered = reorder_genomatrixstream_rows(genos,['s2','s1'])
  >>> reordered.loci
  ('l1', 'l2', 'l3')
  >>> for row in reordered:
  ...   print row
  ('s2', [('T', 'T'), ('G', 'G'), ('C', 'C')])
  ('s1', [('A', 'A'), ('A', 'G'), ('C', 'T')])
  ('s3', [('A', 'A'), ('G', 'G'), ('T', 'T')])

  >>> reordered = reorder_genomatrixstream_rows(genos,['s2','s2'])
  Traceback (most recent call last):
       ...
  ValueError: Duplicated row label: s2
  '''
  genos    = genos.transformed(repack=True).materialize()
  rowset   = set(genos.rows)
  labelset = set(labels)
  extras   = sorted(r for r in genos.rows if r not in labelset)
  order    = [ l for l in labels if l in rowset ] + extras
  remap    = dict( (r,i) for i,r in enumerate(genos.rows))

  try:
    indices = [ remap.pop(l) for l in order ]
  except KeyError:
    raise ValueError('Duplicated row label: %s' % l)

  def _reorder():
    data = genos.use_stream()
    for i in indices:
      yield data[i]

  if genos.format=='ldat':
    models = genos.models
    if models is not None:
      models = pick(models,indices)
    new_genos = genos.clone(_reorder(),loci=order,models=models,unique=True,materialized=False)
  else:
    new_genos = genos.clone(_reorder(),samples=order,unique=True,materialized=False)

  return new_genos


def rename_genomatrixstream_row(genos,rowmap,warn=False):
  '''
  Rename the row label for the genotype matrix data.
  If there is no mapping for a row label, then the original
  label will be used.

  @param   genos: genomatrix stream
  @type    genos: sequence
  @param  rowmap: map between the current row label and a new label
  @type   rowmap: dict
  @return       : genotype matrix with renamed row labels
  @rtype        : generator

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)

  >>> rowmap = {'l1':'L1','l2':'L2'}
  >>> genos = rename_genomatrixstream_row(genos,rowmap)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('L1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('L2', [('A', 'A'), ('A', 'T'), ('T', 'T')])

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples).as_sdat()

  >>> rowmap = {'s1':'S1','s2':'S2','s3':'S3'}
  >>> genos = rename_genomatrixstream_row(genos,rowmap)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('S1', [('A', 'A'), ('A', 'A')])
  ('S2', [('A', 'G'), ('A', 'T')])
  ('S3', [('G', 'G'), ('T', 'T')])
  '''
  if genos.rows is not None:
    rows = tuple(rowmap.get(label,label) for label in genos.rows)
  else:
    rows = None

  if genos.format=='ldat':
    genome = Genome()

    recode = [False]
    def _rename():
      for locus,row in genos:
        new_locus = rowmap.get(locus,locus)
        recode[0] |= _genome_merge_loci(genos.genome, locus, genome, new_locus, warn)

        yield new_locus,row

    new_genos = genos.clone(_rename(),loci=rows,genome=genome,materialized=False)

    if recode[0]:
      new_genos = new_genos.transformed(recode_models=genome)

  else:
    phenome = Phenome()

    def _rename():
      for sample,row in genos:
        new_name = rowmap.get(sample,sample)
        _phenome_merge_individuals(genos.phenome, sample, phenome, new_name, warn)

        yield new_name,row

    new_genos = genos.clone(_rename(),samples=rows,phenome=phenome,materialized=False)

  return new_genos


def filter_genomatrixstream_by_row(genos,rowset,exclude=False):
  '''
  Filter the genotype matrix by a row set.
  Depending on the value of the exclude flag, the row set will
  be used either for inclusion(exclude=False) or exclusion(exclude=True)
  purpose.

  @param   genos: genomatrix stream
  @type    genos: sequence
  @param  rowset: set of the row names
  @type   rowset: set
  @param exclude: flag to exclude rather than include items in rowset
  @type  exclude: bool
  @return       : filtered genotype matrix
  @rtype        : sequence

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples).materialize()
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('A', 'A'), ('A', 'T'), ('T', 'T')])

  >>> rowset = set(['l1'])
  >>> new_genos = filter_genomatrixstream_by_row(genos,rowset)
  >>> new_genos.samples
  ('s1', 's2', 's3')
  >>> for row in new_genos:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])

  >>> new_genos = filter_genomatrixstream_by_row(genos,rowset,exclude=True)
  >>> new_genos.samples
  ('s1', 's2', 's3')
  >>> for row in new_genos:
  ...   print row
  ('l2', [('A', 'A'), ('A', 'T'), ('T', 'T')])

  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples).as_sdat().materialize()
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), ('A', 'A')])
  ('s2', [('A', 'G'), ('A', 'T')])
  ('s3', [('G', 'G'), ('T', 'T')])

  >>> rowset = set(['s1'])
  >>> new_genos = filter_genomatrixstream_by_row(genos,rowset)
  >>> new_genos.loci
  ('l1', 'l2')
  >>> for row in new_genos:
  ...   print row
  ('s1', [('A', 'A'), ('A', 'A')])
  >>> new_genos = filter_genomatrixstream_by_row(genos,rowset,exclude=True)
  >>> new_genos.loci
  ('l1', 'l2')
  >>> for row in new_genos:
  ...   print row
  ('s2', [('A', 'G'), ('A', 'T')])
  ('s3', [('G', 'G'), ('T', 'T')])
  '''
  rowset = as_set(rowset)

  rows = genos.rows
  if rows is not None:
    if exclude:
      rows = tuple(label for label in rows if label not in rowset)
    else:
      rows = tuple(label for label in rows if label     in rowset)

    if rows == genos.rows:
      return genos

  if genos.format=='sdat':
    if exclude:
      def _filter():
        for subject,row in genos:
          if subject not in rowset:
            yield subject,row
    else:
      def _filter():
        for subject,row in genos:
          if subject in rowset:
            yield subject,row

    return genos.clone(_filter(),samples=rows,materialized=False)

  else:
    models = []
    if exclude:
      def _filter():
        for (locus,row),model in izip_exact(genos,genos.models):
          if locus not in rowset:
            models.append(model)
            yield locus,row
    else:
      def _filter():
        for (locus,row),model in izip_exact(genos,genos.models):
          if locus in rowset:
            models.append(model)
            yield locus,row

    return genos.clone(_filter(),loci=rows,models=models,materialized=False)


def transpose_generator(columns, rows, m=32):
  '''
  Transpose a matrix of row labels and row data generating one column of
  data at a time.  Return a generator of the columns of data stored in rows,
  yielding sucessive column labels and column data.

  Requires the input data to be fully materialized, but allows results to be
  streamed.

  @param  columns: matrix column names
  @type   columns: sequence of strs
  @param     rows: genotype matrix data
  @type      rows: sequence
  @return        : tuple of column labels and generator of tuples of row labels and row data
  @rtype         : tuple

  >>> r = [('r1','abc'),
  ...      ('r2','def'),
  ...      ('r3','ghi')]
  >>> rowlabels,c = transpose_generator(['c1','c2','c3'],r)
  >>> rowlabels
  ('r1', 'r2', 'r3')
  >>> for clabel,cdata in c:
  ...   print clabel,cdata
  c1 ['a', 'd', 'g']
  c2 ['b', 'e', 'h']
  c3 ['c', 'f', 'i']
  '''
  n = len(rows)

  if not n or not columns:
    return [],[]

  rowlabels,rows = zip(*rows)

  if m <= 1:
    def _transpose_generator():
      for j,column in enumerate(columns):
        yield column,pick_columns(rows,j)
  else:
    def _transpose_generator():
      n = len(columns)
      for i in xrange(n//m+bool(n%m)):
        cslice = slice(m*i, m*(i+1))
        cols = pick_columns(rows, cslice)
        for label,col in izip(columns[cslice],cols):
          yield label,col

  return rowlabels,_transpose_generator()


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
