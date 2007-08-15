# -*- coding: utf-8 -*-
'''
File:          streams.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU genotype data management objects

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from   types             import NoneType
from   operator          import itemgetter, getitem
from   collections       import defaultdict
from   itertools         import izip,ifilter,imap,chain,groupby,repeat

from   glu.lib.utils     import pick,tally,unique
from   glu.lib.fileutils import load_list,load_map
from   glu.lib.imerge    import imerge
from   glu.lib.xtab      import xtab,rowsby

from   reprs             import snp
from   transform         import GenoTransform, prove_unique_transform
from   merge             import UniqueMerger, VoteMerger, mergefunc_transpose_adapter
from   genoarray         import UnphasedMarkerModel,GenotypeArrayDescriptor,GenotypeArray,Genotype, \
                                model_from_alleles

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
    @rtype :  sequence of sample, locus, and genotype
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
      raise RuntimeError, 'Genotriple stream already used'

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

  def __init__(self, triples, samples=None, loci=None, models=None, order=None,
                              unique=False, materialized=False):
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

    @param      triples: genotriple sequence
    @type       triples: sequence of tuples of sample, locus, and genotype
    @param      samples: optional set of samples refered to by the triples
    @type       samples: sequence, set, or None
    @param         loci: optional set of samples refered to by the triples
    @type          loci: sequence, set, or None
    @param        order: ordering of the triple stream, 'sample' or 'locus', or None
    @type         order: str or None
    @param       unique: flag indicating if repeated elements do not exist within the stream
    @type        unique: bool
    @param materialized: flag indicating if genos is a materialized
                         (multiply iterable) data type
    @type  materialized: bool
    '''
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
    self.models       = models
    self.unique       = bool(unique)
    self.materialized = materialized or isinstance(triples, (list,tuple))

  @staticmethod
  def from_streams(genos, mergefunc=None, order=None):
    '''
    Combine multiple genostreams into one genotriple stream

    @param     genos: genostreams
    @type      genos: list
    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @param     order: ordering of the triple stream, 'sample' or 'locus', or None
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
    ('s2', 'l2', ('T', 'T'))
    ('s3', 'l1', ('T', 'T'))
    ('s3', 'l2', (None, None))
    '''
    if not genos:
      raise ValueError('empty stream list')

    # FIXME: Consider merging before converting to triples if inputs are all
    #        sdat/ldat
    triples = [ g.as_genotriples() for g in genos ]

    if len(triples) == 1:
      return triples[0]

    if mergefunc is not None and order is None:
      order = 'sample'

    if order is not None:
      triples = [ t.sorted(order) for t in triples ]

    if order is None and mergefunc is None:
      triples = combine_unsorted_genotriple_list(triples)

    else:
      triples = combine_sorted_genotriple_list(triples)

    if mergefunc is not None:
      triples = triples.merged(mergefunc)

    assert isinstance(triples,GenotripleStream)

    return triples

  @staticmethod
  def from_tuples(triples, samples=None, loci=None, order=None, unique=False, modelmap=None):
    '''
    Alternate constructor that builds a new GenotripleStream object from a
    sequence of triples with genotypes in tuple
    '''
    triples = encode_genotriples_from_tuples(triples, modelmap)
    return GenotripleStream(triples, samples=samples, loci=loci, order=order, models=modelmap, unique=unique)

  @staticmethod
  def from_strings(triples, genorepr, samples=None, loci=None, order=None, unique=False, modelmap=None):
    '''
    Alternate constructor that builds a new GenotripleStream object from a
    sequence of triples with genotypes in a string format
    '''
    triples = encode_genotriples_from_strings(triples, genorepr, modelmap)
    return GenotripleStream(triples, samples=samples, loci=loci, order=order, models=modelmap, unique=unique)

  def clone(self, triples, **kwargs):
    '''
    Alternative constructor that builds a new GenotripleStream object
    with attributes based on self, but updated with the specified keyword
    arguments.

    For example:

    stream.clone(triples, materialized=False)

     is equivalent to

    GenotripleStream(triples, samples=stream.samples, loci=stream.loci,
                              order=stream.order, unique=stream.unique,
                              materialized=False)
    '''
    kwargs.setdefault('samples',      self.samples)
    kwargs.setdefault('loci',         self.loci)
    kwargs.setdefault('order',        self.order)
    kwargs.setdefault('models',       self.models)
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

    return self.clone(list(self.use_stream()), materialized=True)

  def transformed(self, transform=None, mergefunc=None, **kwargs):
    '''
    Apply filtering and renaming transformations to a genotriple stream.
    Transformations may be specified by key-word arguments or a
    transformation object, though not both.  Inclusions and exclusions are
    always performed before renaming operations.  If a merge function is
    specified when performing renaming transformations, the resulting
    genotriple stream will be guaranteed to contain unique rows and columns
    (see merged function).

    @param           genos: genotriple stream
    @type            genos: sequence
    @param       transform: transformation object (optional)
    @type        transform: GenoTransform object
    @param       mergefunc: function to merge multiple genotypes into a consensus genotype
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
    >>> for row in GenotripleStream(triples).transformed(include_loci=['l1'],exclude_samples=['s3']):
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s2', 'l1', ('G', 'T'))
    '''
    if transform is not None and kwargs:
      raise ValueError, 'specification of both a transform object and keyword arguments is not supported'
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

    if transform.rename_alleles is not None:
      triples = rename_genotriples_alleles(triples, transform.rename_alleles)

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
    >>> for row in GenotripleStream(triples).sorted():
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

    @param  mergefunc: function to merge multiple genotypes into a consensus genotype
    @type   mergefunc: callable
    @param      order: sort order, either 'samples', 'locus'.  Default is 'sample'
    @type       order: str
    @return          : sorted and merged genotriple stream
    @rtype           : GenotripleStream

    >>> triples = [('s1','l1', ('T', 'T')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('A', 'A')),
    ...            ('s1','l1', ('G', 'G')),('s2','l2', ('A', 'A'))]
    >>> triples = GenotripleStream.from_tuples(triples)
    >>> for row in GenotripleStream(triples).merged(VoteMerger(),order='locus'):
    ...   print row
    ('s1', 'l1', (None, None))
    ('s2', 'l1', ('G', 'T'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l2', ('A', 'A'))
    '''
    if self.unique:
      return self

    if self.order not in ('sample','locus'):
      return self.sorted(order).merged(mergefunc)

    return merge_sorted_genotriples(self,mergefunc)

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

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> triples = [('s1','l1',('G', 'G')),('s1','l2',('A', 'A')),
    ...            ('s2','l1',('G', 'T')),('s2','l2',('T', 'T')),
    ...            ('s3','l1',('G', 'G')),('s3','l2',('A', 'A'))]
    >>> triples = list(GenotripleStream.from_tuples(triples))

    >>> merge = VoteMerger()
    >>> ldat = GenotripleStream(triples).as_ldat(merge)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = GenotripleStream(triples,loci=['l1','l2']).as_ldat(merge)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = GenotripleStream(triples,samples=['s1','s2','s3']).as_ldat(merge)
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

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> triples = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
    ...            ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
    ...            ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
    >>> triples = list(GenotripleStream.from_tuples(triples))

    >>> merge = VoteMerger()
    >>> sdat = GenotripleStream(triples).as_sdat(merge)
    >>> sdat.loci
    ('l1', 'l2')
    >>> sdat.samples
    ('s1', 's2', 's3')
    >>> for row in sdat:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])

    >>> sdat = GenotripleStream(triples,loci=['l1','l2']).as_sdat(merge)
    >>> sdat.loci
    ('l1', 'l2')
    >>> sdat.samples
    ('s1', 's2', 's3')
    >>> for row in sdat:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('G', 'G'), ('A', 'A')])

    >>> sdat = GenotripleStream(triples,samples=['s1','s2','s3']).as_sdat(merge)
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
                     unique=True, materialized=False, packed=False):
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

    @param        genos: genomatrix sequence
    @type         genos: genomatrix
    @param       format: format of input genomatrix, either 'ldat' or 'sdat'
    @type        format: str
    @param      samples: list of samples, required for ldat
    @type       samples: list
    @param         loci: list of loci, required for sdat
    @type          loci: list
    @param       unique: flag indicating if repeated row or column elements do not exist
    @type        unique: bool
    @param materialized: flag indicating if this stream is materialized and
                         allows iteration multiple times
    @type  materialized: bool
    @param       packed: flag indicating if genotypes are packed into a
                         compressed array format
    @type        packed: bool
    '''
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
    self.unique       = unique
    self.materialized = materialized
    self.packed       = packed

  @staticmethod
  def from_streams(genos, format, mergefunc=None):
    '''
    Combine multiple genostreams into one genomatrix stream

    @param     genos: genostreams
    @type      genos: list
    @param    format: format of input genomatrix, either 'ldat' or 'sdat'
    @type     format: str
    @param mergefunc: function to merge multiple genotypes into a consensus genotype
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
      raise ValueError, "Invalid genomatrix format '%s'.  Must be either sdat or ldat" % format

    formats  = set(g.format for g in genos)
    informat = formats.pop() if len(formats)==1 else None
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
      return genos.as_ldat(mergefunc)
    else:
      return genos.as_sdat(mergefunc)

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
  def from_tuples(genos, format, samples=None, loci=None, unique=True, modelmap=None):
    '''
    Alternate constructor that builds a new GenomatrixStream object from a
    genotype matrix stream with genotypes in a string format.
    '''
    if format=='ldat':
      columns = samples
    elif format=='sdat':
      columns = loci
    else:
      raise ValueError('Invalid genotype matrix format')

    columns,models,genos = encode_genomatrixstream_from_tuples(columns,genos,format,modelmap=modelmap)
    return GenomatrixStream(genos, format, samples=samples, loci=loci, models=models,
                                   unique=unique, packed=True, materialized=False)

  @staticmethod
  def from_strings(genos, format, genorepr, samples=None, loci=None, unique=True,
               packed=False, modelmap=None):
    '''
    Alternate constructor that builds a new GenomatrixStream object from a
    genotype matrix stream with genotypes in a string format.
    '''
    if format=='ldat':
      columns = samples
    elif format=='sdat':
      columns = loci
    else:
      raise ValueError('Invalid genotype matrix format')

    columns,models,genos = encode_genomatrixstream_from_strings(columns,genos,format,genorepr,modelmap=modelmap)
    return GenomatrixStream(genos, format, samples=samples, loci=loci, models=models,
                                   unique=unique, packed=True, materialized=False)

  def clone(self, genos, **kwargs):
    '''
    Alternative constructor that builds a new GenomatrixStream object
    with a new data stream and attributes based on self, but updated with
    specified keyword arguments.

    For example:

      stream.clone(genos, materialized=False)

    is equivalent to

      GenomatrixStream(genos, format=stream.format, samples=stream.samples,
                              loci=stream.loci, unique=stream.unique,
                              materialized=False, packed=stream.packed)
    '''
    kwargs.setdefault('format',       self.format)
    kwargs.setdefault('samples',      self.samples)
    kwargs.setdefault('loci',         self.loci)
    kwargs.setdefault('models',       self.models)
    kwargs.setdefault('unique',       self.unique)
    kwargs.setdefault('materialized', self.materialized)
    kwargs.setdefault('packed',       self.packed)

    return GenomatrixStream(genos, **kwargs)

  rows    = property(_get_rows,   _set_rows)
  columns = property(_get_columns,_set_columns)

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
    @param       mergefunc: function to merge multiple genotypes into a consensus genotype
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
    if transform is not None and kwargs:
      raise ValueError, 'specification of both a transform object and keyword arguments is not supported'
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

    if transform.rename_alleles is not None:
      genos = rename_genomatrixstream_alleles(genos,transform.rename_alleles)

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
    '''
    if self.unique:
      return self

    merge_rows = self.rows    is None or len(self.rows)    != len(set(self.rows))
    merge_cols = self.columns is None or len(self.columns) != len(set(self.columns))

    if not merge_rows and not merge_cols:
      return self

    if self.format == 'ldat':
      mergefunc = mergefunc_transpose_adapter(mergefunc)

    if merge_rows:
      genos = merge_genomatrixstream(self, mergefunc)
    else:
      genos = merge_genomatrixstream_columns(self, mergefunc)

    return genos

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
      format = 'sdat'
    else:
      format = 'ldat'

    return genos.clone(tgenos, format=format, packed=False, materialized=False)

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

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
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

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
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


def encode_genomatrixstream_from_tuples(columns, genos, format, modelmap=None):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix
  @type         genos: genomatrix generator
  @param     new_repr: internal representation of genotypes to be transformed to
  @type      new_repr: UnphasedMarkerRepresentation or similar object
  @return            : tuple of columns and a genomatrix generator in packed format
  @rtype             : 2-tuple of list of str and genomatrix generator

  >>> defmodel = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> modelmap = defaultdict(lambda: defmodel)

  With modelmap:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',modelmap)
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
  >>> loci,models,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  No modelmap:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat')
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
  >>> loci,models,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True
  '''
  if modelmap is None:
    modelmap = {}

  models = []

  if format=='ldat':
    def _load():
      n = len(columns)

      def _genokey(genos):
        gset = set(tuple(sorted(g)) for g in genos)
        gset.discard( (None,None) )
        return tuple(sorted(gset))

      modelcache = dict( (_genokey(m.genotypes),m) for m in models )
      descrcache = {}

      for label,row in genos:
        model = modelmap.get(label)
        if not model:
          key = _genokey(row)
          model = modelcache.get(key)
          if not model:
            model = modelcache[key] = UnphasedMarkerModel()
            for g in key:
              model.add_genotype(g)

        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        models.append(model)
        yield label,GenotypeArray(descr,row)

  elif format=='sdat':
    n = len(columns)
    updates = []

    for i,column in enumerate(columns):
      if column in modelmap:
        models.append(modelmap[column])
      else:
        model = UnphasedMarkerModel()
        models.append(model)
        updates.append( (i,model.add_genotype) )

    def _load():
      descr = GenotypeArrayDescriptor(models)

      if not updates:
        for label,row in genos:
          yield label,GenotypeArray(descr,row)
      else:
        for label,row in genos:
          for i,add in updates:
            add(row[i])
          yield label,GenotypeArray(descr,row)
  else:
    raise ValueError('Uknown format')

  return columns,models,_load()


def encode_genomatrixstream_from_strings(columns,genos,format,genorepr,modelmap=None):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix
  @type         genos: genomatrix generator

  >>> from reprs import snp
  >>> defmodel  = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> modelmap = defaultdict(lambda: defmodel)

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,modelmap)
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
  >>> loci,models,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp,modelmap)
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
  >>> samples,models,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp)
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
  >>> loci,models,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in new_rows:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])
  '''
  if modelmap is None:
    modelmap = {}

  models = []

  if format=='ldat':
    def _load():
      n = len(columns)

      modelcache  = {}
      descrcache  = {}
      strcache    = {}
      to_string    = genorepr.to_string
      from_strings = genorepr.from_strings

      for label,row in genos:
        model = modelmap.get(label)
        if not model:
          key = tuple(sorted(from_strings(set(row))))
          model = modelcache.get(key)
          if not model:
            model = modelcache[key] = UnphasedMarkerModel()
            for g in key:
              model.add_genotype(g)

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
          yield label,GenotypeArray(descr,imap(getitem, repeat(cache), row))
        except KeyError:
          gset = set(row)
          cache.update( (g,model[r]) for g,r in izip(gset,from_strings(gset)) )
          yield label,GenotypeArray(descr,imap(getitem, repeat(cache), row))

  elif format=='sdat':
    def _load():
      n = len(columns)
      updates   = []
      cachemap  = {}
      cachelist = []

      to_string    = genorepr.to_string
      from_strings = genorepr.from_strings

      for i,column in enumerate(columns):
        if column in modelmap:
          model = modelmap[column]
          models.append(model)
          update = model.get_genotype
        else:
          model = UnphasedMarkerModel()
          update = model.add_genotype
          models.append(model)

        cache = cachemap.get(model)
        if cache is None:
          cache = cachemap[model] = dict( (to_string(g),g) for g in model.genotypes )
          for g in genorepr.missing_geno_strs:
            cache[g] = model[None,None]

        cachelist.append(cache)
        updates.append( (update,cache) )

      repr  = genorepr.from_string
      descr = GenotypeArrayDescriptor(models)

      for label,row in genos:
        try:
          row = GenotypeArray(descr,imap(getitem, cachelist, row) )
        except KeyError:
          for (update,cache),gstr,g in izip(updates,row,from_strings(row)):
            cache[gstr] = update(g)
          row = GenotypeArray(descr,imap(getitem, cachelist, row) )

        yield label,row

  else:
    raise ValueError('Uknown format')

  return columns,models,_load()


def encode_genotriples_from_tuples(triples,modelmap=None):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: genomatriple
  @type       triples: genomatriple generator
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> triples = [('s3', 'l1', ('G', 'G')),('s3', 'l2', ('A', 'A')),
  ...            ('s2', 'l3', ('G', 'T')),('s1', 'l1', ('T', 'T')),
  ...            ('s1', 'l1', ('G', 'G')),('s2', 'l2', ('A', 'A'))]
  >>> for row in encode_genotriples_from_tuples(triples):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  if modelmap is None:
    modelmap = {}

  for sample,locus,geno in triples:
    model = modelmap.get(locus)
    if not model:
      modelmap[locus] = model = UnphasedMarkerModel()
    yield sample,locus,model.add_genotype(geno)


def encode_genotriples_from_strings(triples,genorepr,modelmap=None):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: genomatriple
  @type       triples: genomatriple generator
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> from reprs import snp
  >>> triples = [('s3', 'l1', 'GG'),('s3', 'l2', 'AA'),
  ...            ('s2', 'l3', 'GT'),('s1', 'l1', 'TT'),
  ...            ('s1', 'l1', 'GG'),('s2', 'l2', 'AA')]
  >>> for row in encode_genotriples_from_strings(triples,snp):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  repr = genorepr.from_string
  if modelmap is None:
    modelmap = {}
  for sample,locus,geno in triples:
    model = modelmap.get(locus)
    if not model:
      modelmap[locus] = model = UnphasedMarkerModel()
    yield sample,locus,model.add_genotype(repr(geno))


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

  @param   triples: genotriple stream
  @type    triples: sequence
  @param     order: sort order, either 'samples', 'locus'
  @type      order: str
  @param maxincore: maximum number of triples to process in core (not currently honored)
  @type  maxincore: int or None
  @return         : samples, loci, sorted genotriples
  @rtype          : tuple of list, list, genotriple sequence

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
    raise ValueError, "Unknown ordering specified: '%s'" % order

  # In memory sort should eventually be augmented by an offline multiphase
  # merge sort for very large sets.
  striples = sorted(triples,key=keyfunc)

  # Extract sample and locus sets since this can be done quickly with materialized lists
  samples = list(unique(imap(itemgetter(0), striples)))
  loci    = list(unique(imap(itemgetter(1), striples)))

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
  >>> triplelist = [GenotripleStream(trip1),GenotripleStream(trip2)]
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
    raise TypeError, 'cannot combine a single genotriple stream'
  elif len(triplelist) == 1:
    return triplelist[0]

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
  return GenotripleStream(triples, samples=samples, loci=loci, order=None, unique=unique)


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
  >>> trip1 = GenotripleStream.from_tuples(trip1)
  >>> trip2 = GenotripleStream.from_tuples(trip2)
  >>> triplelist = [GenotripleStream(trip1,order='sample'),GenotripleStream(trip2,order='sample')]
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
    raise TypeError, 'cannot combine a single genotriple stream'
  elif len(triplelist) == 1:
    return triplelist[0]

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
    raise ValueError, 'Cannot merge triplestreams in disparate orders'

  order = order.pop()

  if order == 'sample':
    key = itemgetter(0,1)
  elif order == 'locus':
    key = itemgetter(1,0)
  else:
    raise ValueError, 'Invalid order specified'

  # Perform the merge, marking all of the inputs as used
  triples = imerge(triplelist,key=key)

  # Return a new baby triplestream object
  return GenotripleStream(triples, samples=samples, loci=loci, order=order, unique=unique)


#######################################################################################


def merge_sorted_genotriples(triples,mergefunc):
  '''
  Merge genotypes from a genotriple list sorted by both sample and locus
  (although it does not matter which is primary and secondary key).  The
  supplied merge function is used to produce a consensus genotype for each
  sample and locus.

  @param   triples: sorted genotriple stream
  @type    triples: sequence
  @param mergefunc: merge function taking sample, locus, and list of genotypes
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
    for (sample,locus),trips in groupby(triples, itemgetter(0,1)):
      yield sample,locus,mergefunc(sample,locus,map(itemgetter(2),trips))

  return triples.clone(_merge(),unique=True,materialized=False)


def merge_genomatrixstream_columns(genos, mergefunc):
  '''
  Merge genotypes in columns with identical labels from a genomatrix stream
  using the specified merge function.  Results are in an unpacked list
  representation that may need to be repacked.

  A ValueError will be raised if this function is applied to a genomatrix
  with non-unique rows, since many merge algorithms and merge statistics
  cannot be performed in a step-wise fashion (ie., on columns and rows
  separately).

  The supplied merge function will be called assuming samples are on the
  rows and loci on the columns, for the purpose of collecting merge
  statistics.  If this is not the case, please use the
  glu.lib.genomerge.mergefunc_transpose_adapter.

  @param    columns: matrix column names
  @type     columns: sequence of strs
  @param      genos: genomatrix stream
  @type       genos: sequence
  @param  mergefunc: merge function taking sample, locus, and list of genotypes
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  >>> loci = ('l1','l2','l3')
  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'G') ]),
  ...         ('s2',[(None, None), (None, None), (None, None)]),
  ...         ('s3',[ ('A', 'A'),  (None, None), (None, None)]),
  ...         ('s4',[ ('A', 'T'),  (None, None),  ('T', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> geons = merge_genomatrixstream_columns(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('G', 'G')])
  ('s2', [(None, None), (None, None), (None, None)])
  ('s3', [('A', 'A'), (None, None), (None, None)])
  ('s4', [('A', 'T'), (None, None), ('T', 'T')])

  >>> loci = ('l1','l2','l1')
  >>> rows = [('s1',[('A', 'A'),   (None, None),  ('G', 'G')]),
  ...         ('s2',[ (None, None),  ('A', 'C'),   (None, None)]),
  ...         ('s3',[('A', 'A'),  ('A', 'A'),   (None, None)]),
  ...         ('s4',[('G', 'A'),   (None, None),  ('G', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> genos = merge_genomatrixstream_columns(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), (None, None)])
  ('s2', [(None, None), ('A', 'C')])
  ('s3', [('A', 'A'), ('A', 'A')])
  ('s4', [(None, None), (None, None)])
  '''
  assert mergefunc is not None

  merge_indices = defaultdict(list)
  merge_models  = defaultdict(list)
  new_columns   = []

  for i,column in enumerate(genos.columns):
    if column not in merge_indices:
      new_columns.append(column)
    merge_indices[column].append(i)
  new_columns = tuple(new_columns)

  # Trivial path: no merge is needed
  if new_columns == genos.columns:
    return genos

  # Non-trivial path: one or more columns must be merged
  def _merger():
    rows_seen = set()
    for row_label,row in genos:
      if row_label in rows_seen:
        raise ValueError('row labels required to be unique')

      rows_seen.add(row_label)

      merge_columns = ([row[i] for i in merge_indices[col_label]] for col_label in new_columns)
      new_row = list(imap(mergefunc,repeat(row_label),genos.columns,merge_columns))

      yield row_label,new_row

  if genos.format=='ldat':
    new_genos = genos.clone(_merger(),samples=new_columns,packed=False,materialized=False,unique=True)
  else:
    new_genos = genos.clone(_merger(),loci=new_columns,packed=False,materialized=False,unique=True)

  return new_genos


def merge_genomatrixstream_rows(genos, mergefunc):
  '''
  Merge genotypes in rows with identical labels from a genomatrix stream
  using the specified merge function.  The algorithm used requires a full
  materialization of the genomatrix stream, since row labels must be known.
  Results are in an unpacked list representation that may need to be
  repacked.

  A ValueError will be raised if this function is applied to a genomatrix
  with non-unique columns, since many merge algorithms and merge statistics
  cannot be performed in a step-wise fashion (ie., on columns and rows
  separately).

  The supplied merge function will be called assuming samples are on the
  rows and loci on the columns, for the purpose of collecting merge
  statistics.  If this is not the case, please use the
  glu.lib.genomerge.mergefunc_transpose_adapter.

  @param    columns: matrix column names
  @type     columns: sequence of strs
  @param      genos: genomatrix stream
  @type       genos: sequence
  @param  mergefunc: merge function taking sample, locus, and list ofgenotypes
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  >>> loci = ('l1','l2','l3')
  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'G')]),
  ...         ('s2',[(None, None), (None, None), (None, None)]),
  ...         ('s3',[ ('A', 'A'),  (None, None), (None, None)]),
  ...         ('s4',[ ('A', 'T'),  (None, None),  ('T', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> genos = merge_genomatrixstream_rows(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('G', 'G')])
  ('s2', [(None, None), (None, None), (None, None)])
  ('s3', [('A', 'A'), (None, None), (None, None)])
  ('s4', [('A', 'T'), (None, None), ('T', 'T')])

  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'T')]),
  ...         ('s2',[(None, None),  ('A', 'C'),  (None, None)]),
  ...         ('s1',[ ('A', 'A'),   ('A', 'A'),  (None, None)]),
  ...         ('s1',[ ('A', 'T'),  (None, None),  ('G', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)
  >>> genos = merge_genomatrixstream_rows(genos,VoteMerger())
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), ('A', 'C'), (None, None)])
  '''
  assert mergefunc is not None

  non_unique = any( n>1 for n in tally(genos.columns).itervalues() )

  if non_unique:
    raise ValueError('column labels required to be unique')

  merge_rows = defaultdict(list)
  new_rows   = []

  for row_label,row in genos:
    if row_label not in merge_rows:
      new_rows.append(row_label)
    merge_rows[row_label].append(row)

  def _merger():
    for row_label in new_rows:
      merge_columns = izip(*merge_rows.pop(row_label))
      new_row = list(imap(mergefunc,repeat(row_label),genos.columns,merge_columns))
      yield row_label,new_row

  if genos.format=='ldat':
    new_genos = genos.clone(_merger(),loci=new_rows,packed=False,materialized=False,unique=True)
  else:
    new_genos = genos.clone(_merger(),samples=new_rows,packed=False,materialized=False,unique=True)

  return new_genos


def merge_genomatrixstream(genos, mergefunc):
  '''
  Merge genotypes with identical row or column labels from a genomatrix
  stream using the specified merge function.  The algorithm used requires a
  full materialization of the genomatrix stream, since row labels must be
  known.  Results are in an unpacked list representation that may need to be
  repacked.

  The supplied merge function will be called assuming samples are on the
  rows and loci on the columns, for the purpose of collecting merge
  statistics.  If this is not the case, please use the
  glu.lib.genomerge.mergefunc_transpose_adapter.

  @param    columns: matrix column names
  @type     columns: sequence of strs
  @param      genos: genomatrix stream
  @type       genos: sequence
  @param  mergefunc: merge function taking sample, locus, and list of genotypes
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  >>> loci = ('l1','l2','l3')
  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'G')]),
  ...         ('s2',[(None, None), (None, None), (None, None)]),
  ...         ('s3',[ ('A', 'A'),  (None, None), (None, None)]),
  ...         ('s4',[ ('A', 'T'),  (None, None),  ('T', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)

  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('G', 'G')])
  ('s2', [(None, None), (None, None), (None, None)])
  ('s3', [('A', 'A'), (None, None), (None, None)])
  ('s4', [('A', 'T'), (None, None), ('T', 'T')])
  >>> sorted(merger.samplestats.iteritems())

  >>> sorted(merger.locusstats.iteritems())

  >>> rows = [('s1',[ ('A', 'A'),  (None, None),  ('G', 'T') ]),
  ...         ('s2',[(None, None), ('A', 'C'),   (None, None)]),
  ...         ('s1',[ ('A', 'A'),  ('A', 'A'),   (None, None)]),
  ...         ('s1',[ ('A', 'T'), (None, None),  ('G', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)

  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2', 'l3')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), ('A', 'C'), (None, None)])
  >>> sorted(merger.samplestats.iteritems())

  >>> sorted(merger.locusstats.iteritems())

  >>> loci = ('l1','l2','l1')
  >>> rows = [('s1',[(None, None), (None, None),  ('C', 'T')]),
  ...         ('s2',[(None, None),  ('A', 'G'),   ('T', 'T')]),
  ...         ('s1',[ ('C', 'C'),   ('A', 'A'),  (None, None)]),
  ...         ('s1',[ ('C', 'T'),  (None, None),  ('C', 'T')])]
  >>> genos = GenomatrixStream.from_tuples(rows,'sdat',loci=loci)

  >>> merger=VoteMerger()
  >>> genos = merge_genomatrixstream(genos,merger)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [(None, None), ('A', 'A')])
  ('s2', [('T', 'T'), ('A', 'G')])
  >>> sorted(merger.samplestats.iteritems())

  >>> sorted(merger.locusstats.iteritems())

  '''
  assert mergefunc is not None

  merge_indices = defaultdict(list)

  columns = genos.columns
  new_columns = []
  for i,column in enumerate(columns):
    if column not in merge_indices:
      new_columns.append(column)
    merge_indices[column].append(i)

  new_columns = tuple(new_columns)

  # FIXME: Check if merging is needed

  # Repack, since we're about to materialize
  genos = genos.transformed(repack=True)

  new_rows   = []

  def _merger():
    merge_rows = defaultdict(list)

    for row_label,row in genos:
      if row_label not in merge_rows:
        new_rows.append(row_label)
      merge_rows[row_label].append(row)

    if new_columns == columns:
      # Optimize case where columns are unique as in merge_genomatrix_rows
      for row_label in new_rows:
        # Form list of all genotypes at each new column, freeing used rows
        merge_columns = izip(*merge_rows.pop(row_label))

        # Merge genotypes
        new_row = list(imap(mergefunc,repeat(row_label),columns,merge_columns))

        # Yield new row
        yield row_label,new_row

    else:
      # Otherwise, apply fully general merge over duplicate rows and columns
      for row_label in new_rows:
        # Form list of all genotypes at each new column
        merge_columns = ([ row[i] for row in merge_rows[row_label] for i in merge_indices[col_label] ]
                                  for col_label in new_columns)

        # Merge genotypes
        new_row = list(imap(mergefunc,repeat(row_label),new_columns,merge_columns))

        # Free used rows
        del merge_rows[row_label]

        # Yield new row
        yield row_label,new_row

  if genos.format=='ldat':
    new_genos = genos.clone(_merger(),samples=new_columns,loci=new_rows,packed=False,materialized=False,unique=True)
  else:
    new_genos = genos.clone(_merger(),loci=new_columns,samples=new_rows,packed=False,materialized=False,unique=True)

  return new_genos


def merge_genomatrixstream_list(genos, mergefunc):
  '''
  Take a sequence of genotype matrix streams and merge all genotypes
  identical row or column labels into a single genotype matrix stream using
  the specified merge function.

  All input matricies must meet the following requirements:
    1) share the same orientation (either sdat or ldat)
    2) share the same genotype representation

  The algorithm used requires a full materialization of the genomatrix
  streams, since row labels must be known.  Results are in an unpacked list
  representation that may need to be repacked.

  The supplied merge function will be called assuming samples are on the
  rows and loci on the columns, for the purpose of collecting merge
  statistics.  If this is not the case, please use the
  glu.lib.genomerge.mergefunc_transpose_adapter.

  @param    columns: list matrix column names
  @type     columns: list sequence of strs
  @param      genos: sequence of genomatrix streams
  @type       genos: sequence
  @param  mergefunc: merge function taking sample, locus, and list of genotypes
  @type   mergefunc: callable
  @return          : merged unpacked genomatix stream
  @rtype           : generator

  Test slow-path for heterogeneous schema:

  >>> samples1 =        ('s1',       's2',       's3')
  >>> rows1 = [('l1',[('G', 'G'),   ('A', 'A'),  ('A', 'G')]),
  ...          ('l2',[('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> samples2 =         ('s1',       's3',        's4')
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples1).materialize()

  >>> rows2 = [('l1',[(None, None), ('A', 'G'),  ('A', 'G')]),
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

  Test fast-path for homogeneous schema:

  >>> samples =          ('s1',         's2',        's3')
  >>> rows1 = [('l1',[ ('G', 'G'),   ('G', 'T'),  ('T', 'T')]),
  ...          ('l2',[ ('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> rows2 = [('l1',[(None, None),  ('A', 'T'),  ('T', 'T')]),
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
  '''
  assert mergefunc is not None

  if not genos:
    raise ValueError('Invalid empty genotype list')

  # Fastpath for identical schema: Concatenate all genotype rows and merge
  # using the single-matrix merge function.  Also optimizes the single input
  # case.
  format = genos[0].format
  if not all(g.format==format for g in genos):
    raise ValueError('Input genotypes must all be in same format')

  columns = [ g.columns for g in genos ]
  if all(columns[0]==c for c in columns):
    # Pass-through from merge_genomatrix
    if format=='ldat':
      genos = genos[0].clone(chain(*genos),loci=None,packed=False,materialized=False,unique=False)
    else:
      genos = genos[0].clone(chain(*genos),samples=None,packed=False,materialized=False,unique=False)

    return merge_genomatrixstream(genos, mergefunc)

  new_rows      = {}
  new_columns   = {}
  merge_columns = defaultdict(list)
  merge_rows    = defaultdict(lambda: defaultdict(list))

  for i,g in enumerate(genos):
    # Collect columns and form mappings from old schemas to new
    for j,column in enumerate(g.columns):
      k=new_columns.setdefault(column,len(new_columns))
      merge_columns[i].append( (j,k) )

    # Collect row labels and materialize all rows for later merging
    for row_label,row in g:
      new_rows.setdefault(row_label,len(new_rows))
      merge_rows[row_label][i].append(row)

  # Invert row and column dictionaries to recover insertion orderings
  new_columns = tuple(imap(itemgetter(0),sorted(new_columns.iteritems(),key=itemgetter(1))))
  new_rows    = tuple(imap(itemgetter(0),sorted(new_rows.iteritems(),   key=itemgetter(1))))
  n = len(new_columns)

  def _merger():
    # Fully general merge over duplicate rows and columns
    for row_label in new_rows:
      # Form lists of all genotypes at each new column for all rows with row_label
      new_row = [ [] for i in xrange(n) ]

      # Iterate over input rows and schema, find the cooresponding column
      # mappings, and append the relevant genotypes
      for i,rows in merge_rows.pop(row_label).iteritems():
        for j,k in merge_columns.get(i,[]):
           new_row[k] += (row[j] for row in rows)

      # Merge genotypes
      new_row = list(imap(mergefunc,repeat(row_label),new_columns,new_row))

      # Yield new row
      yield row_label,new_row

  if format=='ldat':
    new_genos = genos[0].clone(_merger(),samples=new_columns,loci=new_rows,packed=False,materialized=False,unique=True)
  else:
    new_genos = genos[0].clone(_merger(),loci=new_columns,samples=new_rows,packed=False,materialized=False,unique=True)

  return new_genos


#######################################################################################


def build_genomatrixstream_from_genotriples(triples, format, mergefunc):
  '''
  Build genomatrix from genotriples using either the xtab or the rowsby
  function.  The rowsby function would be chosen over xtab if and only if:

    1. the necessary columns are given, and
    2. triples have been ordered appropriately (specified by the order argument).

  @param     genos: genotriples
  @type      genos: sequence
  @param    format: format genomatrix
  @type     format: string
  @param mergefunc: function to merge multiple genotypes into a consensus genotype
  @type  mergefunc: callable
  @param   samples: optional sequence of samples if known, otherwise None
  @type    samples: sequence of str or None
  @param      loci: optional sequence of samples if known, otherwise None
  @type       loci: sequence of str or None
  @param     order: order of genotriples ('sample' or 'locus')
  @type      order: str
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
  >>> merge = VoteMerger()
  >>> triples.samples = None
  >>> triples.loci=['l1','l2']
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merge)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A')])
  ('s2', [('G', 'T'), ('T', 'T')])
  ('s3', [('G', 'G'), ('A', 'A')])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merge = VoteMerger()
  >>> triples.loci=None
  >>> triples.samples=['s1','s2','s3']
  >>> genos = build_genomatrixstream_from_genotriples(triples,'sdat',merge)
  >>> genos.loci
  ('l1', 'l2')
  >>> for row in genos:
  ...   print row
  ('s1', [('G', 'G'), ('A', 'A')])
  ('s2', [('G', 'T'), ('T', 'T')])
  ('s3', [('G', 'G'), ('A', 'A')])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merge = VoteMerger()
  >>> genos = build_genomatrixstream_from_genotriples(triples,'ldat',merge)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
  ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]
  '''
  columns = None
  if format == 'sdat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(0),itemgetter(1),itemgetter(2)
    columns = tuple(sorted(triples.loci)) if triples.loci else None
  elif format == 'ldat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(1),itemgetter(0),itemgetter(2)
    mergefunc = mergefunc_transpose_adapter(mergefunc)
    columns = tuple(sorted(triples.samples)) if triples.samples else None
  else:
    raise NotImplementedError,'triple to %s format conversion is not supported' % format

  order = triples.order
  if format == 'sdat' and order != 'sample':
    order = False
  elif format == 'ldat' and order != 'locus':
    order = False

  if columns is None or not order:
    columns,rows,data = xtab(triples, rowkeyfunc, colkeyfunc, valuefunc, mergefunc)
    genos = tuple(izip(rows,data))
    rows = tuple(rows)
  else:
    columns,genos = rowsby(triples, columns, rowkeyfunc, colkeyfunc, valuefunc, mergefunc)
    # FIXME: One day we should be able to trust triples.rows
    rows = None

  columns = tuple(columns)
  if format=='ldat':
    new_genos = GenomatrixStream(genos,format,samples=columns,loci=rows,unique=True)
  else:
    new_genos = GenomatrixStream(genos,format,loci=columns,samples=rows,unique=True)

  return new_genos


#######################################################################################


def pack_genomatrixstream(genos,modelmap=None):
  '''
  Transform a genomatrix into an internal packed representation

  @param    columns: matrix column names
  @type     columns: sequence of strs
  @param      genos: genomatrix
  @type       genos: genomatrix generator
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
      genos = iter(genos)
      try:
        label,row = genos.next()
      except StopIteration:
        return

      descr = GenotypeArrayDescriptor([g.model for g in row])
      yield label,GenotypeArray(descr,row)

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

  return genos.clone(_pack(genos),packed=True)


def rename_genomatrixstream_alleles(genos, rename_alleles):
  '''
  Returns a new genomatrix with the alleles renamed

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix
  @type         genos: genomatrix generator
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
  # FIXME: There is likely a better way that doesn't involve reencoding from scratch
  def rename():
    if genos.format=='ldat':
      for label,row in genos:
        if label in rename_alleles:
          r   = rename_alleles[label]
          row = [ ((r[g[0]],r[g[1]]) if g else g.alleles()) for g in row ]
        yield label,row

    elif genos.format=='sdat':
      remaps = [ rename_alleles.get(h) for h in genos.columns ]
      for label,row in genos:
        row = [ ((r[g[0]],r[g[1]]) if g and r else g.alleles()) for g,r in izip(row,remaps) ]
        yield label,row
    else:
      raise ValueError('Matrix format must be specified when renaming alleles')

  return GenomatrixStream.from_tuples(rename(),genos.format,samples=genos.samples,loci=genos.loci,
                                      unique=genos.unique)


def rename_genotriples_alleles(triples, rename_alleles):
  '''
  Returns a new genotriple stream with the alleles renamed

  @param      triples: genomatriple
  @type       triples: genomatriple generator
  @return            : genotriple with the new internal format
  @rtype             : genotriple generator

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
  # FIXME: There is likely a better way that doesn't involve reencoding from scratch
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

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param   genos: genomatrix
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
  def _filter():
    for lname,row in genos:
      if any(row):
        yield lname,row

  if genos.format=='ldat':
    new_genos = genos.clone(_filter(),loci=None,materialized=False)
  else:
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
  columns_notseen = pick(new_genos.columns, sorted(columns_notseen))
  new_genos = new_genos.clone(rows)
  return filter_genomatrixstream_by_column(new_genos,columns_notseen,exclude=True)


def filter_genotriples_missing(triples):
  '''
  Filter out the genotriple if its genotype is missing.

  @param triples: sequnce of triplets. e.g.
                  [('locus1','sampleA',genotype),
                   ('locus2','sampleB',genotype),...]
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
  @param rows: genotype matrix
  @type  rows: sequence
  @rtype:      generator
  @returns:    a genotype triplet stream

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
      for locus,row in genos:
        for sample,geno in izip(samples, row):
          yield sample,locus,geno
  else:
    order = 'sample'
    def _build():
      loci = genos.loci
      for sample,row in genos:
        for locus,geno in izip(loci, row):
          yield sample,locus,geno

  # Unset the result order if any of the required ordering constraints
  # cannot be verified
  if (genos.samples is None or genos.loci is None
                           or sorted(genos.samples) != list(genos.samples)
                           or sorted(genos.loci)    != list(genos.loci)):
    order = None

  return GenotripleStream(_build(), unique=genos.unique, samples=genos.samples, loci=genos.loci, order=order)


def rename_genotriples(triples,samplemap,locusmap):
  '''
  Rename the sample and locus for each genotriple according to the samplemap
  and locusmap. If there is not mapping for a particular sample or locus,
  the original name will be used.

  @param   triples: a sequence of genotriples. e.g.
                    ('s1','l1','AA'),...
  @type    triples: sequence
  @param samplemap: map between the current sample name
                    and a new name
  @type  samplemap: dict
  @param  locusmap: map between the current locus name
                    and a new name
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
  def _rename():
    for sample,locus,geno in triples:
      if samplemap is not None:
        sample = samplemap.get(sample,sample)
      if locusmap is not None:
        locus = locusmap.get(locus,locus)
      yield sample,locus,geno

  samples = triples.samples
  if samples and samplemap:
    samples = set(samplemap.get(s,s) for s in samples)

  loci = triples.loci
  if loci and locusmap:
    loci = set(locusmap.get(l,l) for l in loci)

  return triples.clone(_rename(),samples=samples,loci=loci,order=None,unique=False)


def filter_genotriples(triples,sampleset,locusset,exclude=False):
  '''
  Filter out the genotriples against the sampelset and locusset.
  Depending on the value of exclude flag, both sets will be treated
  as either inclusion sets(exclude=False) or exclusion sets(exclude=True).

  @param   triples: sequence of genotriples. e.g.
                    ('s1','l1','AA'),...
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
      samples = set(s for s in samples if l not in sampleset)
    else:
      samples = set(s for s in samples if l in sampleset)
  elif samples is None and sampleset is not None and not exclude:
    samples = sampleset

  return triples.clone(_filter(),loci=loci,samples=samples,order=None)


def remap_genotriples(triples,samplemap,locusmap):
  '''
  Remap the genotriples according to the samplemap and locusmap.
  This function will rename the sample and locus names in the genotriple
  and will also filter out the genotriple if either sample or locus is
  not in the samplemap or locusmap.

  @param   triples: sequence of genotriples. e.g.
                    ('s1','l1','AA'),...
  @type    triples: sequence
  @param samplemap: map between the two set of sample names
  @type  samplemap: dict
  @param  locusmap: map between the two set of locus names
  @type   locusmap: dict
  @return         : remapped genotriples
  @rtype          : generator
  '''
  triples = filter_genotriples(triples,samplemap,locusmap)
  return rename_genotriples(triples,samplemap,locusmap)


class NonUniqueError(ValueError): pass


def rename_genomatrixstream_column(genos,colmap):
  '''
  Rename the columns for the genotype matrix data
  according to a name mapping. If the name of the column
  is not in the mapping, its original name will be used.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  colmap: map of the column names
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
  '''
  columns = tuple(colmap.get(name,name) for name in genos.columns)
  if columns == genos.columns:
    return genos

  if genos.format=='ldat':
    genos = genos.clone(genos.use_stream(), samples=columns)
  else:
    genos = genos.clone(genos.use_stream(), loci=columns)

  return genos


def filter_genomatrixstream_by_column(genos,colset,exclude=False):
  '''
  Filter the genotype matrix data by a column set.  Depending on the
  value of the exclude flag, the column set will be used either for
  inclusion (exclude=False) or exclusion (exclude=True).

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  colset: set of the column names
  @type   colset: set
  @param exclude: flag to exclude colset instead of including
  @type  exclude: boolean
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

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> genos = GenomatrixStream.from_strings(rows,'ldat',snp,samples=samples)
  >>> genos = filter_genomatrixstream_by_column(genos,colset,exclude=True)
  >>> genos.samples
  ('s2',)
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'G')])
  ('l2', [('A', 'T')])
  '''
  columns = genos.columns
  if exclude:
    columns = tuple((name,i) for i,name in enumerate(columns) if name not in colset)
  else:
    columns = tuple((name,i) for i,name in enumerate(columns) if name     in colset)

  if columns:
    columns,indices = izip(*columns)
  else:
    indices = []

  if columns == genos.columns:
    return genos

  def _filter():
    for label,row in genos:
      yield label,pick(row[:],indices)

  if genos.format=='ldat':
    new_genos=genos.clone(_filter(),samples=columns,packed=False,materialized=False)
  else:
    models = genos.models
    if models is not None:
      models = pick(models,indices)
    new_genos=genos.clone(_filter(),loci=columns,models=models,packed=False,materialized=False)

  return new_genos


def remap_genomatrixstream_column(genos,colmap):
  '''
  Rename and filter column labels for a genotype matrix data.  If there is
  no mapping for a column label, then the original label will be used.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  colmap: map between the two set of names
  @type   colmap: dict
  @rtype:         generator
  @return:        genotype matrix with renamed and filtered column labels
  '''
  genos = filter_genomatrixstream_by_column(genos,colmap)
  return rename_genomatrixstream_column(genos,colmap)


def reorder_genomatrixstream_columns(genos,labels):
  '''
  Reorder and filter the genotype matrix columns to match the sequence of
  labels provided.  If not all labels appear, the remainder are retained at
  the end of the list in lexicographical order.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being column metadata
  @type     rows: sequence
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
    raise ValueError, 'Duplicated column label: %s' % l

  new_columns = pick(genos.columns,indices)
  rows = iter(genos)
  def _reorder():
    for label,row in rows:
      yield label,pick(row[:], indices)

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

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being column metadata
  @type     rows: sequence
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
    raise ValueError, 'Duplicated row label: %s' % l

  # FIXME: Update models
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


def rename_genomatrixstream_row(genos,rowmap):
  '''
  Rename the row label for the genotype matrix data.
  If there is no mapping for a row label, then the original
  label will be used.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  rowmap: map between the two set of names
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
  '''
  def _rename():
    for label,row in genos:
      yield rowmap.get(label,label),row

  if genos.rows is not None:
    rows = tuple(rowmap.get(label,label) for label in genos.rows)
  else:
    rows = None

  if genos.format=='ldat':
    new_genos = genos.clone(_rename(),loci=rows,materialized=False)
  else:
    new_genos = genos.clone(_rename(),samples=rows,materialized=False)

  return new_genos


def filter_genomatrixstream_by_row(genos,rowset,exclude=False):
  '''
  Filter the genotype matrix by a row set.
  Depending on the value of the exclude flag, the row set will
  be used either for inclusion(exclude=False) or exclusion(exclude=True)
  purpose.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  rowset: set of the row names
  @type   rowset: set
  @param exclude: flag to exclude items in rowset instead of including
  @type  exclude: boolean
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
  '''
  # FIXME: Implement fast path when rows are known
  # FIXME: Update models
  def _filter():
    if exclude:
      for label,row in genos:
        if label not in rowset:
          yield label,row
    else:
      for label,row in genos:
        if label in rowset:
          yield label,row

  rows = genos.rows
  if rows is not None:
    if exclude:
      rows = tuple(label for label in rows if label not in rowset)
    else:
      rows = tuple(label for label in rows if label     in rowset)

  if genos.format=='ldat':
    new_genos=genos.clone(_filter(),loci=rows,materialized=False)
  else:
    new_genos=genos.clone(_filter(),samples=rows,materialized=False)

  return new_genos


def remap_genomatrixstream_row(genos,rowmap):
  '''
  Remap the genotype matrix according to the rowmap.
  This function will rename the row labels in the genotype matrix
  and will also filter out the row if its label is not in the rowmap

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: genotype matrix data with the first row
                  being the column meta-data
  @type     rows: sequence
  @param  rowmap: map between the two set of names
  @type   rowmap: dict
  @return       : genotype matrix with remapped rows
  @rtype        : generator
  '''
  genos = filter_genomatrixstream_by_row(genos,rowmap)
  return rename_genomatrixstream_row(genos,rowmap)


# Does not need to support a genorepr, since only one row is
# materialized at a time.
def transpose_generator(columns, rows, missing=None):
  '''
  Transpose a matrix of row labels and row data generating one column of
  data at a time.  Return a generator of the columns of data stored in rows,
  yielding sucessive column labels and column data.

  Requires the input data to be fully materialized, but allows results to be
  streamed.

  @param columns: matrix column names
  @type  columns: sequence of strs
  @param    rows: iterable sequence of pairs of row labels and row data
  @type     rows: iterable of label and sequence pairs
  @param columns: sequence of column labels corresponding to each row
  @type  columns: sequence of labels
  @return       : tuple of column labels and generator of tuples of row labels and row data
  @rtype        : tuple

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
  m = len(columns)

  if not n or not m:
    return [],[]

  rowlabels,rows = zip(*rows)

  def _transpose_generator():
    for j,column in enumerate(columns):
      newrow = [missing]*n
      for i in xrange(n):
        newrow[i] = rows[i][j]
      yield column,newrow

  return rowlabels,_transpose_generator()


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
