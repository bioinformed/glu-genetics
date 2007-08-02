# -*- coding: utf-8 -*-
'''
File:          genodata.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU genotype data management objects

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from   types       import NoneType
from   operator    import itemgetter
from   collections import defaultdict
from   itertools   import izip,islice,ifilter,imap,chain,groupby,repeat

from   utils       import pick,tally,unique
from   fileutils   import load_list,load_map
from   imerge      import imerge
from   xtab        import xtab,rowsby
from   genoarray   import snp_acgt,snp_marker,get_genorepr
from   genomerge   import UniqueMerger, VoteMerger, mergefunc_transpose_adapter


#FIXME: This function is here to be used as a genorepr when native genotype
#       strings are required.  However, it may move to glu.lib.utils.
def intern_list(seq):
  '''
  A string list functor that automatically interns all elements passed to it

  >>> l  = ['a','b','c','d']
  >>> l1 = intern_list( c+c for c in l )   # Ensure strings are not shared
  >>> l2 = intern_list( c+c for c in l )   # local constants
  >>> all( s1 is s2 is intern(s1) for s1,s2 in izip(l1,l2) )
  True
  >>> l1 = list( c+c for c in l )          # Ensure strings are not shared
  >>> l2 = list( c+c for c in l )          # local constants
  >>> all( s1 is s2 is intern(s1) for s1,s2 in izip(l1,l2) )
  False
  '''
  return list(imap(intern,seq))


#######################################################################################


class GenotypeStream(object):
  __slots__ = []


class GenotripleStream(GenotypeStream):
  '''
  A stream of genotriples with optional metadata
  '''
  format = 'genotriple'

  def __init__(self, triples, samples=None, loci=None, order=None, unique=False, genorepr=snp_acgt, materialized=False):
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
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    @param materialized: flag indicating if genos is a materialized
                         (multiply iterable) data type
    @type  materialized: bool
    '''
    if order not in (None,'locus','sample'):
      raise ValueError

    if isinstance(triples, (list,tuple)):
      materialized = True

    if not isinstance(samples,(NoneType,set)):
      samples = set(samples)

    if not isinstance(loci,(NoneType,set)):
      loci = set(loci)

    self.samples      = samples
    self.loci         = loci
    self.order        = order
    self.unique       = bool(unique)
    self.materialized = materialized
    self.genorepr     = genorepr

    self.__private_triples_do_not_touch = triples

  @staticmethod
  def from_streams(genos, mergefunc=None, order=None, genorepr=snp_acgt):
    '''
    Combine multiple genostreams into one genotriple stream

    @param     genos: genostreams
    @type      genos: list
    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @param     order: ordering of the triple stream, 'sample' or 'locus', or None
    @type      order: str or None
    @param  genorepr: object representing the input/output encoding and
                      internal representation of genotypes
    @type   genorepr: UnphasedMarkerRepresentation or similar object
    @return         : combined genotriple stream
    @rtype          : sequence of sample, locus, and genotype

    >>> trip1 = [('s1','l1', 51),('s1','l2', 17),
    ...          ('s2','l2', 68),('s3','l1', 68)]
    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ 51,   52,   68 ]),
    ...         ('l2', [ 17,   68,   51 ])]
    >>> trip2 = [('s2','l1', 51),('s3','l2', 17)]
    >>> streams = [GenotripleStream(trip1),
    ...            GenomatrixStream(rows,'ldat',samples=samples),
    ...            GenotripleStream(trip2)]
    >>> combined = GenotripleStream.from_streams(streams, genorepr=snp_marker, mergefunc=VoteMerger())
    >>> for row in combined:
    ...   print row
    ('s1', 'l1', ('G', 'G'))
    ('s1', 'l2', ('A', 'A'))
    ('s2', 'l1', 0)
    ('s2', 'l2', ('T', 'T'))
    ('s3', 'l1', ('T', 'T'))
    ('s3', 'l2', 0)
    '''
    if not genos:
      raise ValueError

    if genorepr is None:
      genorepr = genos[0].genorepr

    if mergefunc is not None and order is None:
      order = 'sample'

    triples = [ g.as_genotriples().transformed(genorepr=genorepr) for g in genos ]

    if order is not None:
      triples = [ t.sorted(order) for t in triples ]

    if len(triples) == 1:
      triples = triples[0]

    elif order is None and mergefunc is None:
      triples = combine_unsorted_genotriple_list(triples)

    else:
      triples = combine_sorted_genotriple_list(triples)

    if mergefunc is not None:
      triples = triples.merged(mergefunc)

    assert isinstance(triples,GenotripleStream)

    return triples

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
    return self.__private_triples_do_not_touch is None

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
      self.__private_triples_do_not_touch,triples = None,self.__private_triples_do_not_touch
    else:
      triples = self.__private_triples_do_not_touch

    return triples

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

    return GenotripleStream(list(self), samples=self.samples, loci=self.loci, order=self.order,
                                        genorepr=self.genorepr, materialized=True)

  def recoded(self, genorepr):
    '''
    Returns a new genotriple stream with genotypes recoded to the new representation

    @param genorepr: object representing the input/output encoding and
                     internal representation of genotypes
    @type  genorepr: UnphasedMarkerRepresentation or similar object
    @return        : recoded genotriple stream
    @rtype         : GenotripleStream

    >>> triples = GenotripleStream([('s1','l1',51),('s1','l2',17),('s2','l2',68),('s3','l1',68)])
    >>> recoded = list(triples.recoded(snp_marker))
    >>> recoded
    [('s1', 'l1', ('G', 'G')), ('s1', 'l2', ('A', 'A')), ('s2', 'l2', ('T', 'T')), ('s3', 'l1', ('T', 'T'))]
    >>> list(GenotripleStream(recoded,genorepr=snp_marker).recoded(snp_acgt))
    [('s1', 'l1', 51), ('s1', 'l2', 17), ('s2', 'l2', 68), ('s3', 'l1', 68)]
    '''
    return self.transformed(genorepr=genorepr)

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
    @param        genorepr: object representing the input/output encoding and
                            internal representation of genotypes
    @type         genorepr: UnphasedMarkerRepresentation or similar object
    @return               : transformed genotriple stream
    @rtype                : GenotripleStream

    >>> trips = [('s1','l1', 51),('s1','l2', 17),
    ...          ('s2','l1', 52),('s2','l2', 68),
    ...          ('s3','l1', 51),('s3','l2', 17)]
    >>> for row in GenotripleStream(trips).transformed(include_loci=['l1'],exclude_samples=['s3']):
    ...   print row
    ('s1', 'l1', 51)
    ('s2', 'l1', 52)
    '''

    if transform is not None and kwargs:
      raise ValueError, 'specification of both a transform object and keyword arguments is not supported'
    elif kwargs:
      transform = GenoTransform.from_kwargs(**kwargs)

    if transform is None:
      return self

    # Build a transform object from kwargs if one was not specified
    if transform is None:
      transform = GenoTransform.from_kwargs(**kwargs)

    triples = self.use_stream()
    samples = self.samples
    loci    = self.loci

    # Filter missing genotyes
    if transform.filter_missing_genotypes:
      triples = filter_genotriples_missing(triples)

    # Sample and locus includes
    if transform.samples.include is not None or transform.loci.include is not None:
      if samples is not None and transform.samples.include is not None:
        samples = set(s for s in samples if s in transform.samples.include)
      elif samples is None and transform.samples.include is not None:
        samples = transform.samples.include

      if loci is not None and transform.loci.include is not None:
        loci = set(l for l in loci if l in transform.loci.include)
      elif loci is None and transform.loci.include is not None:
        loci = transform.loci.include

      triples = filter_genotriples(triples,transform.samples.include,transform.loci.include)

    # Sample and locus excludes
    if transform.samples.exclude is not None or transform.loci.exclude is not None:
      if samples is not None and transform.samples.exclude is not None:
        samples = set(s for s in samples if s not in transform.samples.exclude)
      if loci is not None and transform.loci.exclude is not None:
        loci = set(l for l in loci if l not in transform.loci.exclude)

      triples = filter_genotriples(triples,transform.samples.exclude,transform.loci.exclude,exclude=True)

    # Determine if resulting triples will be unique (before renaming samples and loci)
    unique_results = prove_unique_transform(transform=transform,loci=loci,samples=samples,unique=self.unique)

    # Sample and locus renaming
    if transform.samples.rename is not None or transform.loci.rename is not None:
      if samples is not None and transform.samples.rename is not None:
        samples = set(transform.samples.rename.get(s,s) for s in samples)
      if loci is not None and transform.loci.rename is not None:
        loci = set(transform.loci.rename.get(l,l) for l in loci)

      triples = rename_genotriples(triples,transform.samples.rename,transform.loci.rename)

    # Determine result order
    if transform and (transform.samples.rename or transform.loci.rename):
      result_order = None
    else:
      result_order = self.order

    genorepr = self.genorepr
    if transform.genorepr is not None:
      genorepr = transform.genorepr
      if isinstance(genorepr, basestring):
        genorepr = get_genorepr(genorepr)

    if genorepr is not self.genorepr or transform.rename_alleles:
      triples = recode_genotriples(triples, self.genorepr, genorepr, transform.rename_alleles)

    # Build the new triple stream
    triples = GenotripleStream(triples, samples=samples, loci=loci, order=result_order, unique=unique_results,
                               genorepr=genorepr)

    if transform.samples.order is not None and transform.loci.order is not None:
      if transform.samples.order is not None:
        order = 'samples'
      else:
        order = 'loci'
      triples = triples.sorted('samples',sampleorder=transform.samples.order,locusorder=transform.loci.order)

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

    >>> trips = [('s3','l1', 51),('s1','l2', 17),
    ...          ('s2','l1', 52),('s3','l2', 68),
    ...          ('s1','l1', 51),('s2','l2', 17)]
    >>> for row in GenotripleStream(trips).sorted():
    ...   print row
    ('s1', 'l1', 51)
    ('s1', 'l2', 17)
    ('s2', 'l1', 52)
    ('s2', 'l2', 17)
    ('s3', 'l1', 51)
    ('s3', 'l2', 68)
    '''
    if self.order is not None and self.order == order:
      return self

    samples,loci,triples = sort_genotriples(self,order=order,locusorder=locusorder,sampleorder=sampleorder)

    return GenotripleStream(triples, samples=samples, loci=loci, order=order, unique=self.unique,
                                     genorepr=self.genorepr)


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

    >>> trips = [('s1','l1', 68),('s1','l2', 17),
    ...          ('s2','l1', 52),('s2','l2', 17),
    ...          ('s1','l1', 51),('s2','l2', 17)]
    >>> for row in GenotripleStream(trips).merged(VoteMerger(),order='locus'):
    ...   print row
    ('s1', 'l1', 0)
    ('s2', 'l1', 52)
    ('s1', 'l2', 17)
    ('s2', 'l2', 17)
    '''
    if self.unique:
      return self

    if self.order not in ('sample','locus'):
      return self.sorted(order).merged(mergefunc)

    triples = merge_sorted_genotriples(self,mergefunc)

    return GenotripleStream(triples, samples=self.samples, loci=self.loci, order=self.order, unique=True,
                                     genorepr=self.genorepr)

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

    >>> trips = [('s1','l1',51),('s1','l2',17),
    ...          ('s2','l1',52),('s2','l2',68),
    ...          ('s3','l1',51),('s3','l2',17)]
    >>> samples =         ('s1',       's2',       's3')
    >>> rows = [('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')]),
    ...         ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])]
    >>> merge = VoteMerger(missingrepr=snp_marker.missing)
    >>> ldat = GenotripleStream(trips).as_ldat(merge).transformed(genorepr=snp_marker)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = GenotripleStream(trips,loci=['l1','l2']).as_ldat(merge).transformed(genorepr=snp_marker)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])

    >>> ldat = GenotripleStream(trips,samples=['s1','s2','s3']).as_ldat(merge).transformed(genorepr=snp_marker)
    >>> ldat.samples
    ('s1', 's2', 's3')
    >>> for row in ldat:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('G', 'G')])
    ('l2', [('A', 'A'), ('T', 'T'), ('A', 'A')])
    '''
    if mergefunc is None:
      mergefunc = UniqueMerger(missingrepr=self.genorepr.missing)

    samples,genos = build_genomatrixstream_from_triples(self, 'ldat', mergefunc=mergefunc,
                          samples=self.samples, loci=self.loci, order=self.order,
                          genorepr=self.genorepr)

    return GenomatrixStream(genos, 'ldat', samples=samples, unique=True, genorepr=self.genorepr)

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

    >>> trips = [('s1','l1', 51),('s1','l2', 17),
    ...          ('s2','l1', 52),('s2','l2', 68),
    ...          ('s3','l1', 51),('s3','l2', 17)]
    >>> loci = ('l1','l2')
    >>> rows = [('s1', [ 51,   17 ]),
    ...         ('s2', [ 52,   68 ]),
    ...         ('s3', [ 51,   17 ])]
    >>> merge = VoteMerger()
    >>> sdat = GenotripleStream(trips,genorepr=snp_marker).as_sdat(merge)
    >>> assert list(sdat) == rows and sdat.loci == loci
    >>> sdat = GenotripleStream(trips,loci=['l1','l2'],genorepr=snp_marker).as_sdat(merge)
    >>> assert list(sdat) == rows and sdat.loci == loci
    >>> sdat = GenotripleStream(trips,samples=['s1','s2','s3'],genorepr=snp_marker).as_sdat(merge)
    >>> assert list(sdat) == rows and sdat.loci == loci
    '''
    if mergefunc is None:
      mergefunc = UniqueMerger(missingrepr=self.genorepr.missing)

    loci,genos = build_genomatrixstream_from_triples(self, 'sdat', mergefunc=mergefunc,
                                          samples=self.samples, loci=self.loci, order=self.order,
                                          genorepr=self.genorepr)

    return GenomatrixStream(genos, 'sdat', loci=loci, unique=True, genorepr=self.genorepr)


class GenomatrixStream(GenotypeStream):
  '''
  A stream of genomatrix by sample or locus with optional metadata
  '''
  def __init__(self, genos, format, samples=None, loci=None, unique=False,
                     genorepr=snp_acgt, materialized=False, packed=False):
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
    @param     genorepr: object representing the input/output encoding and
                         internal representation of genotypes
    @type      genorepr: UnphasedMarkerRepresentation or similar object
    @param materialized: flag indicating if this stream is materialized and
                         allows iteration multiple times
    @type  materialized: bool
    @param       packed: flag indicating if genotypes are packed into a
                         compressed array format
    @type        packed: bool
    '''

    if format not in ('sdat','ldat'):
      raise ValueError("Invalid genomatrix format '%s'.  Must be either sdat or ldat" % format)

    if loci:
      loci = tuple(loci)

    if samples:
      samples = tuple(samples)

    if format=='ldat' and samples is None:
      raise ValueError('sample metadata are required for ldat streams')

    if format=='sdat' and loci is None:
      raise ValueError('locus metadata are required for ldat streams')


    materialized = isinstance(genos, (list,tuple))

    # If materialized, we can supplement our metadata
    if materialized:
      #FIXME: This restriction can be relaxed at some point
      assert isinstance(genos, (list,tuple)), 'Materialized streams currently support only lists and tuples'

      rowlabels = tuple(imap(itemgetter(0),genos))

      if format == 'sdat':
        if samples is not None and samples!=rowlabels:
          raise ValueError('invalid genomatrixstream metadata')
        samples = rowlabels

      elif format == 'ldat':
        if loci is not None and loci!=rowlabels:
          raise ValueError('invalid genomatrixstream metadata')
        loci    = rowlabels

    self.format       = format
    self.samples      = samples
    self.loci         = loci
    self.unique       = unique
    self.materialized = materialized
    self.genorepr     = genorepr
    self.packed       = packed

    self.__private_genos_do_not_touch = genos

  @staticmethod
  def from_streams(genos, format, mergefunc=None, genorepr=snp_acgt):
    '''
    Combine multiple genostreams into one genomatrix stream

    @param     genos: genostreams
    @type      genos: list
    @param    format: format of input genomatrix, either 'ldat' or 'sdat'
    @type     format: str
    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @param  genorepr: object representing the input/output encoding and
                      internal representation of genotypes
    @type   genorepr: UnphasedMarkerRepresentation or similar object
    @return         : combined genotriple stream
    @rtype          : sequence of sample, locus, and genotype

    >>> trip1 = [('s1','l1', 51),('s1','l2', 17),
    ...          ('s2','l2', 68),('s3','l1', 68)]
    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ 51,   52,   68 ]),
    ...         ('l2', [ 17,   68,   51 ])]
    >>> trip2 = [('s2','l1', 51),('s3','l2', 17)]
    >>> streams = [GenotripleStream(trip1),
    ...            GenomatrixStream(rows,'ldat',samples=samples),
    ...            GenotripleStream(trip2)]
    >>> merger=VoteMerger(snp_marker.missing)
    >>> combined = GenomatrixStream.from_streams(streams, 'ldat', genorepr=snp_marker, mergefunc=merger)
    >>> combined.columns
    ('s1', 's2', 's3')
    >>> for row in combined:
    ...   print row
    ('l1', [('G', 'G'), 0, ('T', 'T')])
    ('l2', [('A', 'A'), ('T', 'T'), 0])

    >>> streams = [GenomatrixStream(rows,'ldat',samples=samples),
    ...            GenomatrixStream(rows,'ldat',samples=samples)]
    >>> combined = GenomatrixStream.from_streams(streams, 'sdat', genorepr=snp_marker, mergefunc=merger)
    >>> combined.columns
    ('l1', 'l2')
    >>> for row in combined:
    ...   print row
    ('s1', [('G', 'G'), ('A', 'A')])
    ('s2', [('G', 'T'), ('T', 'T')])
    ('s3', [('T', 'T'), ('G', 'G')])
    '''
    if format not in ('sdat','ldat'):
      raise ValueError, "Invalid genomatrix format '%s'.  Must be either sdat or ldat" % format

    formats  = set(g.format for g in genos)
    informat = formats.pop() if len(formats)==1 else None
    headers  = [ g.columns for g in genos ] if informat in ('ldat','sdat') else None

    # Single input is trivial -- just transform to the desired genorepr and
    # mergefunc
    if len(genos) == 1:
      genos = genos[0].transformed(genorepr=genorepr,mergefunc=mergefunc)

    # Inputs are all matricies, without knowledge of uniqueness of rows.
    elif formats <= set(['ldat','sdat']):
      genos   = [ (g.as_ldat() if format=='ldat' else g.as_sdat()).transformed(genorepr=genorepr) for g in genos ]
      headers = [ g.columns for g in genos ]
      columns,genos = merge_genomatrixstream_list(headers, genos, mergefunc)

      if format == 'ldat':
        genos = GenomatrixStream(genos,format,samples=columns,unique=True,genorepr=genorepr)
      else:
        genos = GenomatrixStream(genos,format,loci=columns,unique=True,genorepr=genorepr)

    # Very general strategy of converting all inputs to genotriples, sorting
    # by the appopriate clustering, and transforming into ldat or sdat.
    # Unfortunately, very slow.
    else:
      order = 'locus' if format=='ldat' else 'sample'
      genos = GenotripleStream.from_streams(genos, order=order, genorepr=genorepr)

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

  def __iter__(self):
    '''
    Returns the embedded genomatrix stream and marks it as used (ie
    unavailable for further operations) if not already used.  Otherwise,
    raises a RuntimeError exception.

    @return: genomatrix stream
    @rtype : sequence of sample, locus, and genotype
    '''
    return iter(self.use_stream())

  def used_stream(self):
    '''
    Returns True if the genomatrix stream has been used as an iterable and is
    thus no longer available.  Otherwise, returns False.

    @return: availability of the genotriple stream
    @rtype : bool
    '''
    return self.__private_genos_do_not_touch is None

  def use_stream(self):
    '''
    Returns the embedded genomatrix stream and marks it as used (ie
    unavailable for further operations) if not already used.  Otherwise,
    raises a RuntimeError exception.

    @return: genomatrix stream
    @rtype : sequence of sample, locus, and genotype
    '''
    if self.used_stream():
      raise RuntimeError, 'Genomatrix stream already used'

    if not self.materialized:
      self.__private_genos_do_not_touch,genos = None,self.__private_genos_do_not_touch
    else:
      genos = self.__private_genos_do_not_touch

    return genos

  rows    = property(_get_rows,   _set_rows)
  columns = property(_get_columns,_set_columns)
  genos   = property(use_stream)

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

    return GenomatrixStream(list(genos), format=genos.format, samples=genos.samples, loci=genos.loci,
                                         unique=genos.unique, genorepr=genos.genorepr, materialized=True,
                                         packed=True)

  def recoded(self, genorepr):
    '''
    Returns a new genotriple stream with genotypes recoded to the new representation

    @param genorepr: object representing the input/output encoding and
                     internal representation of genotypes
    @type  genorepr: UnphasedMarkerRepresentation or similar object
    @return        : recoded genotriple stream
    @rtype         : GenotripleStream

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [ 51,   52,   68 ]),
    ...         ('l2', [ 17,   68,   51 ])]
    >>> ldat = GenomatrixStream(rows,'ldat',samples=samples)
    >>> recoded = ldat.recoded(snp_marker)
    >>> recoded.columns
    ('s1', 's2', 's3')
    >>> recoded = list(recoded)
    >>> for row in recoded:
    ...   print row
    ('l1', [('G', 'G'), ('G', 'T'), ('T', 'T')])
    ('l2', [('A', 'A'), ('T', 'T'), ('G', 'G')])
    >>> backagain = GenomatrixStream(recoded,'ldat',samples=samples,genorepr=snp_marker)
    >>> for triple in backagain.recoded(snp_acgt).as_genotriples():
    ...   print triple
    ('s1', 'l1', 51)
    ('s2', 'l1', 52)
    ('s3', 'l1', 68)
    ('s1', 'l2', 17)
    ('s2', 'l2', 68)
    ('s3', 'l2', 51)
    '''
    return self.transformed(genorepr=genorepr)

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
    >>> rows = [('l1', [ 51,   52,   68 ]),
    ...         ('l2', [ 17,   68,   51 ])]
    >>> genos = GenomatrixStream(rows,'ldat',samples=samples).transformed(include_loci=['l1'],exclude_samples=['s3'])
    >>> genos.samples
    ('s1', 's2')
    >>> for row in genos:
    ...   print row
    ('l1', [51, 52])
    '''
    packed = self.packed

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

    genos   = self.use_stream()
    rows    = self.rows
    columns = self.columns

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
    # 5) Adjust genorepr
    # 6) Order and repack

    # Apply row includes and excludes
    if rowtransform.exclude:
      columns,genos = filter_genomatrixstream_by_row(columns,genos,rowtransform.exclude,exclude=True)
      if rows is not None:
        rows = [ r for r in self.rows if r not in rowtransform.exclude ]
    if rowtransform.include is not None:
      columns,genos = filter_genomatrixstream_by_row(columns,genos,rowtransform.include)
      if rows is not None:
        rows = [ r for r in self.rows if r in rowtransform.include ]

    # Apply column includes and excludes
    if coltransform.exclude:
      columns,genos = filter_genomatrixstream_by_column(columns,genos,coltransform.exclude,exclude=True)
      packed = False
    if coltransform.include is not None:
      columns,genos = filter_genomatrixstream_by_column(columns,genos,coltransform.include)
      packed = False

    # Determine if resulting data will be unique (before renaming samples and loci)
    if self.format == 'ldat':
      samples = columns
      loci    = rows
    else:
      loci    = columns
      samples = rows

    unique = prove_unique_transform(transform=transform,loci=loci,samples=samples,unique=self.unique)

    # Apply renamings
    if rowtransform.rename:
      columns,genos = rename_genomatrixstream_row(columns,genos,rowtransform.rename)
      if rows is not None:
        rows = [ rowtransform.rename.get(r,r) for r in rows ]
    if coltransform.rename is not None:
      columns,genos = rename_genomatrixstream_column(columns,genos,coltransform.rename)

    # Filter rows and columns with all missing data
    if transform.filter_missing_genotypes:
      rows    = None
      columns = None
      packed  = False
      columns,genos = filter_genomatrixstream_missing(columns,genos)

    # Update genorepr
    genorepr = self.genorepr
    if transform.genorepr is not None:
      genorepr = transform.genorepr
      if isinstance(genorepr, basestring):
        genorepr = get_genorepr(genorepr)

    if genorepr is not self.genorepr or transform.rename_alleles:
      columns,genos = recode_genomatrixstream(columns,genos,self.genorepr,genorepr,self.format,transform.rename_alleles)
      packed = True

    if self.format == 'sdat':
      samples = rows
      loci    = columns
    else:
      samples = columns
      loci    = rows

    genos = GenomatrixStream(genos, self.format, samples=samples, loci=loci, unique=unique,
                             genorepr=genorepr, packed=packed)

    # Ordering by [] and None are distinct cases: the first will order in lexicographical order
    # as a side-effect, while the latter should not alter order.
    if transform.samples.order is not None or transform.loci.order is not None:
      genos = genos.sorted(locusorder=transform.loci.order, sampleorder=transform.samples.order)

    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if transform.repack and not genos.packed:
      columns,packed = pack_genomatrixstream(genos.columns, genos,genos.genorepr)
      genos = GenomatrixStream(packed, format=genos.format, samples=genos.samples, loci=genos.loci,
                               unique=genos.unique, genorepr=genos.genorepr, materialized=False, packed=True)

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
    >>> rows = [('s1', [ 51,   17 ]),
    ...         ('s2', [ 52,   68 ]),
    ...         ('s3', [ 51,   17 ])]
    >>> genos = GenomatrixStream(rows,'sdat',loci=loci,genorepr=snp_marker)
    >>> genos = genos.sorted(locusorder=['l2'],sampleorder=['s2','s1','s3'])
    >>> genos.loci
    ('l2', 'l1')
    >>> for row in genos:
    ...   print row
    ('s2', [68, 52])
    ('s1', [17, 51])
    ('s3', [17, 51])
    >>> genos = GenomatrixStream(rows,'sdat',loci=loci,genorepr=snp_marker).transposed()
    >>> genos = genos.sorted(locusorder=['l2'],sampleorder=['s2','s1','s3'])
    >>> genos.samples
    ('s2', 's1', 's3')
    >>> for row in genos:
    ...   print row
    ('l2', [68, 17, 17])
    ('l1', [52, 51, 51])
    '''
    if self.format == 'sdat':
      columnorder = locusorder
      roworder    = sampleorder
    else:
      columnorder = sampleorder
      roworder    = locusorder

    genos = self
    if roworder is not None:
      genos = self.transformed(repack=True)

    packed  = genos.packed
    genos   = genos.use_stream()
    columns = self.columns

    if columnorder is not None:
      columns,genos = reorder_genomatrixstream_columns(columns,genos,columnorder)
      packed = False

    # FIXME: Broadcast rows from reorder_genomatrix_rows since this requires
    # materialization
    if roworder is not None:
      columns,genos = reorder_genomatrixstream_rows(columns,genos,roworder)

    if self.format == 'sdat':
      loci    = columns
      samples = None
    else:
      loci    = None
      samples = columns

    return GenomatrixStream(genos, self.format, loci=loci, samples=samples, unique=self.unique,
                                   genorepr=self.genorepr, packed=packed)

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
      # Pack, since merge will materialize
      columns,genos = merge_genomatrixstream(self.columns, self.transformed(repack=True), mergefunc)
    else:
      columns,genos = merge_genomatrixstream_columns(self.columns, self.genos, mergefunc)

    if self.format == 'sdat':
      genos = GenomatrixStream(genos, self.format, loci=columns, unique=True, genorepr=self.genorepr, packed=False)
    else:
      genos = GenomatrixStream(genos, self.format, samples=columns, unique=True, genorepr=self.genorepr, packed=False)

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
    >>> rows = [('s1', [ 51,   17 ]),
    ...         ('s2', [ 52,   68 ]),
    ...         ('s3', [ 51,   17 ])]
    >>> genos = GenomatrixStream(rows,'sdat',loci=loci,genorepr=snp_marker).transposed()
    >>> print genos.format
    ldat
    >>> genos.rows
    ('l1', 'l2')
    >>> genos.columns
    ('s1', 's2', 's3')
    >>> for row in genos:
    ...   print row
    ('l1', [51, 52, 51])
    ('l2', [17, 68, 17])
    '''
    rows,genos = transpose_generator(self.columns, list(self), missing=self.genorepr.missing)

    assert self.rows is None or all(r1==r2 for r1,r2 in izip(self.rows,rows))

    self.rows = tuple(rows)

    if self.format == 'ldat':
      format = 'sdat'
    else:
      format = 'ldat'

    return GenomatrixStream(genos, format, samples=self.samples, loci=self.loci, unique=self.unique,
                                           genorepr=self.genorepr, packed=False)

  def as_genotriples(self):
    '''
    Return the current genomatrix data as a GenotripleStream.
    Ordered triple streams can be transformed without full materialization,
    though unordered streams do require materialization.  A merge function is
    required for non-unique streams.

    @return         : genotriples converted into a genomatrix stream
    @rtype          : GenomatrixStream

    >>> trips = [('s1','l1', 51),('s1','l2', 17),
    ...          ('s2','l1', 52),('s2','l2', 68),
    ...          ('s3','l1', 51),('s3','l2', 17)]
    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ 51,   17 ]),
    ...         ('s2', [ 52,   68 ]),
    ...         ('s3', [ 51,   17 ])]
    >>> triples = GenomatrixStream(rows,'sdat',loci=loci,genorepr=snp_marker).as_genotriples()
    >>> for row in triples:
    ...   print row
    ('s1', 'l1', 51)
    ('s1', 'l2', 17)
    ('s2', 'l1', 52)
    ('s2', 'l2', 68)
    ('s3', 'l1', 51)
    ('s3', 'l2', 17)
    '''

    if self.format == 'ldat':
      triples = build_genotriples_by_locus(self.columns,self)
      order   = 'locus'
    else:
      triples = build_genotriples_by_sample(self.columns,self)
      order   = 'sample'

    # Unset the result order if any of the required ordering constraints
    # cannot be verified
    if (self.samples is None or self.loci is None
                             or sorted(self.samples) != self.samples
                             or sorted(self.loci)    != self.loci):
      order = None

    return GenotripleStream(triples, unique=self.unique, order=order, genorepr=self.genorepr)

  def as_ldat(self, mergefunc=None):
    '''
    Return a genomatrix stream in ldat format.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @return         : ldat genomatrix stream
    @rtype          : GenomatrixStream

    >>> loci = ['l1', 'l2']
    >>> rows = [('s1', [ 51,   17 ]),
    ...         ('s2', [ 52,   68 ]),
    ...         ('s3', [ 51,   17 ])]
    >>> genos = GenomatrixStream(rows,'sdat',loci=loci,genorepr=snp_marker).as_ldat()
    >>> genos.samples
    ('s1', 's2', 's3')
    >>> for row in genos:
    ...   print row
    ('l1', [51, 52, 51])
    ('l2', [17, 68, 17])
    '''
    genos = self
    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if genos.format == 'ldat':
      return genos
    else:
      return genos.transposed()

  def as_sdat(self, mergefunc=None):
    '''
    Return a genomatrix stream in sdat format.

    @param mergefunc: function to merge multiple genotypes into a consensus genotype
    @type  mergefunc: callable
    @return         : sdat genomatrix stream
    @rtype          : GenomatrixStream

    >>> samples = ['s1', 's2', 's3']
    >>> rows = [('l1', [51, 52, 51]),
    ...         ('l2', [17, 68, 17])]
    >>> genos = GenomatrixStream(rows,'ldat',samples=samples,genorepr=snp_marker).as_sdat()
    >>> genos.loci
    ('l1', 'l2')
    >>> for row in genos:
    ...   print row
    ('s1', [51, 17])
    ('s2', [52, 68])
    ('s3', [51, 17])
    '''
    genos = self
    if mergefunc is not None:
      genos = genos.merged(mergefunc)

    if genos.format == 'sdat':
      return genos
    else:
      return genos.transposed()

#######################################################################################

list_type = (NoneType,set,dict,list,tuple)
map_type  = (NoneType,dict)


class GenoTransform(object):
  '''
  Create a GenoTransform object to specify various transformation on the genodata.
  Supported operations: include/exclude/rename samples or loci; optional filter to remove missing genotypes
  '''
  def __init__(self, include_samples, exclude_samples, rename_samples, order_samples,
                     include_loci,    exclude_loci,    rename_loci,    order_loci,
                     rename_alleles=None, filter_missing=False, genorepr=None, repack=False):
    '''
    Create a new GenoTransform object with supplied metadata,
    which are used to specify all the operations of transforming the genostream
    and thus must be accurate or else incorrect results are virtually guaranteed.
    When in doubt, do not specify them, as each algorithm can compensate.

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
    @param  rename_alleles: rename alleles from any loci and allele name to new allele name
    @type   rename_alleles: dict from str -> old_allele str -> new_allele str
    @param  filter_missing: filter missing genotypes from the stream
    @type   filter_missing: bool
    @param        genorepr: new genotype representation to recode
    @type         genorepr: UnphasedMarkerRepresentation or similar object
    @param          repack: trigger repacking of genotypes to ensure that the most compact storage
                            method is used
    @type         genorepr: bool
    @return               : transformed genotriple stream
    @rtype                : GenotripleStream
    '''
    self.samples = GenoSubTransform(include_samples, exclude_samples, rename_samples, order_samples)
    self.loci    = GenoSubTransform(include_loci,    exclude_loci,    rename_loci,    order_loci)

    if not isinstance(rename_alleles, map_type):
      rename_alleles = load_rename_alleles_file(rename_alleles)

    self.filter_missing_genotypes = filter_missing
    self.rename_alleles           = rename_alleles
    self.genorepr                 = genorepr
    self.repack                   = repack

  @staticmethod
  def from_options(options):
    '''
    Create a new GenoTransform object from command line option list

    @return: transformed genotriple stream
    @rtype : GenotripleStream
    '''
    return GenoTransform(options.includesamples, options.excludesamples, options.renamesamples, options.ordersamples,
                         options.includeloci,    options.excludeloci,    options.renameloci,    options.orderloci,
                         rename_alleles=options.renamealleles, filter_missing=options.filtermissing,
                         genorepr=options.outgenorepr)

  @staticmethod
  def from_kwargs(**kwargs):
    '''
    Create a new GenoTransform object from key word arguments

    @return: transformed genotriple stream
    @rtype : GenotripleStream
    '''
    transform = GenoTransform(include_samples=kwargs.pop('include_samples',None),
                              exclude_samples=kwargs.pop('exclude_samples',None),
                               rename_samples=kwargs.pop('rename_samples', None),
                                order_samples=kwargs.pop('order_samples',  None),
                                 include_loci=kwargs.pop('include_loci',   None),
                                 exclude_loci=kwargs.pop('exclude_loci',   None),
                                  rename_loci=kwargs.pop('rename_loci',    None),
                                   order_loci=kwargs.pop('order_loci',     None),
                               filter_missing=kwargs.pop('filter_missing', False),
                                     genorepr=kwargs.pop('genorepr',       None),
                                       repack=kwargs.pop('repack',         False),
                               rename_alleles=kwargs.pop('rename_alleles', None))

    if kwargs:
      raise TypeError, "'%s' is an invalid keyword argument for this function" % kwargs.popitem()[0]

    return transform


class GenoSubTransform(object):
  '''
  A GenoSubTransform object with metadata related to samples or loci transformation
  '''
  def __init__(self, include, exclude, rename, order):
    '''
    Create a new GenoSubTransform object

    @param include: filter samples/loci such that they must appear in the set (optional)
    @type  include: set
    @param exclude: filter samples/loci such that they must not appear in the set (optional)
    @type  exclude: set
    @param  rename: rename any samples/loci that appear in the supplied dictionary to the
                            associated value (optional)
    @type   rename: dict from str -> str
    '''
    if not isinstance(include, list_type):
      include = set(load_list(include))

    if not isinstance(exclude, list_type):
      exclude = set(load_list(exclude))

    if not isinstance(rename, map_type):
      rename  = load_map(rename)

    if not isinstance(order, list_type):
      order = load_list(order)

    self.include = include
    self.exclude = exclude
    self.rename  = rename
    self.order   = order


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

  >>> triples = [('s3','l1', 51),('s3','l2', 17),
  ...            ('s2','l2', 52),('s1','l1', 68),
  ...            ('s1','l1', 51),('s2','l2', 17)]
  >>> samples,loci,striples = sort_genotriples(triples,order='sample',sampleorder=['s1','s2','s3'],
  ...                                                                  locusorder=['l1','l2','l3'])
  >>> samples
  ['s1', 's2', 's3']
  >>> loci
  ['l1', 'l2']
  >>> striples
  [('s1', 'l1', 51), ('s1', 'l1', 68), ('s2', 'l2', 17), ('s2', 'l2', 52), ('s3', 'l1', 51), ('s3', 'l2', 17)]
  >>> samples,loci,striples = sort_genotriples(triples,order='sample',locusorder=['l3','l2','l1'])
  >>> samples
  ['s1', 's2', 's3']
  >>> loci
  ['l1', 'l2']
  >>> striples
  [('s1', 'l1', 51), ('s1', 'l1', 68), ('s2', 'l2', 17), ('s2', 'l2', 52), ('s3', 'l2', 17), ('s3', 'l1', 51)]
  >>> samples,loci,striples = sort_genotriples(triples,order='sample',sampleorder=['s3','s2','s1'])
  >>> samples
  ['s3', 's2', 's1']
  >>> loci
  ['l1', 'l2']
  >>> striples
  [('s3', 'l1', 51), ('s3', 'l2', 17), ('s2', 'l2', 17), ('s2', 'l2', 52), ('s1', 'l1', 51), ('s1', 'l1', 68)]
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
  triples = sorted(triples,key=keyfunc)

  # Extract sample and locus sets since this can be done quickly with materialized lists
  samples = list(unique(imap(itemgetter(0), triples)))
  loci    = list(unique(imap(itemgetter(1), triples)))

  return samples,loci,triples


def combine_unsorted_genotriple_list(triplelist):
  '''
  Combine multiple genotriples into one unsorted triple stream

  @param triplelist: list of genotriples
  @type  triplelist: list
  @return          : combined genotriple stream
  @rtype           : sequence of sample, locus, and genotype

  >>> trip1 = [('s3','l1', 51),('s3','l2', 17),
  ...          ('s2','l2', 52),('s1','l1', 68)]
  >>> trip2 = [('s1','l1', 51),('s2','l2', 17)]
  >>> triplelist = [GenotripleStream(trip1),GenotripleStream(trip2)]
  >>> combined = combine_unsorted_genotriple_list(triplelist)
  >>> for row in combined:
  ...   print row
  ('s3', 'l1', 51)
  ('s3', 'l2', 17)
  ('s2', 'l2', 52)
  ('s1', 'l1', 68)
  ('s1', 'l1', 51)
  ('s2', 'l2', 17)
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
  reprs   = [ triples.genorepr for triples in triplelist ]

  # Check that all representations align
  reprs = set(reprs)

  if len(reprs) != 1:
    raise ValueError, 'Cannot merge triplestreams in disparate representations'

  genorepr = reprs.pop()

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
  return GenotripleStream(triples, samples=samples, loci=loci, order=None, unique=unique, genorepr=genorepr)


def combine_sorted_genotriple_list(triplelist):
  '''
  Combine multiple sorted genotriples into one sorted triple stream

  @param triplelist: list of genotriples
  @type  triplelist: list
  @return          : combined genotriple stream
  @rtype           : sequence of sample, locus, and genotype

  >>> trip1 = [('s1','l1', 51),('s1','l2', 17),
  ...          ('s2','l2', 52),('s3','l1', 68)]
  >>> trip2 = [('s2','l1', 51),('s3','l2', 17)]
  >>> triplelist = [GenotripleStream(trip1,order='sample'),GenotripleStream(trip2,order='sample')]
  >>> combined = combine_sorted_genotriple_list(triplelist)
  >>> for row in combined:
  ...   print row
  ('s1', 'l1', 51)
  ('s1', 'l2', 17)
  ('s2', 'l1', 51)
  ('s2', 'l2', 52)
  ('s3', 'l1', 68)
  ('s3', 'l2', 17)
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
  reprs   = [ triples.genorepr for triples in triplelist ]

  # Check that all representations align
  reprs = set(reprs)

  if len(reprs) != 1:
    raise ValueError, 'Cannot merge triplestreams in disparate representations'

  genorepr = reprs.pop()

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
  return GenotripleStream(triples, samples=samples, loci=loci, order=order, unique=unique, genorepr=genorepr)


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

  >>> t = [('l1','s1','AA'),('l1','s1',   0),('l1','s2','AB'),('l2','s1','AA'),
  ...      ('l2','s1','AA'),('l3','s1','BB'),('l3','s1','BB'),('l3','s1','AB')]
  >>> list(merge_sorted_genotriples(t,VoteMerger()))
  [('l1', 's1', 'AA'), ('l1', 's2', 'AB'), ('l2', 's1', 'AA'), ('l3', 's1', 0)]
  '''
  for (sample,locus),trips in groupby(triples, itemgetter(0,1)):
    yield sample,locus,mergefunc(sample,locus,map(itemgetter(2),trips))


def merge_genomatrixstream_columns(columns, genos, mergefunc):
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
  >>> rows = [('s1',[17,   0,  51]),
  ...         ('s2',[ 0,   0,   0]),
  ...         ('s3',[17,   0,   0]),
  ...         ('s4',[52,   0,  68])]
  >>> loci,rows = merge_genomatrixstream_columns(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in rows:
  ...   print row
  ('s1', [17, 0, 51])
  ('s2', [0, 0, 0])
  ('s3', [17, 0, 0])
  ('s4', [52, 0, 68])

  >>> loci = ('l1','l2','l1')
  >>> rows = [('s1',[17,   0,  51]),
  ...         ('s2',[ 0,  11,   0]),
  ...         ('s3',[17,  17,   0]),
  ...         ('s4',[52,   0,  52])]
  >>> loci,rows = merge_genomatrixstream_columns(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2')
  >>> for row in rows:
  ...   print row
  ('s1', [0, 0])
  ('s2', [0, 11])
  ('s3', [17, 17])
  ('s4', [52, 0])
  '''
  assert mergefunc is not None

  merge_indices = {}
  new_columns   = []

  for i,column in enumerate(columns):
    if column not in merge_indices:
      new_columns.append(column)
      merge_indices[column] = [i]
    else:
      merge_indices[column].append(i)

  new_columns = tuple(new_columns)

  # Trivial path: no merge is needed
  if new_columns == columns:
    return columns,genos

  # Non-trivial path: one or more columns must be merged
  def _merger():
    rows_seen = set()
    for row_label,row in genos:
      if row_label in rows_seen:
        raise ValueError('row labels required to be unique')

      rows_seen.add(row_label)

      merge_columns = ([row[i] for i in merge_indices[col_label]] for col_label in new_columns)
      new_row = list(imap(mergefunc,repeat(row_label),columns,merge_columns))

      yield row_label,new_row

  return new_columns,_merger()


def merge_genomatrixstream_rows(columns, genos, mergefunc):
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
  >>> rows = [('s1',[17,   0,  51]),
  ...         ('s2',[ 0,   0,   0]),
  ...         ('s3',[17,   0,   0]),
  ...         ('s4',[52,   0,  68])]
  >>> loci,rows = merge_genomatrixstream_rows(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in rows:
  ...   print row
  ('s1', [17, 0, 51])
  ('s2', [0, 0, 0])
  ('s3', [17, 0, 0])
  ('s4', [52, 0, 68])

  >>> rows = [('s1',[17,   0,  52]),
  ...         ('s2',[ 0,  11,   0]),
  ...         ('s1',[17,  17,   0]),
  ...         ('s1',[52,   0,  52])]
  >>> loci,rows = merge_genomatrixstream_rows(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in rows:
  ...   print row
  ('s1', [0, 17, 52])
  ('s2', [0, 11, 0])
  '''
  assert mergefunc is not None

  non_unique = any( n>1 for n in tally(columns).itervalues() )

  if non_unique:
    raise ValueError('column labels required to be unique')

  def _merger():
    merge_rows = defaultdict(list)
    new_rows   = []

    for row_label,row in genos:
      if row_label not in merge_rows:
        new_rows.append(row_label)
      merge_rows[row_label].append(row)

    for row_label in new_rows:
      merge_columns = izip(*merge_rows.pop(row_label))
      new_row = list(imap(mergefunc,repeat(row_label),columns,merge_columns))
      yield row_label,new_row

  return columns,_merger()


def merge_genomatrixstream(columns, genos, mergefunc):
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
  >>> rows = [('s1',[17,   0,  51]),
  ...         ('s2',[ 0,   0,   0]),
  ...         ('s3',[17,   0,   0]),
  ...         ('s4',[52,   0,  68])]
  >>> loci,genos = merge_genomatrixstream(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in rows:
  ...   print row
  ('s1', [17, 0, 51])
  ('s2', [0, 0, 0])
  ('s3', [17, 0, 0])
  ('s4', [52, 0, 68])

  >>> rows = [('s1',[17,   0,  52]),
  ...         ('s2',[ 0,  11,   0]),
  ...         ('s1',[17,  17,   0]),
  ...         ('s1',[52,   0,  52])]
  >>> loci,rows = merge_genomatrixstream(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in rows:
  ...   print row
  ('s1', [0, 17, 52])
  ('s2', [0, 11, 0])

  >>> loci = ('l1','l2','l1')
  >>> rows = [('s1',[ 0,   0,  52]),
  ...         ('s2',[ 0,  11,   7]),
  ...         ('s1',[17,  17,   0]),
  ...         ('s1',[52,   0,  52])]
  >>> loci,rows = merge_genomatrixstream(loci,rows,VoteMerger())
  >>> loci
  ('l1', 'l2')
  >>> for row in rows:
  ...   print row
  ('s1', [0, 17])
  ('s2', [7, 11])
  '''
  assert mergefunc is not None

  merge_indices = defaultdict(list)

  new_columns = []
  for i,column in enumerate(columns):
    if column not in merge_indices:
      new_columns.append(column)
    merge_indices[column].append(i)

  new_columns = tuple(new_columns)

  def _merger():
    merge_rows = defaultdict(list)
    new_rows   = []

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
        new_row = list(imap(mergefunc,repeat(row_label),columns,merge_columns))

        # Free used rows
        del merge_rows[row_label]

        # Yield new row
        yield row_label,new_row

  return new_columns,_merger()


def merge_genomatrixstream_list(columns, genos, mergefunc):
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
  >>> rows1 = [('l1',[('G', 'G'), ('G', 'T'), ('T', 'T')]),
  ...          ('l2',[('A', 'A'), ('T', 'T'), ('G', 'G')])]
  >>> samples2 =        ('s1',       's3',       's4')
  >>> rows2 = [('l1',[   None,    ('A', 'T'), ('C', 'C')]),
  ...          ('l3',[('A', 'A'),  None,      ('A', 'T')])]
  >>> columns = [samples1,samples2,samples1]
  >>> genos   = [rows1,rows2,rows1]
  >>> samples,rows = merge_genomatrixstream_list(columns,genos,VoteMerger(missingrepr=None))
  >>> samples
  ('s1', 's2', 's3', 's4')
  >>> for row in rows:
  ...   print row
  ('l1', [('G', 'G'), ('G', 'T'), None, ('C', 'C')])
  ('l2', [('A', 'A'), ('T', 'T'), ('G', 'G'), None])
  ('l3', [('A', 'A'), None, None, ('A', 'T')])

  Test fast-path for homogeneous schema:

  >>> columns =        [('s1',       's2',       's3')]*3
  >>> rows1 = [('l1',[('G', 'G'), ('G', 'T'), ('T', 'T')]),
  ...          ('l2',[('A', 'A'), ('T', 'T'), ('G', 'G')])]
  >>> rows2 = [('l1',[   None,    ('A', 'T'), ('T', 'T')]),
  ...          ('l3',[('A', 'A'),  None,      ('A', 'T')])]
  >>> genos = [rows1,rows2,rows1]
  >>> samples,rows = merge_genomatrixstream_list(columns,genos,VoteMerger(missingrepr=None))
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', [('G', 'G'), None, ('T', 'T')])
  ('l2', [('A', 'A'), ('T', 'T'), ('G', 'G')])
  ('l3', [('A', 'A'), None, ('A', 'T')])
  '''
  assert mergefunc is not None

  if not genos:
    raise ValueError('Invalid empty genotype list')

  # Fastpath for identical schema: Concatenate all genotype rows and merge
  # using the single-matrix merge function.  Also optimizes the single input
  # case.
  if all(columns[0]==c for c in columns):
    # Pass-through from merge_genomatrix
    return merge_genomatrixstream(columns[0], chain(*genos), mergefunc)

  new_rows      = {}
  new_columns   = {}
  merge_columns = defaultdict(list)
  merge_rows    = defaultdict(lambda: defaultdict(list))

  for i,(g,c) in enumerate(izip(genos,columns)):
    # Collect columns and form mappings from old schemas to new
    for j,column in enumerate(c):
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

  return new_columns,_merger()


#######################################################################################


def prove_bijective_mapping(items,transform):
  '''
  Construct the minimal sample reverse map by removing excluded items
  to verify that no two map to the same identifier.

  @param     items: sequence of samples/loci if known, otherwise None
  @type      items: sequence of str or None
  @param transform: transformation object
  @type  transform: GenoTransform object
  @return         : uniqueness of the mapping
  @rtype          : bool

  >>> samples = ['s1', 'ns1', 's2','s3']
  >>> loci = ['l1','l2','l3','l4']
  >>> rename_samples = {'s1':'ns1','s2':'ns2','s3':'ns3'}
  >>> include_samples = ['s1', 'ns1', 's2']
  >>> rename_loci = {'l1':'nl1','l2':'nl2','l3':'nl3','l4':'nl4'}
  >>> include_loci = ['l1','l2','l3']
  >>> transform = GenoTransform(include_samples, None, rename_samples,None,include_loci,None,rename_loci,None)
  >>> prove_bijective_mapping(samples,transform.samples)
  False
  >>> prove_bijective_mapping(loci,transform.loci)
  True
  '''
  # Construct the minimal sample reverse map by removing excluded items
  # and verify that no two items map to the same identifier.
  if not transform.rename:
    return True

  # Cannot prove uniqueness when a renaming without knowing the universe of
  # possible items
  if items is None and transform.include is None:
    return False
  elif items is not None and transform.include is not None:
    items = set(items) & set(transform.include)
  elif transform.include is not None:
    items = set(transform.include)

  # Construct the minimal sample reverse map by removing excluded items
  # to verify that no two map to the same identifier.
  reverse_map = {}
  for item in items:
    renamed = transform.rename.get(item,item)
    reverse_map.setdefault(renamed,set()).add(item)

  # Mapping is unique if and only if all reverse map values are unique
  return all( len(v)<=1 for v in reverse_map.itervalues() )


def prove_unique_transform(transform=None,samples=None,loci=None,unique=False):
  '''
  Prove uniqueness of transformation operations

  @param transform: transformation object (optional)
  @type  transform: GenoTransform object
  @param   samples: optional set of samples refered to by the triples
  @type    samples: sequence, set, or None
  @param      loci: optional set of samples refered to by the triples
  @type       loci: sequence, set, or None
  @param    unique: flag indicating if repeated elements do not exist within the stream
  @type     unique: bool
  @return         : uniqueness of resulting triples
  @rtype          : bool
  '''

  # If the data aren't unique coming in, then we must assume they will not
  # be after a transformation
  if not unique:
    return False

  if transform is None:
    return True

  # Construct the minimal sample reverse map by removing excluded samples
  # and loci to verify that no two samples map to the same identifier.
  return   (prove_bijective_mapping(samples, transform.samples) and
            prove_bijective_mapping(loci,    transform.loci))


#######################################################################################


#FIXME: Extend xtab to return rows and columns since it has trivial knowledge of them
#FIXME: Convert to build_genomatixstream_from_genotriplestream
def build_genomatrixstream_from_triples(genos, format, mergefunc, samples=None, loci=None, order=None, genorepr=None):
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

  >>> rows = [('s1','l1', 51),('s1','l2', 17),
  ...         ('s2','l1', 52),('s2','l2', 68),
  ...         ('s3','l1', 51),('s3','l2', 17)]
  >>> merge = UniqueMerger()
  >>> loci,new_rows = build_genomatrixstream_from_triples(rows,'sdat',merge,loci=['l1','l2'],order='sample')
  >>> loci
  ('l1', 'l2')
  >>> for row in new_rows:
  ...   print row
  ('s1', [51, 17])
  ('s2', [52, 68])
  ('s3', [51, 17])

  >>> rows = [('s1','l1', 51),('s1','l1', 17),
  ...         ('s2','l1', 52),('s2','l1', 68),
  ...         ('s3','l1', 51),('s3','l1', 17)]
  >>> loci,new_rows = build_genomatrixstream_from_triples(rows,'sdat',merge,loci=['l1','l2'],order='sample')
  >>> for row in new_rows:
  ...   print row
  Traceback (most recent call last):
       ...
  NonUniqueGenotypeError: Non-unique genotype found

  >>> rows = [('s1','l1', 51),('s1','l2', 17),
  ...         ('s2','l1', 52),('s2','l2', 68),
  ...         ('s3','l1', 51),('s3','l2', 17)]
  >>> merge = VoteMerger()
  >>> loci,new_rows = build_genomatrixstream_from_triples(rows,'sdat',merge,loci=['l1','l2'],order='sample')
  >>> loci
  ('l1', 'l2')
  >>> for row in new_rows:
  ...   print row
  ('s1', [51, 17])
  ('s2', [52, 68])
  ('s3', [51, 17])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merge = VoteMerger()
  >>> loci,new_rows = build_genomatrixstream_from_triples(rows,'sdat',merge,samples=['s1','s2','s3'],order='sample')
  >>> loci
  ('l1', 'l2')
  >>> for row in new_rows:
  ...   print row
  ('s1', [51, 17])
  ('s2', [52, 68])
  ('s3', [51, 17])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]

  >>> merge = VoteMerger()
  >>> samples,new_rows = build_genomatrixstream_from_triples(rows,'ldat',merge,samples=['s1','s2','s3'],order='sample')
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', [51, 52, 51])
  ('l2', [17, 68, 17])
  >>> sorted(merge.samplestats.iteritems())
  [('s1', [2, 0, 0, 0, 0]), ('s2', [2, 0, 0, 0, 0]), ('s3', [2, 0, 0, 0, 0])]
  >>> sorted(merge.locusstats.iteritems())
  [('l1', [3, 0, 0, 0, 0]), ('l2', [3, 0, 0, 0, 0])]
  '''
  columns = None
  if format == 'sdat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(0),itemgetter(1),itemgetter(2)
    if not isinstance(loci, (NoneType,list)):
      columns = sorted(loci)
    else:
      columns = loci
  elif format == 'ldat':
    rowkeyfunc,colkeyfunc,valuefunc = itemgetter(1),itemgetter(0),itemgetter(2)
    mergefunc = mergefunc_transpose_adapter(mergefunc)
    if not isinstance(samples, (NoneType,list)):
      columns = sorted(samples)
    else:
      columns = samples
  else:
    raise NotImplementedError,'triple to %s format conversion is not supported' % format

  if format == 'sdat' and order != 'sample':
    order = False
  elif format == 'ldat' and order != 'locus':
    order = False

  if columns is None or not order:
    columns,rows,data = xtab(genos, rowkeyfunc, colkeyfunc, valuefunc, mergefunc)
    genos = izip(rows,data)
  else:
    columns,genos = rowsby(genos, columns, rowkeyfunc, colkeyfunc, valuefunc, mergefunc)

  columns = tuple(columns)
  return columns,genos


#######################################################################################


def pack_genomatrixstream(columns, genos, genorepr):
  '''
  Transform a genomatrix into an internal packed representation

  @param    columns: matrix column names
  @type     columns: sequence of strs
  @param      genos: genomatrix
  @type       genos: genomatrix generator
  @param   genorepr: object representing the input/output encoding and
                     internal representation of genotypes
  @type    genorepr: UnphasedMarkerRepresentation or similar object
  @return          : genomatrix with a packed internal format
  @rtype           : genomatrix generator

  >>> samples = ('s1','s2','s3')
  >>> genos = [('l1',[17,0,51]),
  ...          ('l2',[0,0,0]),
  ...          ('l3',[17,0,0]),
  ...          ('l4',[52,0,68])]
  >>> samples,rows = pack_genomatrixstream(samples,genos,genorepr=snp_acgt)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', array('B', [17, 0, 51]))
  ('l2', array('B', [0, 0, 0]))
  ('l3', array('B', [17, 0, 0]))
  ('l4', array('B', [52, 0, 68]))
  '''
  def _pack():
    for label,geno in genos:
      yield label,genorepr.pack_reps(geno)
  return columns,_pack()


def recode_genomatrixstream(columns, genos, old_genorepr, new_repr, format=None, rename_alleles=None):
  '''
  Returns a new genomatrix with the genotypes recoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix
  @type         genos: genomatrix generator
  @param old_genorepr: internal representation of genotypes to be transformed from
  @type  old_genorepr: Unph[AasedMarkerRepresentation or similar object
  @param     new_repr: internal representation of genotypes to be transformed to
  @type      new_repr: UnphasedMarkerRepresentation or similar object
  @return            : genomatrix with a packed internal format
  @rtype             : genomatrix generator

  >>> samples = ('s1','s2','s3')
  >>> genos = [('l1',[17,0,51]),
  ...          ('l2',[0,0,0]),
  ...          ('l3',[17,0,0]),
  ...          ('l4',[52,0,68])]
  >>> samples,rows = recode_genomatrixstream(samples,genos,snp_acgt,snp_marker)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', [('A', 'A'), None, ('G', 'G')])
  ('l2', [None, None, None])
  ('l3', [('A', 'A'), None, None])
  ('l4', [('G', 'T'), None, ('T', 'T')])

  >>> complement = {'A':'T','T':'A','C':'G','G':'C',None:None}
  >>> rename = {'l1':complement, 'l4':complement}
  >>> samples,rows = recode_genomatrixstream(samples,genos,snp_acgt,snp_marker,'ldat',rename)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', [('T', 'T'), None, ('C', 'C')])
  ('l2', [None, None, None])
  ('l3', [('A', 'A'), None, None])
  ('l4', [('C', 'A'), None, ('A', 'A')])

  >>> loci = ('l1','l2','l3','l4')
  >>> genos = [('s1',[17,0,17,52]),
  ...          ('s2',[0,0,0,0]),
  ...          ('s3',[51,0,0,68])]
  >>> loci,rows = recode_genomatrixstream(loci,genos,snp_acgt,snp_marker,'sdat',rename)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in rows:
  ...   print row
  ('s1', [('T', 'T'), None, ('A', 'A'), ('C', 'A')])
  ('s2', [None, None, None, None])
  ('s3', [('C', 'C'), None, None, ('A', 'A')])
  '''
  def _recode():
    if not rename_alleles:
      for label,row in genos:
        yield label,new_repr.pack_genos(old_genorepr.genos_from_reps(row))

    elif format=='ldat':
      for label,row in genos:
        row = old_genorepr.genos_from_reps(row)
        if label in rename_alleles:
          r   = rename_alleles[label]
          row = ( ((r[g[0]],r[g[1]]) if g else g) for g in row )
        yield label,new_repr.pack_genos(row)

    elif format=='sdat':
      remaps = [ rename_alleles.get(h) for h in columns ]
      for label,row in genos:
        row = old_genorepr.genos_from_reps(row)
        row = [ ((r[g[0]],r[g[1]]) if g and r else g) for g,r in izip(row,remaps) ]
        yield label,new_repr.pack_genos(row)
    else:
      raise ValueError('Matrix format must be specified when renaming alleles')

  return columns,_recode()


def recode_genotriples(triples, old_genorepr, new_genorepr, rename_alleles=None):
  '''
  Returns a new genotriples with the genotypes recoded to the a new internal representation

  @param      triples: genomatriple
  @type       triples: genomatriple generator
  @param old_genorepr: internal representation of genotypes to be transformed from
  @type  old_genorepr: UnphasedMarkerRepresentation or similar object
  @param new_genorepr: internal representation of genotypes to be transformed to
  @type  new_genorepr: UnphasedMarkerRepresentation or similar object
  @return            : genotriple with the new internal format
  @rtype             : genotriple generator

  >>> triples = [('s3','l1', 51),('s3','l2', 17),
  ...            ('s2','l2', 52),('s1','l1', 68),
  ...            ('s1','l1', 51),('s2','l2', 17)]
  >>> for row in recode_genotriples(triples,snp_acgt,snp_marker):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l2', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  >>> rename = {'l1':{'A':'T','T':'A','C':'G','G':'C',None:None}}
  >>> for row in recode_genotriples(triples,snp_acgt,snp_marker,rename):
  ...   print row
  ('s3', 'l1', ('C', 'C'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l2', ('G', 'T'))
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l1', ('C', 'C'))
  ('s2', 'l2', ('A', 'A'))
  '''
  if not rename_alleles:
    for sample,locus,geno in triples:
      geno = new_genorepr.rep_from_geno(old_genorepr.geno_from_rep(geno))
      yield sample,locus,geno
  else:
    for sample,locus,geno in triples:
      geno = old_genorepr.geno_from_rep(geno)
      if locus in rename_alleles:
        remap = rename_alleles[locus]
        geno  = (remap[geno[0]],remap[geno[1]]) if geno else geno
      geno = new_genorepr.rep_from_geno(geno)
      yield sample,locus,geno


def filter_genomatrixstream_missing(columns,genos):
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
  >>> rows = [('l1',[17,0,51]),
  ...         ('l2',[0,0,0]),
  ...         ('l3',[17,0,0]),
  ...         ('l4',[52,0,68])]
  >>> samples,rows = filter_genomatrixstream_missing(samples,rows)
  >>> samples
  ('s1', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', [17, 51])
  ('l3', [17, 0])
  ('l4', [52, 68])

  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[0,0,0]),
  ...         ('l2',[17,68,51]),
  ...         ('l3',[0,0,0]),
  ...         ('l4',[52,0,68])]
  >>> samples,rows = filter_genomatrixstream_missing(samples,rows)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in rows:
  ...   print row
  ('l2', [17, 68, 51])
  ('l4', [52, 0, 68])
  '''

  # Helper generator function used to filter for completely empty rows.
  # Since the columns_notseen set is invariant to these rows, they may be
  # silently filtered.  This is also desirable because if all columns are
  # seen, then only this simple filter is required to process the remaining
  # rows and the rest need not be materialized.  This can be a huge
  # performance improvement for many datasets with no missing columns.
  def _filter_missing_rows(genos):
    for lname,row in genos:
      if any(row):
        yield lname,row

  genos   = _filter_missing_rows(genos)
  rows    = []

  # Set of column indices not yet observed to have data
  columns_notseen = set(range(len(columns)))

  for lname,row in genos:
    # Store row in materialized list
    rows.append( (lname,row) )

    # Remove any column indices with non-missing data
    columns_notseen.difference_update( [ i for i in columns_notseen if row[i] ] )

    # Stop materializing if there are no more indices to remove
    if not columns_notseen:
      return columns,chain(rows,genos)

  # Full materialize was necessary and some columns need to be filtered
  columns_notseen = pick(columns, sorted(columns_notseen))
  return filter_genomatrixstream_by_column(columns,rows,columns_notseen,exclude=True)


def filter_genotriples_missing(triples):
  '''
  Filter out the genotriple if its genotype is missing.

  @param triples: sequnce of triplets. e.g.
                  [('locus1','sampleA',genotype),
                   ('locus2','sampleB',genotype),...]
  @type  triples: sequence
  @return       : iterator with missing genotype in the triples being filtered out.
  @rtype        : iterator

  >>> ts=[('l1','s1',0),('l2','s2',1)]
  >>> print list(filter_genotriples_missing(ts))
  [('l2', 's2', 1)]
  '''
  return ifilter(itemgetter(2), triples)


def build_genotriples_by_locus(samples,rows):
  '''
  Generate genotype triples from the locus major genotype matrix.
  @param rows: genotype matrix
  @type  rows: sequence
  @rtype:      generator
  @returns:    a genotype triplet stream

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> for s,l,g in build_genotriples_by_locus(samples,rows):
  ...   print s,l,g
  s1 l1 AA
  s2 l1 AG
  s3 l1 GG
  s1 l2 AA
  s2 l2 AT
  s3 l2 TT
  '''
  for locus,genos in rows:
    for sample,geno in izip(samples, genos):
      yield sample,locus,geno


def build_genotriples_by_sample(loci,rows):
  '''
  Generate genotype triples from the sample major genotype matrix.
  @param rows: genotype matrix
  @type  rows: sequence
  @rtype:      generator
  @returns:    a genotype triplet stream

  >>> loci =        ('l1','l2','l3')
  >>> rows = [('s1',['AA','AG','GC']),
  ...         ('s2',['AT','GG','CC'])]
  >>> for s,l,g in build_genotriples_by_sample(loci,rows):
  ...   print s,l,g
  s1 l1 AA
  s1 l2 AG
  s1 l3 GC
  s2 l1 AT
  s2 l2 GG
  s2 l3 CC
  '''
  for sample,genos in rows:
    for locus,geno in izip(loci, genos):
      yield sample,locus,geno


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

  >>> samplemap = dict([('s1','S1'),('s2','S2')])
  >>> locmap    = dict([('l1','L1'),('l2','L2')])
  >>> triples = [('s1','l1','AA'),('s1','l2','AA'),
  ...            ('s2','l1','AA'),('s2','l2','AA')]
  >>> for sample,loc,geno in rename_genotriples(triples,samplemap,locmap):
  ...   print sample,loc,geno
  S1 L1 AA
  S1 L2 AA
  S2 L1 AA
  S2 L2 AA
  '''
  for sample,locus,geno in triples:
    if samplemap is not None:
      sample = samplemap.get(sample,sample)
    if locusmap is not None:
      locus = locusmap.get(locus,locus)
    yield sample,locus,geno


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

  >>> sset = set(['s1'])
  >>> lset = set(['l1','l2'])
  >>> triples = [('s1','l1','AA'),('s1','l2','AA'),
  ...            ('s2','l1','AA'),('s2','l2','AA')]
  >>> for s,l,g in filter_genotriples(triples,sset,lset):
  ...   print s,l,g
  s1 l1 AA
  s1 l2 AA
  >>> lset = None
  >>> for s,l,g in filter_genotriples(triples,sset,lset,exclude=True):
  ...   print s,l,g
  s2 l1 AA
  s2 l2 AA
  '''
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


def rename_genomatrixstream_column(columns,rows,colmap):
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
  >>> samples,rows = rename_genomatrixstream_column(samples,rows,colmap)
  >>> samples
  ('S1', 'S2', 'S3')
  >>> for row in rows:
  ...   print row
  ('l1', ['AA', 'AG', 'GG'])
  ('l2', ['AA', 'AT', 'TT'])
  '''
  columns = tuple(colmap.get(name,name) for name in columns)
  return columns,rows


# FIXME: Optimize trivial case
def filter_genomatrixstream_by_column(columns,rows,colset,exclude=False):
  '''
  Filter the genotype matrix data by a column set.
  Depending on the value of the exclude flag, the column set will
  be used either for inclusion(exclude=False) or exclusion(exclude=True)
  purpose.

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
  >>> samples,rows = filter_genomatrixstream_by_column(samples,rows,colset)
  >>> samples
  ('s1', 's3')
  >>> for row in rows:
  ...   print row
  ('l1', ['AA', 'GG'])
  ('l2', ['AA', 'TT'])

  >>> samples =     ('s1','s2','s3')
  >>> rows = [('l1',['AA','AG','GG']),
  ...         ('l2',['AA','AT','TT'])]
  >>> samples,rows = filter_genomatrixstream_by_column(samples,rows,colset,exclude=True)
  >>> samples
  ('s2',)
  >>> for row in rows:
  ...   print row
  ('l1', ['AG'])
  ('l2', ['AT'])
  '''
  if exclude:
    columns = tuple((name,i) for i,name in enumerate(columns) if name not in colset)
  else:
    columns = tuple((name,i) for i,name in enumerate(columns) if name     in colset)

  if columns:
    columns,indices = izip(*columns)
  else:
    indices = []

  def _filter():
    for locus,genos in rows:
      yield locus,pick(genos,indices)

  return columns,_filter()


def remap_genomatrixstream_column(columns,rows,colmap):
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
  columns,rows = filter_genomatrixstream_by_column(columns,rows,colmap)
  return rename_genomatrixstream_column(columns,rows,colmap)


def reorder_genomatrixstream_columns(columns,rows,labels):
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
  >>> new_samples,new_rows = reorder_genomatrixstream_columns(samples,rows,['s2','s1','s3','s4'])
  >>> new_samples
  ('s2', 's1', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', ['AG', 'AA', 'GG'])
  ('l2', ['AT', 'AA', 'TT'])

  >>> new_samples,new_rows = reorder_genomatrixstream_columns(samples,rows,['s2','s1'])
  >>> new_samples
  ('s2', 's1', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', ['AG', 'AA', 'GG'])
  ('l2', ['AT', 'AA', 'TT'])

  >>> new_samples,new_rows = reorder_genomatrixstream_columns(samples,rows,['s2','s2'])
  Traceback (most recent call last):
       ...
  ValueError: Duplicated column label: s2
  '''
  columnset = set(columns)
  labelset  = set(labels)
  extras    = sorted(l for l in columns if l not in labelset)
  order     = [ l for l in labels if l in columnset ] + list(extras)
  remap     = dict( (l,i) for i,l in enumerate(columns) )

  try:
    indices = [ remap.pop(l) for l in order ]
  except KeyError:
    raise ValueError, 'Duplicated column label: %s' % l

  def _reorder():
    for rowlabel,row in rows:
      yield rowlabel,pick(row, indices)

  return pick(columns,indices),_reorder()


def reorder_genomatrixstream_rows(columns, rows, labels):
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

  >>> loci,new_rows = reorder_genomatrixstream_rows(loci,rows,['s2','s1','s3','s4'])
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in new_rows:
  ...   print row
  ('s2', ['TT', 'GG', 'CC'])
  ('s1', ['AA', 'AG', 'CT'])
  ('s3', ['AA', 'GG', 'TT'])

  >>> loci,new_rows = reorder_genomatrixstream_rows(loci,rows,['s2','s1'])
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in new_rows:
  ...   print row
  ('s2', ['TT', 'GG', 'CC'])
  ('s1', ['AA', 'AG', 'CT'])
  ('s3', ['AA', 'GG', 'TT'])

  >>> loci,new_rows = reorder_genomatrixstream_rows(loci,rows,['s2','s2'])
  >>> loci
  ('l1', 'l2', 'l3')
  >>> for row in new_rows:
  ...   print row
  Traceback (most recent call last):
       ...
  ValueError: Duplicated row label: s2
  '''
  def _reorder():
    data      = list(rows)
    rowset    = set(l for l,data in data)
    labelset  = set(labels)
    extras    = sorted(l for l,data in data if l not in labelset)
    order     = [ l for l in labels if l in rowset ] + list(extras)
    remap     = dict( (l,i) for i,(l,data) in enumerate(data) )

    try:
      indices = [ remap.pop(l) for l in order ]
    except KeyError:
      raise ValueError, 'Duplicated row label: %s' % l

    for i in indices:
      yield data[i]

  return columns,_reorder()


def rename_genomatrixstream_row(columns,rows,rowmap):
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
  >>> rowmap = {'l1':'L1','l2':'L2'}
  >>> samples,new_rows = rename_genomatrixstream_row(samples,rows,rowmap)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('L1', ['AA', 'AG', 'GG'])
  ('L2', ['AA', 'AT', 'TT'])
  '''
  def _rename():
    for rowid,genos in rows:
      yield rowmap.get(rowid,rowid),genos
  return columns,_rename()


def filter_genomatrixstream_by_row(columns,rows,rowset,exclude=False):
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
  >>> rowset = set(['l1'])
  >>> samples,new_rows = filter_genomatrixstream_by_row(samples,rows,rowset)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', ['AA', 'AG', 'GG'])

  >>> samples,new_rows = filter_genomatrixstream_by_row(samples,rows,rowset,exclude=True)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l2', ['AA', 'AT', 'TT'])
  '''
  def _filter():
    if exclude:
      for rowid,genos in rows:
        if rowid not in rowset:
          yield rowid,genos
    else:
      for rowid,genos in rows:
        if rowid in rowset:
          yield rowid,genos

  return columns,_filter()


def remap_genomatrixstream_row(columns,rows,rowmap):
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
  columns,rows = filter_genomatrixstream_by_row(columns,rows,rowmap)
  return rename_genomatrixstream_row(rows,rowmap)


# Does not need to support a genorepr, since only one row is
# materialized at a time.
def transpose_generator(columns, rows, missing=0):
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


def transpose_matrix(columns, rows, genorepr=snp_acgt.pack_reps, missing=0):
  '''
  Transpose a matrix of row labels and row data in memory.

  Requires the input data to be fully materialized, and produces a fully
  materialized transpose.

  @param  columns: matrix column names
  @type   columns: sequence of strs
  @param     rows: sequence of pairs of row labels and row data
  @type      rows: sequence of label and sequence pairs
  @param  columns: sequence of column labels corresponding to each row
  @type   columns: sequence of labels
  @param genorepr: function to convert list genotype strings to desired
                   internal representation
  @type  genorepr: unary function
  @return        : tuple of column labels and generator of list of row labels and row data
  @rtype         :  tuple

  >>> r = [('r1','abc'),
  ...      ('r2','def'),
  ...      ('r3','ghi')]
  >>> rowlabels,c = transpose_matrix(['c1','c2','c3'],r,list)
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

  e = [missing]*n
  newrows = [ genorepr(e) for i in xrange(m) ]

  for i,genos in enumerate(rows):
    for j,geno in enumerate(genos):
      newrows[j][i] = geno

  return rowlabels,zip(columns,newrows)


def _test_genodata():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test_genodata()


###################################################
#                                                 #
# IMPORT GENOIO MODULE FOR BACKWARD COMPATIBILITY #
from genoio import *                              #
#                                                 #
###################################################
