# -*- coding: utf-8 -*-

__abstract__  = 'genotype transformation objects'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import optparse

from   types             import NoneType
from   collections       import defaultdict
from   itertools         import izip

from   glu.lib.utils     import as_set
from   glu.lib.fileutils import namefile,get_arg,trybool,load_list,load_map,load_table


seq_type = (NoneType,set,dict,list,tuple)
map_type = (NoneType,dict)


def load_rename_alleles_file(filename):
  '''
  Load an allele renameing file

  @param filename: a file name or file object
  @type  filename: str or file object

  >>> from StringIO import StringIO
  >>> data = StringIO('l1\\tA,C,G,T\\tT,G,C,A\\nl3\\tA\\tC\\nl5\\tA,B\\tC,T')
  >>> for lname,alleles in sorted(load_rename_alleles_file(data).iteritems()):
  ...   print lname,sorted(alleles.iteritems())
  l1 [(None, None), ('A', 'T'), ('C', 'G'), ('G', 'C'), ('T', 'A')]
  l3 [(None, None), ('A', 'C')]
  l5 [(None, None), ('A', 'C'), ('B', 'T')]
  '''
  rows = load_table(filename)

  rename = {}
  for i,row in enumerate(rows):
    if not row:
      continue

    if len(row) < 3:
      raise ValueError('Incomplete allele rename record %d in %s' % (i+1,namefile(filename)))

    lname,old_alleles,new_alleles = row[:3]

    lname = intern(lname.strip())

    old_alleles = [ intern(a.strip()) for a in old_alleles.split(',') ]
    new_alleles = [ intern(a.strip()) for a in new_alleles.split(',') ]

    if not lname and not old_alleles and not new_alleles:
      continue

    if not lname:
      raise ValueError('Missing locus name in allele rename record %d in %s' % (i+1,namefile(filename)))

    if len(old_alleles) != len(new_alleles) or '' in old_alleles or '' in new_alleles:
      raise ValueError('Invalid allele rename record %d for %s in %s' % (i+1,lname,namefile(filename)))

    locus_rename = dict( izip(old_alleles,new_alleles) )
    locus_rename[None] = None

    if lname in rename and rename[lname] != locus_rename:
      raise ValueError('Inconsistent rename record %d for %s in %s' % (i+1,lname,namefile(filename)))

    rename[lname] = locus_rename

  return rename


class GenoTransform(object):
  '''
  Create a GenoTransform object to specify various transformation on the genodata.
  Supported operations: include/exclude/rename samples or loci; optional filter to remove missing genotypes
  '''
  def __init__(self, include_samples, exclude_samples, rename_samples, order_samples,
                     include_loci,    exclude_loci,    rename_loci,    order_loci,
                     recode_models=None, rename_alleles=None, repack=False, filter_missing=False):
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
    @param   recode_models: recode using a new genome descriptor.  Default is None
    @type    recode_models: Genome instance or None
    @param  rename_alleles: rename alleles for any loci in the supplied dictionary from old allele name to new allele name
    @type   rename_alleles: dict from old_allele str -> new_allele str
    @param  filter_missing: filter missing genotypes from the stream
    @type   filter_missing: bool
    @param          repack: trigger repacking of genotypes to ensure that the most compact storage
                            method is used
    @type           repack: bool
    @return               : transformed genotriple stream
    @rtype                : GenotripleStream
    '''
    self.samples = GenoSubTransform(include_samples, exclude_samples, rename_samples, order_samples)
    self.loci    = GenoSubTransform(include_loci,    exclude_loci,    rename_loci,    order_loci)

    if not isinstance(rename_alleles, map_type):
      rename_alleles = load_rename_alleles_file(rename_alleles)

    self.recode_models            = recode_models
    self.rename_alleles           = rename_alleles
    self.repack                   = repack
    self.filter_missing_genotypes = filter_missing

  @staticmethod
  def from_options(options):
    '''
    Create a new GenoTransform object from command line option list

    @return: transformed genotriple stream
    @rtype : GenotripleStream
    '''
    return GenoTransform(include_samples=getattr(options,'includesamples',None),
                         exclude_samples=getattr(options,'excludesamples',None),
                          rename_samples=getattr(options,'renamesamples', None),
                           order_samples=getattr(options,'ordersamples',  None),
                            include_loci=getattr(options,'includeloci',   None),
                            exclude_loci=getattr(options,'excludeloci',   None),
                             rename_loci=getattr(options,'renameloci',    None),
                              order_loci=getattr(options,'orderloci',     None),
                          rename_alleles=getattr(options,'renamealleles', None),
                          filter_missing=getattr(options,'filtermissing', None))

  @staticmethod
  def from_kwargs(extra_args=None,**kwargs):
    '''
    Create a new GenoTransform object from key word arguments

    @return: transformed genotriple stream
    @rtype : GenotripleStream
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    transform = GenoTransform(include_samples=get_arg(args,['include_samples','includesamples']),
                              exclude_samples=get_arg(args,['exclude_samples','excludesamples']),
                               rename_samples=get_arg(args,['rename_samples', 'renamesamples' ]),
                                order_samples=get_arg(args,['order_samples',  'ordersamples'  ]),
                                 include_loci=get_arg(args,['include_loci',   'includeloci'   ]),
                                 exclude_loci=get_arg(args,['exclude_loci',   'excludeloci'   ]),
                                  rename_loci=get_arg(args,['rename_loci',    'renameloci'    ]),
                                   order_loci=get_arg(args,['order_loci',     'orderloci'     ]),
                                recode_models=get_arg(args,['recode_models',  'recodemodels'  ]),
                               rename_alleles=get_arg(args,['rename_alleles', 'renamealleles' ]),
                                       repack=trybool(get_arg(args,['repack'])),
                               filter_missing=trybool(get_arg(args,['filter_missing','filtermissing'])))

    if extra_args is None and args:
      raise ValueError("'%s' is an invalid keyword argument for this function" % kwargs.popitem()[0])

    return transform


  @staticmethod
  def from_object(transform):
    if isinstance(transform, optparse.Values) or transform.__class__ is optparse.Values:
      transform = GenoTransform.from_options(transform)
    elif isinstance(transform, dict):
      transform = GenoTransform.from_kwargs(transform)
    elif not isinstance(transform, GenoTransform):
      raise ValueError('Invalid genotype transformation specification')

    return transform


  def merge(self, other):
    samples = self.samples.merge(other.samples)
    loci    = self.loci.merge(other.loci)

    rename_alleles           = _merge_rename_alleles(self.rename_alleles, other.rename_alleles)
    recode_models            = _merge_genome(self.recode_models, other.recode_models)
    rename_alleles           = _merge_map(self.rename_alleles, other.rename_alleles)
    repack                   = self.repack or other.repack
    filter_missing_genotypes = self.filter_missing_genotypes or other.filter_missing_genotypes

    return GenoTransform(include_samples=samples.include,
                         exclude_samples=samples.exclude,
                          rename_samples=samples.rename,
                           order_samples=samples.order,
                            include_loci=loci.include,
                            exclude_loci=loci.exclude,
                             rename_loci=loci.rename,
                              order_loci=loci.order,
                           recode_models=recode_models,
                          rename_alleles=rename_alleles,
                                  repack=repack,
                          filter_missing=filter_missing_genotypes)


def _merge_rename_alleles(amap1, amap2):
  if amap1 is None:
    return amap2
  if amap2 is None:
    return amap1

  raise NotImplementedError


def _merge_genome(genome1, genome2):
  if genome1 is None:
    return genome2
  if genome2 is None:
    return genome1

  raise NotImplementedError


def _merge_map(map1, map2):
  if map1 is None:
    return map2
  if map2 is None:
    return map1

  raise NotImplementedError


def _merge_set(set1,set2):
  if set1 is None:
    return set2
  if set2 is None:
    return set1

  raise as_set(set1)|as_set(set2)


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
    @param   order: sort order, either 'samples', 'locus'
    @type    order: str
    '''
    if not isinstance(include, seq_type):
      include = set(load_list(include))
    elif include is not None:
      include = as_set(include)

    if not isinstance(exclude, seq_type):
      exclude = set(load_list(exclude))
    elif exclude is not None:
      exclude = as_set(exclude)

    if not isinstance(rename, map_type):
      rename = load_map(rename)

    if not isinstance(order, seq_type):
      order = load_list(order)

    self.include = include
    self.exclude = exclude
    self.rename  = rename
    self.order   = order

  def merge(self, other):
    include = _merge_set(self.include, other.include)
    exclude = _merge_set(self.exclude, other.exclude)
    rename  = _merge_map(self.rename,  other.rename)
    order   = self.order or other.order

    return GenoSubTransform(include, exclude, rename, order)


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
    items = as_set(items) & as_set(transform.include)
  elif transform.include is not None:
    items = as_set(transform.include)

  # Construct the minimal sample reverse map by removing excluded items
  # to verify that no two map to the same identifier.
  reverse_map = defaultdict(set)
  for item in items:
    renamed = transform.rename.get(item,item)
    reverse_map[renamed].add(item)

  # Mapping is unique if and only if all reverse map values are unique
  return all( len(v)<=1 for v in reverse_map.itervalues() )


def prove_unique_transform(transform=None,samples=None,loci=None,unique=False):
  '''
  Prove uniqueness of transformation operations

  @param transform: transformation object (optional)
  @type  transform: GenoTransform object
  @param   samples: optional set of samples refered to by the triples
  @type    samples: sequence, set, or None
  @param      loci: optional set of loci refered to by the triples
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


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
