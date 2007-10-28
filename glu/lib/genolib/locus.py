# -*- coding: utf-8 -*-
'''
File:          locus.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-10-15

Abstract:      GLU locus model input and output

Requires:      Python 2.5

Revision:      $Id$
'''

from collections               import defaultdict

from glu.lib.fileutils         import namefile,load_table,get_arg,tryint,parse_augmented_filename
from glu.lib.genolib.genoarray import model_from_alleles


def _get_header(filename,hmap,header,required=False):
  index = hmap.get(header,[])

  if not index:
    if required:
      raise ValueError('Cannot find %s header in locus file %s' % (header,namefile(filename)))
    return None

  elif len(index) != 1:
    raise ValueError('Ambiguous %s header in locus file %s' % (header,namefile(filename)))

  return index[0]


class Locus(object):
  __slots__ = ('name','chromosome','location','model')

  def __init__(self, name, model, chromosome=None, location=None):
    self.name       = name
    self.model      = model
    self.chromosome = chromosome
    self.location   = location


def load_locus_records(filename,extra_args=None,modelcache=None,**kwargs):
  '''
  Load a locus file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      alleles: alleles
  @param      alleles: sequence of objects
  @param  max_alleles: maximum number of alleles allowed. Default is 2
  @type   max_alleles: int or None
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @return            : default_max_alleles, default_alleles, sequence of
                       locus models of the form: locus name, max_alleles, alleles

  >>> from StringIO import StringIO
  >>> l = """LOCUS\\tMAX_ALLELES\\tALLELES\\tCHROMOSOME\\tLOCATION
  ... l1\\t2\\tA,C\\t1\\t0
  ... l2
  ... l3\\t\\tA,G\\t1
  ... l4\\t45\\t\\tX\\t1234
  ... l5\\t45\\tA,B,C,D,E,F,G,H,I,J,K\\tM\\t5541
  ... """
  >>> max_alleles,alleles,loci = load_locus_records(StringIO(l))
  >>> max_alleles
  2
  >>> alleles
  []
  >>> for lname,max_alleles,alleles,chromosome,location in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 2 [] None None
  l3 2 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> max_alleles,alleles,loci = load_locus_records(StringIO(l),max_alleles=4,alleles=['A','C','G','T'])
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles,chromosome,location in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 4 ['A', 'C', 'G', 'T'] None None
  l3 2 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write(l)
  >>> f.flush()
  >>> max_alleles,alleles,loci = load_locus_records(f.name + ':max_alleles=4:alleles=A,C,G,T')
  >>> for lname,max_alleles,alleles,chromosome,location in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 4 ['A', 'C', 'G', 'T'] None None
  l3 2 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> max_alleles,alleles,loci = load_locus_records(':max_alleles=4:alleles=A,C,G,T')
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles,chromosome,location in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  '''
  if isinstance(filename,basestring) and filename.startswith(':'):
    filename = parse_augmented_filename(filename,kwargs)
    assert not filename
  else:
    rows = load_table(filename,want_header=True,extra_args=kwargs,**kwargs)
    header = rows.next()

  default_alleles     = get_arg(kwargs, ['alleles'], None)
  default_max_alleles = int(get_arg(kwargs, ['max_alleles'], 2)) or None

  if kwargs and extra_args is not None:
    extra_args.update(kwargs)
  elif kwargs:
    raise ValueError('Unexpected filename arguments: %s' % str(args))

  if default_alleles is None:
    default_alleles = []
  elif isinstance(default_alleles,basestring):
    default_alleles = [ a.strip() for a in default_alleles.split(',') ]

  if not filename:
    return default_max_alleles,default_alleles,[]

  hmap = defaultdict(list)
  for i,h in enumerate(header):
    hmap[h].append(i)

  lname_index       = _get_header(filename,hmap,'LOCUS',required=True)
  max_alleles_index = _get_header(filename,hmap,'MAX_ALLELES')
  alleles_index     = _get_header(filename,hmap,'ALLELES')
  chromosome_index  = _get_header(filename,hmap,'CHROMOSOME')
  location_index    = _get_header(filename,hmap,'LOCATION')

  def _loci():
    for i,row in enumerate(rows):
      if not row:
        continue

      n = len(row)

      if lname_index >= n:
        lname = ''
      else:
        lname = intern(row[lname_index])

      if not lname:
        raise ValueError('Invalid locus file record %d in %s (blank locus name)' % (i+1,namefile(filename)))

      informative = False

      if max_alleles_index is not None and max_alleles_index<n:
        max_alleles = int(row[max_alleles_index] or 0)
        informative = max_alleles > 0
      else:
        max_alleles = default_max_alleles

      if alleles_index is not None and alleles_index<n:
        alleles = ( a.strip() for a in row[alleles_index].split(',') )
        alleles = [ intern(a) for a in alleles if a ]
      elif informative:
        alleles = []
      else:
        alleles = default_alleles

      max_alleles = max(max_alleles or 0,len(alleles)) or None

      if chromosome_index is not None and chromosome_index<n:
        chromosome = intern(row[chromosome_index].strip())
      else:
        chromosome = None

      if location_index is not None and location_index<n:
        location = int(row[location_index].strip())
      else:
        location = None

      yield lname,max_alleles,alleles,chromosome,location

  return default_max_alleles,default_alleles,_loci()


def load_models(filename,modelcache=None,**kwargs):
  '''
  Return the default model and a sequence of Locus objects from an augmented
  locus description file

  FIXME: Add docstring and doctests
  '''
  if modelcache is None:
    modelcache = {}

  default_max_alleles,default_alleles,loci = load_locus_records(filename,**kwargs)

  defmodel = model_from_alleles(default_alleles,max_alleles=default_max_alleles)

  def _loci():
    for lname,max_alleles,alleles,chromosome,location in loci:
      key   = (max_alleles,alleles)
      model = modelcache.get(key)
      if model is None:
        model = modelcache[key] = model_from_alleles(alleles,max_alleles=max_alleles)

      yield Locus(lname, chromosome, location, model)

  return defmodel,_loci()


def load_modelmap(filename,**kwargs):
  '''
  Create a modelmap dictionary from an augmented locus description file

  FIXME: Add docstring and doctests
  '''
  defmodel,loci = load_models(filename,**kwargs)
  modelmap = defaultdict(lambda: defmodel)
  modelmap.update( (locus.name,locus.model) for locus in loci )
  return modelmap


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
