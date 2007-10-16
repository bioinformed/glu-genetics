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

from glu.lib.fileutils         import namefile,load_table,get_arg,parse_augmented_filename
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


def load_locus_file(filename,extra_args=None,**kwargs):
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
  >>> l = """LOCUS\\tMAX_ALLELES\\tALLELES
  ... l1\\t2\\tA,C
  ... l2
  ... l3\\t\\tA,G
  ... l4\\t45
  ... l5\\t45\\tA,B,C,D,E,F,G,H,I,J,K
  ... """
  >>> max_alleles,alleles,loci = load_locus_file(StringIO(l))
  >>> max_alleles
  2
  >>> alleles
  []
  >>> for lname,max_alleles,alleles in loci:
  ...   print lname,max_alleles,alleles
  l1 2 ['A', 'C']
  l2 2 []
  l3 2 ['A', 'G']
  l4 45 []
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

  >>> max_alleles,alleles,loci = load_locus_file(StringIO(l),max_alleles=4,alleles=['A','C','G','T'])
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles in loci:
  ...   print lname,max_alleles,alleles
  l1 2 ['A', 'C']
  l2 4 ['A', 'C', 'G', 'T']
  l3 2 ['A', 'G']
  l4 45 []
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write(l)
  >>> f.flush()
  >>> max_alleles,alleles,loci = load_locus_file(f.name + ':max_alleles=4:alleles=A,C,G,T')
  >>> for lname,max_alleles,alleles in loci:
  ...   print lname,max_alleles,alleles
  l1 2 ['A', 'C']
  l2 4 ['A', 'C', 'G', 'T']
  l3 2 ['A', 'G']
  l4 45 []
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

  >>> max_alleles,alleles,loci = load_locus_file(':max_alleles=4:alleles=A,C,G,T')
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles in loci:
  ...   print lname,max_alleles,alleles
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
        alleles = [ intern(a.strip()) for a in row[alleles_index].split(',') ]
      elif informative:
        alleles = []
      else:
        alleles = default_alleles

      max_alleles = max(max_alleles or 0,len(alleles)) or None

      yield lname,max_alleles,alleles

  return default_max_alleles,default_alleles,_loci()


def load_modelmap(filename,extra_args=None,**kwargs):
  '''
  Create a modelmap dictionary from an augmented locus description file

  FIXME: Add docstring and doctests
  '''
  default_max_alleles,default_alleles,loci = load_locus_file(filename,extra_args,**kwargs)
  defmodel = model_from_alleles(default_alleles,max_alleles=default_max_alleles)
  modelmap = defaultdict(lambda: defmodel)
  modelmap.update( (lname,model_from_alleles(alleles,max_alleles=max_alleles))
                       for lname,max_alleles,alleles in loci )
  return modelmap


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
