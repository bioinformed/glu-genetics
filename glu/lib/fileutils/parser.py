# -*- coding: utf-8 -*-

__abstract__  = 'file related utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os

from   glu.lib.utils import is_str


__all__ = ['parse_augmented_name','parse_augmented_filename','get_arg',
           'tryfloat','tryint','tryint1','trybool']


def parse_augmented_name(name,args):
  '''
  Retrieve option-value pairs from the name delimited by colon.  Unlike
  parse_augmented_filename, name need not be a file and no filesystem
  checking is performed.

  @param  name: potentially augmented name
  @type   name: str
  @param  args: option dictionary
  @type   args: dict
  @return     : name
  @rtype      : str

  >>> args = {}
  >>> parse_augmented_name('foo',args)
  'foo'
  >>> sorted(args.iteritems())
  []

  >>> args = {}
  >>> parse_augmented_name('foo:',args)
  'foo'
  >>> sorted(args.iteritems())
  []

  >>> args = {}
  >>> parse_augmented_name('hello:foo=bar',args)
  'hello'
  >>> sorted(args.iteritems())
  [('foo', 'bar')]

  >>> args = {}
  >>> parse_augmented_name('goodbye:foo=bar:baz:=spoo',args)
  'goodbye'
  >>> sorted(args.iteritems())
  [('', 'spoo'), ('baz', ''), ('foo', 'bar')]
  '''
  name_save = name

  while ':' in name:
    name,arg = name.rsplit(':',1)
    kv = arg.split('=')

    if len(kv) > 2:
      raise ValueError("Invalid augmented name argument in '%s'" % name_save)
    elif len(kv) == 1:
      kv.append('')

    if kv != ['','']:
      args[kv[0]] = kv[1]

  return name


def parse_augmented_filename(filename,args):
  '''
  Retrieve option-value pairs from the filename delimited by colon.  Options and values are added to the args dictionary.

  @param  filename: potentially augmented file name string
  @type   filename: str
  @param      args: option dictionary
  @type       args: dict
  @return         : file name
  @rtype          : str

  >>> args = {}
  >>> parse_augmented_filename('file.gz',args)
  'file.gz'
  >>> sorted(args.iteritems())
  []

  >>> args = {}
  >>> parse_augmented_filename('file.gz:',args)
  'file.gz'
  >>> sorted(args.iteritems())
  []

  >>> args = {}
  >>> parse_augmented_filename('file.gz:foo=bar',args)
  'file.gz'
  >>> sorted(args.iteritems())
  [('foo', 'bar')]

  >>> args = {}
  >>> parse_augmented_filename('file.gz:foo=bar:baz:=spoo',args)
  'file.gz'
  >>> sorted(args.iteritems())
  [('', 'spoo'), ('baz', ''), ('foo', 'bar')]
  '''
  if not is_str(filename):
    return filename

  filename_save = filename
  while not os.path.exists(filename) and ':' in filename:
    filename,arg = filename.rsplit(':',1)
    kv = arg.split('=')

    if len(kv) > 2:
      raise ValueError("Invalid augmented filename argument in '%s'" % filename_save)
    elif len(kv) == 1:
      kv.append('')

    if kv != ['','']:
      args[kv[0]] = kv[1]

  return filename


def get_arg(args, names, default=None):
  '''
  Retrieve argument which exists in the supplied list

  @param     args: supplied argument list or set
  @type      args: list or set
  @param    names: list or set of potential argument names
  @type     names: list or set
  @param  default: set the default return to None
  @type   default: str or None
  @return        : argument which exists in the supplied list
  @rtype         : str or None
  '''
  for name in names:
    if name in args:
      return args.pop(name)
  return default


def tryfloat(n):
  '''
  Try to coerce an arbitrary object to a float, otherwise return the
  original value.  Existing integer objects are returned as-is without
  conversion to a floating point value.

  @param s: arbitrary item
  @type  s: object
  @return : integer coerced value or the original value
  @rtype  : float or object

  >>> tryfloat(1)
  1
  >>> tryfloat(1.65)==1.65
  True
  >>> tryfloat(1L)
  1L
  >>> tryfloat(None)
  >>> tryfloat('1e-1')==1e-1
  True
  >>> tryfloat(' 1.2 ')
  1.2
  >>> tryfloat([1,2,3])
  [1, 2, 3]
  '''
  if n is None:
    return None
  elif isinstance(n, (int,long,float)):
    return n

  try:
    return float(n)
  except (TypeError,ValueError):
    return n


def tryint(s):
  '''
  Try to coerce an arbitrary object to an integer, otherwise return the
  original value

  @param s: arbitrary item
  @type  s: object
  @return : integer coerced value or the original value
  @rtype  : int or object

  >>> tryint(1)
  1
  >>> tryint(1L)
  1
  >>> tryint(None)
  >>> tryint('1')
  1
  >>> tryint(' 1 ')
  1
  >>> tryint([1,2,3])
  [1, 2, 3]
  '''
  try:
    return int(s)
  except (ValueError,TypeError):
    return s


def tryint1(s):
  '''
  Try to coerce an arbitrary object to an integer, otherwise return the
  original value.  Values provided as strings are assumed to be from a
  1-based indexing scheme, so they are decremented by one upon return.

  @param s: arbitrary item
  @type  s: object
  @return : integer coerced value or the original value
  @rtype  : int or object

  >>> tryint1(1)
  1
  >>> tryint1(1L)
  1
  >>> tryint1(None)
  >>> tryint1('1')
  0
  >>> tryint1(' 1 ')
  0
  >>> tryint1([1,2,3])
  [1, 2, 3]
  >>> tryint1(0)
  0
  >>> tryint1('0')
  Traceback (most recent call last):
     ...
  ValueError: Index must be greater than zero
  '''
  try:
    ss = int(s)
  except (ValueError,TypeError):
    return s

  if isinstance(s, basestring):
    if ss==0:
      raise ValueError('Index must be greater than zero')
    ss -= 1

  return ss


def trybool(s):
  if isinstance(s,str):
    s = s.lower()
    if s in ('t','true','y','yes','1'):
      s = True
    elif s in ('f','false','n','no','0'):
      s = False
  return bool(s)


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
