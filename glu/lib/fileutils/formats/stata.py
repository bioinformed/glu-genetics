# -*- coding: utf-8 -*-

__abstract__  = 'fileutils Stata file support'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

'''
This code is based on the statsmodels project
(http://statsmodels.sourceforge.net/), which is in turn based on code from
the PyDTA package (http://presbrey.mit.edu/PyDTA).
'''

from   struct                   import unpack, calcsize
from   itertools                import izip

from   glu.lib.fileutils.auto   import autofile, hyphen
from   glu.lib.fileutils.parser import parse_augmented_filename, get_arg


__all__ = ['table_reader_stata']


class StataMissingValue(object):
  '''
  An observation's missing value.

  Parameters
  -----------
  offset
  value

  Attributes
  ----------
  string
  value

  Notes
  -----
  More information: <http://www.stata.com/help.cgi?missing>
  '''

  def __init__(self, offset, value):
    self._value = value
    # FIXME[kbj]: The "and '.' or" below is clearly wrong.  I'm not sure
    #             what this is supposed to do.
    if isinstance(value, (int,long)):
      self._str = value-offset is 1 and '.' or ('.' + chr(value-offset+96))
    else:
      self._str = '.'

  string = property(lambda self: self._str,
                    doc="The Stata representation of the missing value: '.', '.a'..'.z'")

  value = property(lambda self: self._value,
                   doc='The binary representation of the missing value.')

  def __str__(self):
    return self._str
  __str__.__doc__ = string.__doc__


class StataVariable(object):
  '''
  A dataset variable.  Not intended for public use.

  Parameters
  ----------
  variable_data

  Attributes
  -----------
  format : str
      Stata variable format.  See notes for more information.
  index : int
      Zero-index column index of variable.
  label : str
      Data Label
  name : str
      Variable name
  type : str
      Stata data type.  See notes for more information.
  value_format : str
      Value format.

  Notes
  -----
  More information: http://www.stata.com/help.cgi?format
  '''
  def __init__(self, variable_data):
      self._data = variable_data

  def __int__(self):
      return self.index

  def __str__(self):
    return self.name

  index = property(lambda self: self._data[0],
                   doc='the variable\'s index within an observation')

  type = property(lambda self: self._data[1],
                  doc='the data type of variable\n\nPossible types are:\n'
                      '{1..244:string, b:byte, h:int, l:long, f:float, d:double)')
  name = property(lambda self: self._data[2], doc='the name of the variable')

  format = property(lambda self: self._data[4],
                    doc='the variable\'s Stata format')

  value_format = property(lambda self: self._data[5],
                   doc='the variable\'s value format')

  label = property(lambda self: self._data[6], doc='the variable\'s label')
  __int__.__doc__ = index.__doc__
  __str__.__doc__ = name.__doc__


class StataReader(object):
  '''
  Stata .dta file reader.

  Provides methods to return the metadata of a Stata .dta file and
  a generator for the data itself.

  Parameters
  ----------
  file : file-like
      A file-like object representing a Stata .dta file.

  missing_values : bool
      If missing_values is True, parse missing_values and return a
      Missing Values object instead of None.

  Notes
  -----

  This is known only to work on file formats 113 (Stata 8/9) and 114 (Stata
  10/11).  It seems to work on file format 111, but this requires additional
  validation.  Needs to be tested on older versions.  Known not to work on
  format 104, 108, 110.

  For more information about the .dta format see
  http://www.stata.com/help.cgi?dta
  http://www.stata.com/help.cgi?dta_113
  '''

  _header = {}
  _data_location = 0
  _col_sizes = ()
  _has_string_data = False
  _missing_values = False
  TYPE_MAP = range(251)+list('bhlfd')
  MISSING_VALUES = { 'b': (-127,100), 'h': (-32767, 32740), 'l':
          (-2147483647, 2147483620), 'f': (-1.701e+38, +1.701e+38), 'd':
          (-1.798e+308, +8.988e+307) }

  def __init__(self, fname, missing_values=False, missing_code=None):
    self._missing_values = missing_values
    self._missing_code = missing_code
    self._parse_header(fname)

  def file_headers(self):
    '''
    Returns all .dta file headers.

    out: dict
        Has keys typlist, data_label, lbllist, varlist, nvar, filetype,
        ds_format, nobs, fmtlist, vlblist, time_stamp, srtlist, byteorder
    '''
    return self._header

  def file_format(self):
    '''
    Returns the file format.

    Returns
    -------
    out : int

    Notes
    -----
    Format 111: Stata ???
    Format 113: Stata 8/9
    Format 114: Stata 10/11
    '''
    return self._header['ds_format']

  def file_label(self):
    '''
    Returns the dataset's label.

    Returns
    ------
    out: string
    '''
    return self._header['data_label']

  def file_timestamp(self):
    '''
    Returns the date and time Stata recorded on last file save.

    Returns
    -------
    out : str
    '''
    return self._header['time_stamp']

  def variables(self):
    '''
    Returns a list of the dataset's StataVariables objects.
    '''
    return map(StataVariable, izip(range(self._header['nvar']),
                                         self._header['typlist'],
                                         self._header['varlist'],
                                         self._header['srtlist'],
                                         self._header['fmtlist'],
                                         self._header['lbllist'],
                                         self._header['vlblist']))

  def dataset(self, as_dict=False):
    '''
    Returns a Python generator object for iterating over the dataset.


    Parameters
    ----------
    as_dict : bool, optional
        If as_dict is True, yield each row of observations as a dict.
        If False, yields each row of observations as a list.

    Returns
    -------
    Generator object for iterating over the dataset.  Yields each row of
    observations as a list by default.

    Notes
    -----
    If missing_values is True during instantiation of StataReader then
    observations with StataMissingValue(s) are not filtered and should
    be handled by your application.
    '''

    try:
      self._file.seek(self._data_location)
    except IOError:
      pass

    if as_dict:
      vars = map(str, self.variables())
      for i in xrange(len(self)):
        yield dict(zip(vars, self._next()))
    else:
      for i in xrange(self._header['nobs']):
        yield self._next()

  ### Python special methods

  def __iter__(self):
    for i in xrange(self._header['nobs']):
        yield self._next()

  def __len__(self):
    '''
    Return the number of observations in the dataset.

    This value is taken directly from the header and includes observations
    with missing values.
    '''
    return self._header['nobs']

  def __getitem__(self, k):
    '''
    Seek to an observation indexed k in the file and return it, ordered
    by Stata's output to the .dta file.

    k is zero-indexed.  Prefer using R.data() for performance.
    '''
    if not (type(k) is int or type(k) is long) or k < 0 or k > len(self)-1:
      raise IndexError(k)
    loc = self._data_location + sum(self._col_size()) * k
    if self._file.tell() != loc:
      self._file.seek(loc)
    return self._next()

  ### Private methods

  def _null_terminate(self, s):
    try:
      return s.lstrip('\x00')[:s.index('\x00')]
    except Exception:
      return s

  def _parse_header(self, file_object):
    self._file = file_object

    # parse headers
    self._header['ds_format'] = unpack('b', self._file.read(1))[0]

    if self._header['ds_format'] not in [111,113,114]:
      raise ValueError('Unsupported Stata format: %d' % self._header['ds_format'])

    byteorder = self._header['byteorder'] = unpack('b',self._file.read(1))[0]==0x1 and '>' or '<'

    self._header['filetype'] = unpack('b', self._file.read(1))[0]
    self._file.read(1)
    nvar = self._header['nvar'] = unpack(byteorder+'h',self._file.read(2))[0]

    if self._header['ds_format'] < 114:
      self._header['nobs'] = unpack(byteorder+'i', self._file.read(4))[0]
    else:
      self._header['nobs'] = unpack(byteorder+'i', self._file.read(4))[0]

    self._header['data_label'] = self._null_terminate(self._file.read(81))
    self._header['time_stamp'] = self._null_terminate(self._file.read(18))

    # parse descriptors
    self._header['typlist'] = [self.TYPE_MAP[ord(self._file.read(1))]    for i in range(nvar)]
    self._header['varlist'] = [self._null_terminate(self._file.read(33)) for i in range(nvar)]
    self._header['srtlist'] = unpack(byteorder+('h'*(nvar+1)),self._file.read(2*(nvar+1)))[:-1]

    if self._header['ds_format'] <= 113:
      self._header['fmtlist'] = [self._null_terminate(self._file.read(12)) for i in range(nvar)]
    else:
      self._header['fmtlist'] = [self._null_terminate(self._file.read(49)) for i in range(nvar)]

    self._header['lbllist'] = [self._null_terminate(self._file.read(33)) for i in range(nvar)]
    self._header['vlblist'] = [self._null_terminate(self._file.read(81)) for i in range(nvar)]

    # ignore expansion fields
    # When reading, read five bytes; the last four bytes now tell you
    # the size of the next read, which you discard.  You then continue
    # like this until you read 5 bytes of zeros.

    # TODO: The way I read this is that they both should be zero, but
    #       that's not what we get.
    while True:
      data_type = unpack(byteorder+'b', self._file.read(1))[0]
      data_len  = unpack(byteorder+'i', self._file.read(4))[0]
      if data_type == 0:
        break
      self._file.read(data_len)

    # other state vars
    try:
      self._data_location = self._file.tell()
    except IOError:
      self._data_location = None

    self._has_string_data = len(filter(lambda x: type(x) is int, self._header['typlist'])) > 0
    self._col_size()

  def _calcsize(self, fmt):
    return isinstance(fmt,int) and fmt or calcsize(self._header['byteorder']+fmt)

  def _col_size(self, k=None):
    '''Calculate size of a data record.'''
    if len(self._col_sizes) == 0:
      self._col_sizes = map(lambda x: self._calcsize(x),self._header['typlist'])
    if k is None:
      return self._col_sizes
    else:
      return self._col_sizes[k]

  def _unpack(self, fmt, byt):
    d = unpack(self._header['byteorder']+fmt, byt)[0]
    if fmt[-1] in self.MISSING_VALUES:
      nmin, nmax = self.MISSING_VALUES[fmt[-1]]
      if d < nmin or d > nmax:
        if self._missing_values:
          return StataMissingValue(nmax, d)
        else:
          return self._missing_code
    return d

  def _next(self):
    typlist = self._header['typlist']
    if self._has_string_data:
      data = [None]*self._header['nvar']
      for i in range(len(data)):
        if type(typlist[i]) is int:
          data[i] = self._null_terminate(self._file.read(typlist[i]))
        else:
          data[i] = self._unpack(typlist[i], self._file.read(self._col_size(i)))
      return data
    else:
      return map(lambda i: self._unpack(typlist[i],self._file.read(self._col_size(i))),range(self._header['nvar']))


def table_reader_stata(filename, extra_args=None, **kwargs):
  '''
  Return a configured delimited table reader
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name = parse_augmented_filename(filename,args)
  hyin = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  dta = StataReader(autofile(hyphen(name,hyin)), missing_code='')

  def _rows(dta):
    header = dta.file_headers()
    yield header['varlist']
    for row in dta:
      yield map(str,row)

  return _rows(dta)
