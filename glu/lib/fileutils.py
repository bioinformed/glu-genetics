# -*- coding: utf-8 -*-
'''
File:          fileutils.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-02-21

Abstract:      file related utility functions

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys
import csv

from   itertools import islice

__all__ = ['autofile','namefile','hyphen',
           'guess_format','related_file','guess_related_file',
           'load_list','load_map','load_table']


COMPRESSED_SUFFIXES = set(['gz','Z'])


# Create more standard aliases for Python CSV module dialects.  'excel' is
# now available as 'csv' and 'excel-tab' is now available as 'tsv'.
csv.register_dialect('csv', csv.get_dialect('excel'))
csv.register_dialect('tsv', csv.get_dialect('excel-tab'))


def OSGzipFile_popen(filename, mode):
  '''
  Spawn gzip process and connect to input/output pipes

  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @return         : file object to read from or write to
  @rtype          : file object

  '''
  if os.environ.get('GLU_NOSPAWN'):
    raise OSError('Spawning external processes disabled by GLU_NOSPAWN')

  # Shell escaping is out of the scope of this function, so ensure that no
  # illegal characters appear
  if set('\'"\\') & set(filename):
    raise OSError('Invalid characters in filename for this method')

  gzipexe = os.environ.get('GLU_GZIP','gzip')

  #FIXME: Add additional sanity checks to the executable name
  if set('\'"\\;') & set(gzipexe):
    raise OSError('Invalid characters in gzip executable')

  if 'w' in mode:
    f = os.popen("%s -c > '%s'" % (gzipexe,filename), 'w', 10240)
  else:
    f = os.popen("%s -d -c '%s'" % (gzipexe,filename), 'r', 10240)

  return f


def OSGzipFile_subprocess(filename, mode):
  '''
  Spawn a subprocess to run gzip and connect to input/output pipes

  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @return         : file object to read from or write to
  @rtype          : file object
  '''

  from subprocess import Popen,PIPE

  if os.environ.get('GLU_NOSPAWN'):
    raise OSError('Spawning external processes disabled by GLU_NOSPAWN')

  gzipexe = os.environ.get('GLU_GZIP','gzip')

  #FIXME: Add additional sanity checks to the executable name
  if set('\'"\\;') & set(gzipexe):
    raise OSError('Invalid characters in gzip executable')

  if 'w' in mode:
    out = file(filename,mode)
    cmd = [gzipexe,'-c']
    f   = Popen(cmd, stdin=PIPE, stdout=out).stdin
  else:
    cmd = [gzipexe,'-d','-c',filename]
    f   = Popen(cmd, stdout=PIPE, universal_newlines='U' in mode).stdout

  return f


def hyphen(filename,defaultfile):
  '''
  Return the default if input is '-', otherwise itself

  @param  defaultfile: default file name or file object
  @type   defaultfile: str or file object
  @param     filename: file name or file object
  @type      filename: str or file object
  @return            : file name or file object to read from or write to
  @rtype             : str or file object
  '''
  if filename == '-':
    return defaultfile
  else:
    return filename


def isstr(s):
  '''
  Return whether the input is a string

  @return   : indicate if the input is a string
  @rtype    : bool
  '''
  return isinstance(s, basestring)


def namefile(filething):
  '''
  Return a human-comprehensible name for the file or file object provided.
  Recognizes file objects with the 'name' attribute, including sys.stdin,
  sys.stdout, and sys.stderr.

  @param  filething: file or filename
  @type   filething: file object or string
  @return:           human-comprehensible file name

  >>> namefile(file('/dev/null'))
  '/dev/null'
  >>> namefile(sys.stdin)
  '<stdin>'
  >>> namefile(sys.stderr)
  '<stderr>'
  >>> namefile('/dev/null')
  '/dev/null'
  '''
  if isstr(filething):
    return filething
  elif getattr(filething, 'name', None) is not None:
    return filething.name
  else:
    return repr(filething)


_autofile_errors = (ImportError,ValueError,OSError)


def compressed_filename(filename):
  '''
  Determine if the input file is in or needs to be in compressed format

  @param  filename: file name or file object
  @type   filename: str or file object
  @return         : suffix of the filename if the filename indicates a compress format,
                    otherwise, empty string
  @rtype          : str

  >>> compressed_filename('subjects.sdat')
  ''
  >>> compressed_filename('subjects.sdat.gz')
  'gz'
  >>> compressed_filename('../subjects.sdat')
  ''
  >>> compressed_filename('../subjects.sdat.gz')
  'gz'
  '''
  if not isstr(filename):
    return ''

  filename = os.path.expanduser(filename)
  parts = os.path.basename(filename).split('.')
  if parts and parts[-1] in COMPRESSED_SUFFIXES:
    return parts[-1]
  return ''


def autofile(filename, mode='r'):
  '''
  Return a file object in the correct compressed format as specified, which is ready to read from or write to

  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @return         : file object to read from or write to
  @rtype          : file object
  '''
  # Pass non-string filename objects back as file-objects
  if not isstr(filename):
    return filename

  filename = os.path.expanduser(filename)
  if compressed_filename(filename):
    try:
      f = OSGzipFile_popen(filename, mode)
    except _autofile_errors:
      try:
        f = OSGzipFile_subprocess(filename, mode)
      except _autofile_errors:
        import gzip
        f = gzip.GzipFile(filename, mode)

  else:
    f = file(filename, mode)

  return f


def guess_format(filename, formats):
  '''
  Return the guessed format indicated by the filename itself

  @param  filename: file name
  @type   filename: str
  @param    format: any expected file format
  @type     format: list of strs
  @return         : file format guessed from the filename
  @rtype          : str

  >>> f = ['sdat','ldat']
  >>> guess_format('subjects.sdat', f)
  'sdat'
  >>> guess_format('subjects.sdat.gz', f)
  'sdat'
  >>> guess_format('../subjects.sdat', f)
  'sdat'
  >>> guess_format('../subjects.sdat.gz', f)
  'sdat'
  >>> guess_format('../subjects.sdat.gz:format=ldat', f)
  'ldat'
  '''

  if not isstr(filename):
    return None

  # Parse to remove augmented arguments
  args = {}
  filename = parse_augmented_filename(filename,args)

  if 'format' in args:
    return args['format']

  parts = os.path.basename(filename).split('.')

  if parts and parts[-1] in COMPRESSED_SUFFIXES:
    parts.pop()

  if parts and parts[-1] in formats:
    return parts[-1]

  return None


def related_file(filename,extension):
  '''
  Return a filename with the extension provided

  @param   filename: base filename with extension
  @type    filename: str
  @param extensions: new extension (without a '.')
  @type  extensions: str
  @return          : new filename
  @rtype           : str

  >>> related_file('foo.ldat', 'sdat')
  'foo.sdat'
  >>> related_file('foo', 'sdat')
  'foo.sdat'
  '''
  if not isstr(filename):
    return None

  prefix = os.path.splitext(filename)[0]

  if not prefix:
    raise ValueError('invalid filename')

  return '%s.%s' % (prefix,extension)


def guess_related_file(filename,extensions):
  '''
  Find a related file with the same name except different prefix.  Only
  files that exist are returned.

  @param   filename: base filename with extension
  @type    filename: str
  @param extensions: list of alternative extensions (without a '.')
  @type  extensions: list of str
  @return          : new filename
  @rtype           : str

  >>> guess_related_file('fileutils.dat',['py'])
  'fileutils.py'
  '''
  if not isstr(filename):
    return None

  prefix,ext = os.path.splitext(filename)

  if not prefix:
    raise ValueError('invalid filename')

  for new_ext in extensions:
    testfile = '%s.%s' % (prefix,new_ext)
    if os.path.isfile(testfile):
      return testfile

  return None


def parse_augmented_filename(filename,args):
  '''
  Retrieve option-value pairs from the filename delimited by colon

  @param  filename: file name
  @type   filename: str
  @return         : filename and tuples of appended additional option and its value
  @rtype          : list

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
  if not isstr(filename):
    return filename

  filename_save = filename
  while not os.path.exists(filename) and ':' in filename:
    filename,arg = filename.rsplit(':',1)
    kv = arg.split('=')

    if len(kv) > 2:
      raise ValueError("Invalid extended filename argument in '%s'" % filename_save)
    elif len(kv) == 1:
      kv.append('')

    if kv != ['','']:
      args[kv[0]] = kv[1]

  return filename


DIALECT_KWARGS = ['dialect','delimiter','doublequote','escapechar','lineterminator',
                  'quotechar','quoting','skipinitialspace','strict']


_unescape_chars = [('\\r','\r'),('\\t','\t'),('\\n','\n'),('\\\\','\\')]

def unescape(s):
  for e,c in _unescape_chars:
    s = s.replace(e,c)
  return s


def tryint(s):
  try:
    return int(s)
  except ValueError:
    return s


def get_csv_dialect(args,default_dialect='tsv'):
  '''
  Retrieve the standardized csv argument list

  @param             args: csv arguments
  @type              args: dict
  @param  default_dialect: set the default value for dialect to 'tsv'
  @type   default_dialect: str
  @return                : list of tuples of (option,value)
  @rtype                 : list of tuples

  >>> args = dict(dialect='csv',quoting='none',skipinitialspace='y',doublequote='false',unrelated='baz',lineterminator='\\n',
  ...             delimiter='\\t')
  >>> sorted(get_csv_dialect(args).iteritems())
  [('delimiter', '\\t'), ('dialect', 'csv'), ('doublequote', False), ('lineterminator', '\\n'), ('quoting', 3), ('skipinitialspace', True)]
  >>> sorted(args.iteritems())
  [('unrelated', 'baz')]
  '''
  dargs = {'dialect':default_dialect}
  for darg in DIALECT_KWARGS:
    if darg in args:
      dargs[darg] = unescape(args.pop(darg))

  if 'quoting' in dargs:
    q = 'QUOTE_' + dargs['quoting'].upper()
    dargs['quoting'] = getattr(csv, q, q)

  # Parse truth values
  for tf in ('skipinitialspace','doublequote'):
    if tf in dargs:
      s = dargs[tf].lower()
      if s in ('t','true','y','1'):
        s = True
      elif s in ('f','false','n','0'):
        s = False
      else:
        raise ValueError
      dargs[tf] = s

  return dargs


def _literal_list(filename):
  return isstr(filename) and not os.path.exists(filename) and filename.startswith(':')


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


def load_list(filename, extra_args=None, **kwargs):
  '''
  Load list of values from a file or literal list. By default, the specified
  file is interpreted as a tab-delimited file with only the first field
  considered.  Blank lines are skipped, whitespace is stripped from the
  beginning and end of each field selected.  These settings can be modified
  via function arguments or by appending additional options to the filename.

  The dialect argument can be used to specify a Python csv module dialect
  name or Dialect object.  In addition, the following Python csv module
  options can be specified by appending information to the specified
  filename: dialect, delimiter, doublequote, escapechar, lineterminator,
  quotechar, and quoting.  The syntax is
  ':option=value', appended to the end of the filename.

  A skip parameter can be specified is used to ignore a certain number of
  lines (e.g., headers) at the beginning of the file.  It can be specified
  as either a function argument or by appending ':skip=n', where n is an
  integer, to the filename.

  The index of the element to be considered can be overridden with the index
  parameter (0 based) or by appending ':index=i' to the end of the filename
  (the latter having higher precedence).  Simple parsing of column headers
  is supported.  If i is not an integer then is assumed to be a column
  heading name in the first row of data (after skipping rows) No escaping or
  quoting is allowed.

  If the filename is a string that begins with ':', it is interpreted as a
  literal list of comma separated values, ignoring the skip and index
  parameters.  No escaping or quoting is allowed.

  The following parameters and aliases are accepted as part of the augmented
  filename:
       skip,  s: number of header lines to skip
       index, i: index of field to select or name of header
        dialect: csv module dialect name (typically 'csv' or 'tsv')
      delimiter: single field delimiter character
    doublequote: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
     escapechar: one-character string that removes any special meaning
                 from the following character.
      quotechar: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
        quoting: Controls when quotes characters are recognised.  It can
                 take on any of the cav QUOTE_ constants (without the QUOTE_
                 prefix).

  @param         filename: file name, file object, or literal string
  @type          filename: str or file object
  @param             skip: number of header lines to skip
  @type              skip: int
  @param            index: index of field to select or name of header.
  @type             index: int or str
  @param          dialect: csv module dialect name or dialect object
  @type           dialect: str or csv.Dialect
  @param       extra_args: optional dictionary to store extraneous arguments, instead of
                           raising an error.
  @type        extra_args: dict
  @return:                 list of values selected
  @rtype:                  list or str

  >>> from StringIO import StringIO
  >>> load_list(StringIO("loc1\\tsample1\\nloc2\\n\\n loc3\\n"))
  ['loc1', 'loc2', 'loc3']
  >>> load_list(StringIO("loc1,sample1\\nloc2\\n\\n loc3\\n"),dialect='csv')
  ['loc1', 'loc2', 'loc3']
  >>> load_list(StringIO("loc1\\tsample1\\nloc2\\n\\n loc3\\n"),index=1)
  ['sample1']
  >>> load_list(':loc1, loc2, loc3')
  ['loc1', 'loc2', 'loc3']
  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1\\tH2\\nv11\\tv12\\nv21\\n\\tv32')
  >>> f.flush()
  >>> load_list(f.name + ':skip=1')
  ['v11', 'v21']
  >>> load_list(f.name + ':skip=1:index=0')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=\\t:index=H1')
  ['v11', 'v21']
  >>> load_list(f.name + ':dialect=excel-tab:skip=1:index=1')
  ['v12', 'v32']
  >>> load_list(f.name + ':index=H2')
  ['v12', 'v32']
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1,H2\\nv11,v12\\nv21\\n,v32')
  >>> f.flush()
  >>> load_list(f.name + ':delimiter=,:skip=1')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=,:skip=1:index=0')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=,:index=H1')
  ['v11', 'v21']
  >>> load_list(f.name + ':dialect=excel:skip=1:index=1')
  ['v12', 'v32']
  >>> load_list(f.name + ':dialect=excel:index=H2')
  ['v12', 'v32']
  '''
  # Process literal syntax
  if _literal_list(filename):
    return [ intern(l.strip()) for l in filename[1:].split(',') if l ]

  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  # Otherwise, handle file name or file object cases
  name      = parse_augmented_filename(filename,args)
  dialect   = get_csv_dialect(args, 'tsv')
  skip      = int(get_arg(args,     ['skip', 's'], 0))
  index     = tryint(get_arg(args,  ['index','i'], 0))

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lines = csv.reader(autofile(name),**dialect)

  if skip:
    lines = islice(lines,skip,None)

  # Resolve column heading name into indicex
  if isstr(index):
    try:
      header = map(str.strip,lines.next())
    except StopIteration:
      return []

    try:
      index = header.index(index)
    except ValueError:
      raise ValueError("Cannot find header '%s'" % index)

  values = (intern(l[index].strip()) for l in lines if len(l)>index)
  return [ v for v in values if v ]


# Anonymous object
_nothing = object()

def load_map(filename,unique=True,extra_args=None,**kwargs):
  '''
  Creates a dictionary representing a list or mapping from a text file.
  Valid files should be composed of standard ASCII lines of text with tab
  delimited fields.  Only the first and second fields of each line are
  considered. Whitespace is stripped from the beginning and end of every
  field considered.  If the skip parameter is used to ignore a certain
  number of lines (e.g., headers) at the beginning of the file.  A default
  parameter may be specified to assign values to keys with empty or
  non-existant value fields.  Otherwise, the value will be set equal to the
  key.

  If the file or filename is a string that begins with ':', it is
  interpreted as a literal list of comma separated values, ignoring the skip
  and index parameters.  No escaping or quoting is allowed.

  The dialect argument can be used to specify a Python csv module dialect
  name or Dialect object.  In addition, the following Python csv module
  options can be specified by appending information to the specified
  filename: dialect, delimiter, doublequote, escapechar, lineterminator,
  quotechar, and quoting.  The syntax is ':option=value', appended to the
  end of the filename.

  A skip parameter can be specified is used to ignore a certain number of
  lines (e.g., headers) at the beginning of the file.  It can be specified
  as either a function argument or by appending ':skip=n', where n is an
  integer, to the filename.

  The index of the key and value elements can be overridden with the
  key_index and value_index parameters (0 based) or by appending ':key=i'
  and/or ':value=j' to the end of the filename.  Simple parsing of column
  headers is supported.  If either key or value is not an integer then it is
  assumed to be a column heading name in the first row of data (after
  skipping rows). No escaping or quoting is allowed.

  The default behavior (unique=True) requires that all mappings are unique
  and consistent in the sense that every key must map to exactly one value,
  though that mapping may appear more than once in the input.  Otherwise a
  ValueError is raised.  If unique=False, then a dictionary of key to list
  of values is returned.

  Each line is mapped to key,value pairs that are returned in a dictionary
  and treated as follows:

  1) Empty line             --> skipped
  2) Empty key field        --> skipped
  3) Only key field         --> identity mapping, key=key field, value=key field
                                if not default, otherwise default
  4) Empty value field      --> identity mapping, key=key field, value=key field
                                if not default, otherwise default
  5) Non-empty value field  --> key=key field, value=value field
  6) More than two fields   --> additional fields are ignored

  NOTE: This file format is best suited to representing inclusion lists and
  mapping from one value to another.  Due to the way the data are
  parsed, it is impossible to obtain an implicit exclusion list using
  null values.  This is because the second key is always set a
  non-empty string.

  The following parameters and aliases are accepted as part of the augmented
  filename:
       skip,  s: number of header lines to skip
   key_index, i: index of field to select or name of header for keys
 value_index, v: index of field to select or name of header for values
  default,def,d: default value for keys with no or empty value
        dialect: csv module dialect name ('csv' or 'tsv')
      delimiter: single field delimiter character
    doublequote: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
     escapechar: one-character string that removes any special meaning
                 from the following character.
      quotechar: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
        quoting: Controls when quotes characters are recognised.  It can
                 take on any of the cav QUOTE_ constants (without the QUOTE_
                 prefix).

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       unique: require key value pairs to be unique and consistent
  @type        unique: bool
  @param         skip: number of initial lines to be skipped (defaults to 0)
  @type          skip: int
  @param      default: default value for keys with no or empty value, if not
                       specified the value is taken to be the key
  @type       default: object
  @param      dialect: csv module dialect name or dialect object
  @type       dialect: str or csv.Dialect
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @return            : mapping from key string to value string, if
                       unique=True.  Otherwise, mapping from key string to
                       list of value strings.
  @rtype             : dict

  >>> def test(f,**kw): return sorted(load_map(f,**kw).iteritems())
  >>> from StringIO import StringIO
  >>> test(StringIO("loc1\\tlocA \\nloc2\\n\\n\\t\\n loc3 \\tlocC\\n\\tfoo\\nloc4\\t\\n"))
  [('loc1', 'locA'), ('loc2', 'loc2'), ('loc3', 'locC'), ('loc4', 'loc4')]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"))
  [('loc1', 'loc1'), ('loc2', 'loc2')]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tfoo\\nloc2"),unique=False)
  [('loc1', ['loc1', 'foo']), ('loc2', ['loc2'])]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tfoo\\nloc2")) # doctest: +ELLIPSIS
  Traceback (most recent call last):
       ...
  ValueError: Found 1 non-unique mapping in <StringIO.StringIO instance at ...>. The first is "loc1" -> "loc1" and "foo"
  >>> test(StringIO("loc1\\tloc1\\nloc2"),default='missing')
  [('loc1', 'loc1'), ('loc2', 'missing')]
  >>> test(':loc1,loc2')
  [('loc1', 'loc1'), ('loc2', 'loc2')]
  >>> test(':loc1=locA ,loc2,,=, loc3 =locC,=foo,loc4=')
  [('loc1', 'locA'), ('loc2', 'loc2'), ('loc3', 'locC'), ('loc4', 'loc4')]
  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1\\tH2\\nv11\\tv12\\nv21\\n\\tv32')
  >>> f.flush()
  >>> test(f.name + ':skip=1')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':skip=1:key=0:value=1')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':skip=1:key=0')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':key=H1:v=H2')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':value=H2')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':k=H1')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':skip=1:key=1:value=0')
  [('v12', 'v11'), ('v32', 'v32')]
  >>> test(f.name + ':key=H2:value=H1')
  [('v12', 'v11'), ('v32', 'v32')]
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  # Handle :literal syntax
  if _literal_list(filename):
    key_index,value_index = 0,1
    lines = (kv.split('=') for kv in filename[1:].split(',') if kv)
    default = get_arg(args, ['default'], _nothing)

  else:
    name        = parse_augmented_filename(filename,args)
    dialect     = get_csv_dialect(args,'tsv')
    skip        = int(get_arg(args, ['s','skip'], 0))
    default     = get_arg(args, ['d','def','default'], _nothing)
    key_index   = tryint(get_arg(args, ['k','key',  'key_index'  ], 0))
    value_index = tryint(get_arg(args, ['v','value','value_index'], 1))

    lfile = autofile(name)
    lines = csv.reader(lfile, **dialect)

    if skip:
      lines = islice(lines,skip,None)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  # Resolve column heading names into indices
  if isstr(key_index) or isstr(value_index):
    try:
      header = map(str.strip,lines.next())
    except StopIteration:
      return {}

    if isstr(key_index):
      try:
        key_index = header.index(key_index)
      except ValueError:
        raise ValueError("Cannot find key header '%s'" % key_index)

    if isstr(value_index):
      try:
        value_index = header.index(value_index)
      except ValueError:
        raise ValueError("Cannot find value header '%s'" % value_index)

  # Parse the data file
  def _load_map_generator():
    for row in lines:
      n = len(row)

      key = intern(row[key_index].strip()) if key_index<n else None

      if not key:
        continue

      value = intern(row[value_index].strip()) if value_index<n else None

      if not value:
        if default is not _nothing:
          value = default
        else:
          value = key

      yield key,value

  m = {}
  unique_fails = []
  for key,value in _load_map_generator():
    if unique:
      fail = m.get(key,value) != value
      if fail:
        unique_fails.append( (key,value) )
      else:
        m[key] = value
    else:
      l = m.setdefault(key,[])
      if value not in l:
        l.append(value)

  if unique_fails:
    n = len(unique_fails)
    key,value = unique_fails[0]
    raise ValueError('Found %d non-unique mapping in %s. The first is "%s" -> "%s" and "%s"'
                       % (n,namefile(filename),key,m[key],value))

  return m


def get_arg_columns(args):
  cols = get_arg(args, ['c','cols','columns'])

  if not cols or isinstance(cols,list):
    return cols

  cols = cols.split(',')
  columns = []
  for c in cols:
    c = c.strip()
    if '-' in c:
      c = tuple(tryint(f.strip()) for f in c.split('-'))
      if len(c) != 2 or '' in c:
        raise ValueError('Invalid range specified: %s' % str(c))
    else:
      c = tryint(c)
    columns.append(c)

  return columns


def headers_needed(columns):
  '''
  Return whether headers are needed

  @param      columns: indices, names, or ranges of columns to select, comma
                       delimited
  @type       columns: list of strings, integers, or 2-tuples for ranges
  @return            : whether headers are needed
  @rtype             : bool
  '''
  for col in columns:
    if isinstance(col,tuple):
      if isstr(col[0]) or isstr(col[1]):
        return True
    elif isstr(col):
      return True
  return False


def resolve_column_headers(header,columns):
  '''
  @param      columns: indices, names, or ranges of columns to select, comma
                       delimited
  @type       columns: list of strings, integers, or 2-tuples for ranges
  @param       header: header line of the input file
  @type        header: sequence of strs
  @return            : resolved header line
  @rtype             : sequence of strs
  '''

  for i,col in enumerate(columns):
    if isstr(col):
      try:
        yield header.index(col)
      except ValueError:
        raise ValueError("Cannot find header '%s'" % col)
    elif isinstance(col,tuple):
      assert len(col) == 2

      if isstr(col[0]):
        try:
          c1 = header.index(col[0])
        except ValueError:
          raise ValueError("Cannot find header '%s'" % col[0])
      else:
        c1 = col[0]

      if isstr(col[1]):
        try:
          c2 = header.index(col[1])
        except ValueError:
          raise ValueError("Cannot find header '%s'" % col[1])
      else:
        c2 = col[1]

      yield (c1,c2)
    else:
      yield col


def load_table(filename,want_header=False,extra_args=None,**kwargs):
  '''
  Return a table of data read from a delimited file as a list of rows, based
  on several parameters.

  Valid files should be composed of standard ASCII lines of text with tab or
  other delimited fields.  If 'columns' is not specified, then all fields
  are considered, otherwise only certain columns will be returned.
  Otherwise, columns must contain a list of column indices, strings
  representing the name of columns to be selected, or two-tuples of either
  strings or integers representing column range selections.  Whitespace is
  stripped from the beginning and end of every field.  If the skip parameter
  is used to ignore a certain number of lines (e.g., headers) at the
  beginning of the file.

  Results are returned as a sequence of rows, each row a list of field
  values.  Empty strings will be returned for rows with too few columns, so
  that all row results will have an identical number of fields.  If no
  columns are explicitly requested, the resulting rows will be return as-is.

  The dialect argument can be used to specify a Python csv module dialect
  name or Dialect object.  In addition, the following Python csv module
  options can be specified by appending information to the specified
  filename: dialect, delimiter, doublequote, escapechar, lineterminator,
  quotechar, and quoting.  The syntax is ':option=value', appended to the
  end of the filename.

  A skip parameter can be specified is used to ignore a certain number of
  lines (e.g., headers) at the beginning of the file.  It can be specified
  as either a function argument or by appending ':skip=n', where n is an
  integer, to the filename.

  The columns parameter can be overridden with a comma seperated list of
  column indices or column ranges (0 based) by appending, e.g.,
  ':columns=1,0,5-10' to the end of the filename.  Simple parsing of column
  headers is supported.  If any non-integer columns are specified, then the
  first non-skipped row of the file is assumed to contain heading names and
  indices are determined by finding the appropriate headings.  For example,
  ':columns=0,1,PID,SID,0' on a file with the header 'PID,SID' is equivalent
  to specifying ':skip=1:columns=0,1,0,1,0'.  No escaping or quoting is
  allowed.

  The following parameters and aliases are accepted as part of the augmented
  filename:
        skip, s: number of header lines to skip
     columns, c: indices, names, or ranges of columns to select, comma delimited
        dialect: csv module dialect name ('csv' or 'tsv')
      delimiter: single field delimiter character
    doublequote: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
    doublequote: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
     escapechar: one-character string that removes any special meaning
                 from the following character.
      quotechar: one-character string used to quote fields containing
                 special characters, such as the delimiter or quotechar, or
                 which contain new-line characters.
        quoting: Controls when quotes characters are recognised.  It can
                 take on any of the cav QUOTE_ constants (without the QUOTE_
                 prefix).

  @param     filename: file name or file object
  @type      filename: str or file object
  @param  want_header: flag to request header to be returned (if present)
  @type   want_header: bool
  @param         skip: number of initial lines to be skipped (defaults to 0)
  @type          skip: int
  @param      columns: indices, names, or ranges of columns to select, comma
                       delimited
  @type       columns: list of strings, integers, or 2-tuples for ranges
  @param      dialect: csv module dialect name or dialect object
  @type       dialect: str or csv.Dialect
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @return            : sequence of rows containing the columns requested
  @rtype             : generator

  >>> def test(f,**kw): return list(load_table(f,**kw))
  >>> from StringIO import StringIO
  >>> test(StringIO("loc1\\tlocA \\nloc2\\n\\n\\t\\n loc3 \\tlocC\\n\\tfoo\\nloc4\\t\\n"))
  [['loc1', 'locA'], ['loc2', ''], ['', ''], ['', ''], ['loc3', 'locC'], ['', 'foo'], ['loc4', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),skip=2)
  [['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[0,1])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[1,0])
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c1','c2'])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[('c1','c2')])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c2','c1'])
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1\\tH2\\nv11\\tv12\\nv21\\n\\tv32')
  >>> f.flush()
  >>> test(f.name + ':skip=1')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':skip=1:columns=0,1')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':skip=1:columns=0-1')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1,H2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1-H2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1-H2,H2,0,0-1')
  [['v11', 'v12', 'v12', 'v11', 'v11', 'v12'], ['v21', '', '', 'v21', 'v21', ''], ['', 'v32', 'v32', '', '', 'v32']]
  >>> del f
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  dialect = get_csv_dialect(args,'tsv')
  skip    = int(get_arg(args, ['s','skip'], 0))
  columns = get_arg_columns(args)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name)
  lines = csv.reader(lfile, **dialect)

  if skip:
    lines = islice(lines,skip,None)

  # All columns are to be returned
  if not columns:
    try:
      row = lines.next()
    except StopIteration:
      return

    row = map(str.strip,row)
    yield row

    n = len(row)
    for row in lines:
      result  = map(str.strip,row)
      result += ['']*(n-len(result))
      yield result

    return

  # Otherwise, resolve column heading names into indices
  header = None
  if headers_needed(columns):
    try:
      header = map(str.strip,lines.next())
    except StopIteration:
      return

    columns = list(resolve_column_headers(header,columns))

  # Build indices
  indices = []
  for col in columns:
    if isinstance(col,tuple):
      indices.extend( xrange(col[0],col[1]+1) )
    else:
      indices.append(col)

  if want_header and header is not None:
    yield [ header[j] for j in indices ]

  n = len(indices)

  # Build result rows
  for row in lines:
    m = len(row)
    result = [ (row[j].strip() if j<m else '') for j in indices ]
    if header is not None:
      result += ['']*(n-len(result))
    yield result


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
