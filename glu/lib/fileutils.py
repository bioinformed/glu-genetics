# -*- coding: utf-8 -*-

__abstract__  = 'file related utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import csv

from   operator      import itemgetter
from   itertools     import islice,izip,chain

from   glu.lib.utils import is_str,as_set,peekfirst,deprecated_by,unique,tally


__all__ = ['autofile','namefile','hyphen',
           'guess_format','related_file','guess_related_file',
           'list_reader', 'map_reader', 'table_reader', 'table_writer',
           'cook_table', 'sort_table', 'uniq_table', 'subset_variables',
           'create_categorical_variables']


TABLE_FORMATS = set(['xls','csv'])
COMPRESSED_SUFFIXES = {'gz':'gzip', 'Z'  :'gzip',
                       'bz':'bzip2','bz2':'bzip2'}


# Create more standard aliases for Python CSV module dialects.  'excel' is
# now available as 'csv' and 'excel-tab' is now available as 'tsv'.
csv.register_dialect('csv', csv.get_dialect('excel'))
csv.register_dialect('tsv', csv.get_dialect('excel-tab'))


def spawn_compressor(binary, filename, mode, bufsize=-1):
  '''
  Spawn a subprocess to run a compressor like gzip or bzip2 and connect to
  input/output pipes

  @param    binary: executable name
  @type     binary: str
  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @param   bufsize: buffering mode and size.  0=unbuffered, 1=linebuffered, >1 buffer size,
                    -1 default buffering (default)
  @type    bufsize: int
  @return         : file object to read from or write to
  @rtype          : file object
  '''
  if os.environ.get('GLU_NOSPAWN'):
    raise OSError('Spawning external processes disabled by GLU_NOSPAWN')

  from subprocess import Popen,PIPE

  if 'w' in mode:
    out = file(filename,mode)
    cmd = [binary,'-c']
    f   = Popen(cmd, stdin=PIPE, stdout=out, bufsize=bufsize).stdin
  else:
    cmd = [binary,'-d','-c',filename]
    f   = Popen(cmd, stdout=PIPE, universal_newlines='U' in mode, bufsize=bufsize).stdout

  return f


def hyphen(filename,defaultfile,args=None):
  '''
  Return the default if input is '-', otherwise itself

  @param  defaultfile: default file name or file object
  @type   defaultfile: str or file object
  @param     filename: file name or file object
  @type      filename: str or file object
  @param         args: augmented argument output dictionary or None
  @type          args: dict or None
  @return            : file name or file object to read from or write to
  @rtype             : str or file object
  '''
  if not is_str(filename):
    return filename

  if args is None:
    args = {}

  if parse_augmented_filename(filename,args) == '-':
    return defaultfile

  return filename


def namefile(filething):
  '''
  Return a human-comprehensible name for the file or file object provided.
  Recognizes file objects with the 'name' attribute, including sys.stdin,
  sys.stdout, and sys.stderr.

  @param  filething: file or filename
  @type   filething: file object or string
  @return:           human-comprehensible file name

  >>> import sys
  >>> namefile(file('/dev/null'))
  '/dev/null'
  >>> namefile(sys.stdin)
  '<stdin>'
  >>> namefile(sys.stderr)
  '<stderr>'
  >>> namefile('/dev/null')
  '/dev/null'
  '''
  if is_str(filething):
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
  @return         : compression format, if compressed, otherwise an empty string
  @rtype          : str

  >>> compressed_filename('subjects.sdat')
  ''
  >>> compressed_filename('subjects.sdat.gz')
  'gzip'
  >>> compressed_filename('../subjects.sdat')
  ''
  >>> compressed_filename('../subjects.sdat.gz')
  'gzip'
  >>> compressed_filename('../subjects.sdat.bz2')
  'bzip2'
  '''
  if not is_str(filename):
    return ''

  filename = os.path.expanduser(filename)
  parts    = os.path.basename(filename).split('.')
  ext      = parts[-1] if parts else ''
  return COMPRESSED_SUFFIXES.get(ext,'')


def autofile(filename, mode='r'):
  '''
  Return a file object in the correct compressed format as specified, which
  is ready to read from or write to

  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @return         : file object to read from or write to
  @rtype          : file object
  '''
  # Pass non-string filename objects back as file-objects
  if not is_str(filename):
    return filename

  filename = os.path.expanduser(filename)
  comp     = compressed_filename(filename)

  if not comp:
    f = file(filename, mode)
  elif comp == 'gzip':
    try:
      f = spawn_compressor(os.environ.get('GLU_GZIP','gzip'), filename, mode)
    except _autofile_errors:
      import gzip
      f = gzip.GzipFile(filename, mode)
  elif comp == 'bzip2':
    try:
      f = spawn_compressor(os.environ.get('GLU_BZIP2','bzip2'), filename, mode)
    except _autofile_errors:
      import bz2
      f = bz2.BZ2File(filename, mode)

  return f


def guess_format(filename, formats, args=None):
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

  if not is_str(filename):
    return None

  # Parse to remove augmented arguments
  if args is None:
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
  if not is_str(filename):
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

  >>> guess_related_file('fileutils.dat',['py']) # doctest:+SKIP
  'fileutils.py'
  '''
  if not is_str(filename):
    return None

  prefix,ext = os.path.splitext(filename)

  if not prefix:
    raise ValueError('invalid filename')

  for new_ext in extensions:
    testfile = '%s.%s' % (prefix,new_ext)
    if os.path.isfile(testfile):
      return testfile

  return None


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


DIALECT_KWARGS = ['dialect','delimiter','doublequote','escapechar','lineterminator',
                  'quotechar','quoting','skipinitialspace','strict']


_unescape_literal_chars = [('\\r','\r'),('\\t','\t'),('\\n','\n'),('\\\\','\\')]


def unescape_literal(s):
  '''
  Evaluate C/Python character literals \r, \t, \n, and \\ into the
  respective single-character ASCII representation
  '''
  for e,c in _unescape_literal_chars:
    s = s.replace(e,c)
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
      dargs[darg] = unescape_literal(args.pop(darg))

  if 'quoting' in dargs:
    q = 'QUOTE_' + dargs['quoting'].upper()
    dargs['quoting'] = getattr(csv, q, q)

  # Parse truth values
  for tf in ('skipinitialspace','doublequote'):
    if tf in dargs:
      dargs[tf] = trybool(dargs[tf])

  return dargs


def _literal_list(filename):
  return is_str(filename) and not os.path.exists(filename) and filename.startswith(':')


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


def list_reader(filename,extra_args=None,**kwargs):
  '''
  Load list of values from a file or literal list.

  By default, a file is interpreted as a tab-delimited file with only the
  first field considered.  Blank lines are skipped, whitespace is stripped
  from the beginning and end of each field selected.  These settings can be
  modified via function arguments or by appending additional options to the
  filename.  File loading is delegated to the table_reader function, so refer
  to table_reader for details on supported arguments.

  The column number or name may be specified as per the table_reader function.
  However, for backward compatibility the 'index' parameter is also
  supported as a synonym to 'columns'.

  If the filename is a string that begins with ':', it is interpreted as a
  literal list of comma separated values, ignoring all other parameters
  (including skip!).  No escaping or quoting is allowed.

  The following parameters and aliases are accepted as part of the augmented
  filename (all but the first from table_reader):

  [by list_reader]
       index, i: column index or header name

  [by table_reader]
     columns, c: indices, names, or ranges of columns to select, comma delimited
       skip,  s: number of header lines to skip
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
                 take on any of the csv QUOTE_ constants (without the QUOTE_
                 prefix).

  @param         filename: file name, file object, or literal string
  @type          filename: str or file object
  @param           hyphen: if the input filename is '-', then use this file
                           object instead
  @type            hyphen: file object
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
  >>> list_reader(StringIO("loc1\\tsample1\\nloc2\\n\\n loc3\\n"))
  ['loc1', 'loc2', 'loc3']
  >>> list_reader(StringIO("loc1,sample1\\nloc2\\n\\n loc3\\n"),dialect='csv')
  ['loc1', 'loc2', 'loc3']
  >>> list_reader(StringIO("loc1\\tsample1\\nloc2\\n\\n loc3\\n"),index=1)
  ['sample1']
  >>> list_reader(':loc1, loc2, loc3')
  ['loc1', 'loc2', 'loc3']
  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1\\tH2\\nv11\\tv12\\nv21\\n\\tv32')
  >>> f.flush()
  >>> list_reader(f.name + ':skip=1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':skip=1:index=1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':delimiter=\\t:index=H1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':dialect=excel-tab:skip=1:index=2')
  ['v12', 'v32']
  >>> list_reader(f.name + ':index=H2')
  ['v12', 'v32']
  >>> list_reader(f.name + ':c=H2')
  ['v12', 'v32']
  >>> list_reader(f.name + ':c=2:skip=1')
  ['v12', 'v32']
  >>> list_reader(f.name + ':i=2:c=1:skip=1')
  Traceback (most recent call last):
     ...
  ValueError: Invalid specification of both index and columns
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1,H2\\nv11,v12\\nv21\\n,v32')
  >>> f.flush()
  >>> list_reader(f.name + ':delimiter=,:skip=1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':delimiter=,:skip=1:index=1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':delimiter=,:index=H1')
  ['v11', 'v21']
  >>> list_reader(f.name + ':dialect=excel:skip=1:index=2')
  ['v12', 'v32']
  >>> list_reader(f.name + ':dialect=excel:index=H2')
  ['v12', 'v32']
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  # Process literal syntax
  if _literal_list(filename):
    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    return [ intern(l.strip()) for l in filename[1:].split(',') if l ]

  # Otherwise, handle file name or file object cases
  name    = parse_augmented_filename(filename,args)
  columns = get_arg(args, ['c','cols','columns'])
  index   = tryint1(get_arg(args, ['index','i']))

  if index is not None and columns is not None:
    raise ValueError('Invalid specification of both index and columns')

  if index is None and columns is None:
    index = 0
  elif index is None:
    index = columns

  items = table_reader(name, columns=index, extra_args=args)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  return [ intern(i[0]) for i in items if i and i[0] ]


# Anonymous object
_nothing = object()


def map_reader(filename,unique=True,extra_args=None,**kwargs):
  '''
  Creates a dictionary representing a list or mapping from a text file.

  By default, valid files should be composed of standard ASCII lines of text
  with tab delimited fields.  Only the first and second fields of each line
  are considered. Whitespace is stripped from the beginning and end of every
  field considered.  If the skip parameter is used to ignore a certain
  number of lines (e.g., headers) at the beginning of the file.  A default
  parameter may be specified to assign values to keys with empty or
  non-existent value fields.  Otherwise, the value will be set equal to the
  key.  File loading is delegated to the table_reader function, so refer
  to table_reader for details on supported arguments.

  The key and value column number or name may be specified as per the
  table_reader function using the 'columns' argument.  The first and second
  columns returned are taken as keys and values, respectively.  This is the
  recommended mode of operation.  For backward compatibility, the the index
  or name of key and value columns can be overridden with the key_index and
  value_index parameters (0 based indices or column names) or by appending
  ':key=i' and/or ':value=j' to the end of the filename (1 based indices or
  column names).  No escaping or quoting is allowed.

  #FIXME: Document literal syntax

  The default behavior (unique=True) requires that all mappings are unique
  and consistent in the sense that every key must map to exactly one value,
  though that mapping may appear redundantly in the input.  Otherwise a
  ValueError is raised.  If unique=False, then a dictionary of each key to a
  list of unique values is returned.

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

  [by map_reader]
      key_index, i: index of field to select or name of header for keys
    value_index, v: index of field to select or name of header for values
     default,def,d: default value for keys with no or empty value

  [by table_reader]
        columns, c: indices, names, or ranges of columns to select, comma delimited
          skip,  s: number of header lines to skip
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
                    take on any of the csv QUOTE_ constants (without the QUOTE_
                    prefix).

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       hyphen: if the input filename is '-', then use this file
                       object instead
  @type        hyphen: file object
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

  >>> def test(f,**kw): return sorted(map_reader(f,**kw).iteritems())
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
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[0,1])
  [('loc1', 'loc1'), ('loc2', 'loc2')]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[1,0])
  [('loc1', 'loc1')]
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
  >>> test(f.name + ':skip=1:key=1:value=2')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':skip=1:key=1')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':key=H1:v=H2')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':value=H2')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':k=H1')
  [('v11', 'v12'), ('v21', 'v21')]
  >>> test(f.name + ':skip=1:key=2:value=1')
  [('v12', 'v11'), ('v32', 'v32')]
  >>> test(f.name + ':key=H2:value=H1')
  [('v12', 'v11'), ('v32', 'v32')]
  >>> test(f.name + ':c=H2,H1')
  [('v12', 'v11'), ('v32', 'v32')]
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  # Handle :literal syntax
  if _literal_list(filename):
    rows = (kv.split('=') for kv in filename[1:].split(',') if kv)
    default = get_arg(args, ['default'], _nothing)

  else:
    name        = parse_augmented_filename(filename,args)
    default     = get_arg(args, ['d','def','default'], _nothing)
    key_index   = tryint1(get_arg(args, ['k','key',  'key_index'  ]))
    value_index = tryint1(get_arg(args, ['v','value','value_index']))
    columns     = get_arg(args, ['c','cols','columns'])

    if (key_index is not None or value_index is not None) and columns is not None:
      raise ValueError('Invalid specification of both key/value and columns')

    if key_index is None:
      key_index = 0

    if value_index is None:
      value_index = 1

    if columns is None:
      columns = [key_index,value_index]

    rows = table_reader(name,columns=columns,extra_args=args)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  # Parse the data file
  def _map_reader():
    for row in rows:
      n = len(row)

      key = intern(row[0].strip()) if n else None

      if not key:
        continue

      value = intern(row[1].strip()) if n>1 else None

      if not value:
        if default is not _nothing:
          value = default
        else:
          value = key

      yield key,value

  m = {}
  unique_fails = []
  for key,value in _map_reader():
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


def resolve_column_header_atom(header,column):
  '''
  Resolve an atomic column header into a column index.  Atoms can be header
  names, string indices (1-based), or integer indices (0-based). Raises
  ValueError if the column index could not be resolved.

  @param header: optional schema header
  @type  header: list or None
  @param column: column name, index string (1-based), or index number (0-based)
  @type  column: int or str
  @return      : column index
  @rtype       : int

  >>> resolve_column_header_atom(['H1','H2'],0)
  0
  >>> resolve_column_header_atom(['H1','H2'],1)
  1
  >>> resolve_column_header_atom(['H1','H2'],'1')
  0
  >>> resolve_column_header_atom(['H1','H2'],'2')
  1
  >>> resolve_column_header_atom(['H1','H2'],'H1')
  0
  >>> resolve_column_header_atom(['H1','H2'],'H2')
  1
  '''
  if header is not None:
    try:
      return header.index(column)
    except ValueError:
      pass

  col = tryint1(column)

  if isinstance(col,int):
    if col < 0:
      raise ValueError('Invalid negative column index')

    return col

  raise ValueError("Cannot find header '%s'" % (column,))


def resolve_column_header(header,column):
  '''
  Resolve a complex column header into a column index.  Complex columns can
  be atoms, range tuples, or range strings.  Range tuples are 2-tuples of
  atoms.  Range strings are hyphen delimited atoms.  Raises ValueError if
  the column index could not be resolved.

  @param header: optional schema header
  @type  header: list or None
  @param column: atom, range tuple, or range string
  @type  column: int, str, 2-tuple of atoms
  @return      : column index
  @rtype       : int

  >>> resolve_column_header(['H1','H2'],0)
  0
  >>> resolve_column_header(['H1','H2'],1)
  1
  >>> resolve_column_header(['H1','H2'],'1')
  0
  >>> resolve_column_header(['H1','H2'],'2')
  1
  >>> resolve_column_header(['H1','H2'],'H1')
  0
  >>> resolve_column_header(['H1','H2'],'H2')
  1
  >>> resolve_column_header(['H1','H2'],'1-1')
  (0, 0)
  >>> resolve_column_header(['H1','H2'],'1-2')
  (0, 1)
  >>> resolve_column_header(['H1','H2'],'H1-H1')
  (0, 0)
  >>> resolve_column_header(['H1','H2'],'H1-H2')
  (0, 1)
  >>> resolve_column_header(['H1','H2'],'-H1')
  (0, 0)
  >>> resolve_column_header(['H1','H2'],'H2-')
  (1, 1)
  '''
  try:
    return resolve_column_header_atom(header,column)
  except ValueError:
    pass

  if isinstance(column,tuple):
    assert len(column)==2
    try:
      r1 = resolve_column_header_atom(header,column[0])
      r2 = resolve_column_header_atom(header,column[1])
      return r1,r2
    except ValueError:
      pass

  if is_str(column):
    hyphens = column.count('-')
    if hyphens:
      start = 0
      for i in range(hyphens):
        r = column.find('-',start)
        r1,r2 = column[:r].strip(),column[r+1:].strip()
        start = r+1

        if not r1:
          r1 = 0
        if not r2:
          if header is None:
            raise ValueError('Unbounded ranges require headers to be specified')
          r2 = max(0,len(header)-1)

        try:
          r1 = resolve_column_header_atom(header,r1)
          r2 = resolve_column_header_atom(header,r2)
          return min(r1,r2),max(r1,r2)
          break
        except ValueError:
          pass

  # Reraise last error
  raise


def resolve_column_headers(header,include,exclude=None):
  '''
  @param       header: header line of the input file
  @type        header: sequence of strs
  @param      include: indices, names, or ranges of columns to select, comma
                       delimited
  @type       include: list of strings, integers, or 2-tuples for ranges
  @param      exclude: indices, names, or ranges of columns to exclude, comma
                       delimited
  @type       exclude: list of strings, integers, or 2-tuples for ranges
  @return            : resolved header line
  @rtype             : sequence of strs
  '''
  if isinstance(include,int):
    include = [include]
  elif isinstance(include,str):
    include = include.split(',')

  if include is None:
    if header:
      indices = range(len(header))
    else:
      raise ValueError('Cannot resolve indices for unknown header')
  else:
    indices = []
    for column in include:
      col = resolve_column_header(header,column)
      if isinstance(col,tuple):
        indices.extend( xrange(col[0],col[1]+1) )
      else:
        indices.append(col)

  if exclude:
    exclude_indices = set(resolve_column_headers(header,exclude))
    indices = [ i for i in indices if i not in exclude_indices ]

  return indices


def table_columns(rows, columns, drop=None, header=None, want_header=False):
  '''
  Return a rows from a sequence for columns specified by name or number. A
  single row of column headers may be embedded within the row sequence or
  provided out-of-band.  Similarly, the header may be returned in-band or
  not at all.

  This function is typically used to post-process rows of data read from
  various delimited formats to select the number and order of columns.  In
  addition, all rows are padded to be of length of at least the header or
  first row, though longer rows are never truncated.  This simplifies many
  downstream operations.

  The columns parameter can be None, a comma separated list of column
  indices or column names.  None implies that all columns should be
  returned.  Index numbers given as strings are interpreted as 1-based
  indexes, while indices provided as Python objects are treated as 0-based.
  The rationale behind the different bases is that users will typically
  specify the string-based indices and counting from 1 is much more natural.
  Programmatic interfaces will provide native integers and will follow the
  standard 0-based convention.

  @param         rows: input row sequence
  @type          rows: iterable sequence of sequences
  @param      columns: columns to extract
  @type       columns: string or sequence
  @param       header: out-of-band specification of row headers
  @type        header: sequence
  @param  want_header: flag to request header to be returned as part of the
                       output stream
  @type   want_header: bool
  @return            : sequence of rows containing the columns requested
  @rtype             : generator

  >>> list(table_columns([['loc1','locA '],['loc2'],[''],['',''],[' loc3 ','locC'],['','foo'],['loc4','']],None))
  [['loc1', 'locA'], ['loc2', ''], ['', ''], ['', ''], ['loc3', 'locC'], ['', 'foo'], ['loc4', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],[0,1]))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],[1,0]))
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],['c1','c2']))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],None,drop='c2'))
  [['loc1'], ['loc2'], ['loc1'], ['loc2']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],None,drop=[1]))
  [['loc1'], ['loc2'], ['loc1'], ['loc2']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],[('c1','c2')]))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],[(0,'c2')]))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],[('1','c2')]))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['c1','c2'],['loc1'],['loc2'],['loc1','loc1'],['loc2']],['c2','c1']))
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],['c1','c2'],header=['c1','c2']))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],[('c1','c2')],header=['c1','c2']))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],[(0,'c2')],header=['c1','c2'],want_header=True))
  [['c1', 'c2'], ['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],[('1','c2')],header=['c1','c2']))
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> list(table_columns([['loc1'],['loc2'],['loc1','loc1'],['loc2']],['c2','c1'],header=['c1','c2']))
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  '''
  rows = iter(rows)

  # All columns are to be returned
  if not columns and not drop:
    def _table_reader_all(header):
      if header is None:
        try:
          row = rows.next()
        except StopIteration:
          return

        # Process row 1 (may or may not be a header)
        row = map(str.strip,row)
        n   = len(row)
        yield row

      else:
        header = map(str.strip,header)
        n      = len(header)
        if want_header:
          yield header

      for row in rows:
        result  = map(str.strip,row)
        result += ['']*(n-len(result))
        yield result

    return _table_reader_all(header)

  # Otherwise, resolve column heading names into indices
  try:
    indices = resolve_column_headers(header,columns,drop)
  except ValueError:
    if header is not None:
      raise

    try:
      header = map(str.strip,rows.next())
    except StopIteration:
      return

    indices = resolve_column_headers(header,columns,drop)

  if not want_header:
    header = None

  def _table_reader_columns():
    if header is not None:
      yield [ header[j] for j in indices ]

    # Build result rows
    for row in rows:
      m = len(row)
      result  = [ (row[j].strip() if j<m else '') for j in indices ]
      yield result

  return _table_reader_columns()


def table_reader(filename,want_header=False,extra_args=None,**kwargs):
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

  The columns parameter can be overridden with a comma separated list of
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
     columns, c: indices, names, or ranges of columns to select, comma delimited
        skip, s: number of header lines to skip
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
                 take on any of the csv QUOTE_ constants (without the QUOTE_
                 prefix).

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       hyphen: if the input filename is '-', then use this file
                       object instead
  @type        hyphen: file object
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

  >>> def test(f,**kw): return list(table_reader(f,**kw))
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
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[(0,'c2')])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[('1','c2')])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c2','c1'])
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c1','c2'],header=['c1','c2'])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[('c1','c2')],header=['c1','c2'])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[(0,'c2')],header=['c1','c2'],want_header=True)
  [['c1', 'c2'], ['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[('1','c2')],header=['c1','c2'])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("loc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c2','c1'],header=['c1','c2'])
  [['', 'loc1'], ['', 'loc2'], ['loc1', 'loc1'], ['', 'loc2']]

  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1\\tH2\\nv11\\tv12\\nv21\\n\\tv32')
  >>> f.flush()
  >>> test(f.name + ':skip=1')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':skip=1:columns=1,2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':skip=1:columns=1-2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1,H2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1-H2')
  [['v11', 'v12'], ['v21', ''], ['', 'v32']]
  >>> test(f.name + ':columns=H1-H2,H2,1,1-2')
  [['v11', 'v12', 'v12', 'v11', 'v11', 'v12'], ['v21', '', '', 'v21', 'v21', ''], ['', 'v32', 'v32', '', '', 'v32']]
  >>> del f
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  skip    = int(get_arg(args, ['skip','s'], 0))
  format  = get_arg(args, ['format'], guess_format(name, TABLE_FORMATS)) or 'tsv'
  columns = get_arg(args, ['c','cols','columns'])
  drop    = get_arg(args, ['d','drop'])
  header  = get_arg(args, ['h','header'])

  if header is not None and isinstance(header,str):
    header = map(str.strip,header.split(','))

  format = format.lower()
  if format in ('xls','excel'):
    rows = table_reader_excel(name, extra_args=args)
  elif format in ('delimited','tsv','csv'):
    rows = table_reader_delimited(name, extra_args=args)
  else:
    raise NotImplementedError("File format '%s' is not supported" % format)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if skip:
    rows = islice(rows,skip,None)

  return table_columns(rows,columns,drop,header=header,want_header=want_header)


class TableWriter(object):
  '''
  Write selected columns to a lower-level tabular data writer object
  '''
  def __init__(self, writer, columns, drop):
    self.writer  = writer
    self.columns = columns
    self.drop    = drop

    # Try to resolve column headings without a header row
    try:
      self.indices = resolve_column_headers(None,columns,drop)
    except ValueError:
      self.indices = None

  def writerows(self, rows):
    indices = self.indices

    if indices is None:
      # Infer indices from the first row if needed
      try:
        row,rows = peekfirst(rows)
      except IndexError:
        return

      self.indices = indices = resolve_column_headers(row,self.columns,self.drop)

    for row in rows:
      m = len(row)
      row = [ (row[j] if j<m else '') for j in indices ]
      self.writer.writerow(row)

  def writerow(self, row):
    indices = self.indices

    # First row and indices require header information
    if indices is None:
      self.indices = indices = resolve_column_headers(row,self.columns,self.drop)

    m = len(row)
    row = [ (row[j] if j<m else '') for j in indices ]
    self.writer.writerow(row)


##########################################################################
# Preserve deprecated load_* APIs until just before 1.0 release

@deprecated_by('list_reader')
def load_list(*args, **kwargs):
  return list_reader(*args, **kwargs)
load_list.__doc__ = list_reader.__doc__

@deprecated_by('map_reader')
def load_map(*args, **kwargs):
  return map_reader(*args, **kwargs)
load_map.__doc__ = map_reader.__doc__

@deprecated_by('table_columns')
def load_table_rows(*args, **kwargs):
  return table_columns(*args, **kwargs)
load_table_rows.__doc__ = table_columns.__doc__

@deprecated_by('table_reader')
def load_table(*args, **kwargs):
  return table_reader(*args, **kwargs)
load_table.__doc__ = table_reader.__doc__

##########################################################################


def table_writer(filename,extra_args=None,**kwargs):
  '''
  Return an object that can write a table of data to a delimited file based
  on several parameters.

  Output will be composed of standard ASCII lines of text with tab or
  other delimited fields.

  The dialect argument can be used to specify a Python csv module dialect
  name or Dialect object.  In addition, the following Python csv module
  options can be specified by appending information to the specified
  filename: dialect, delimiter, doublequote, escapechar, lineterminator,
  quotechar, and quoting.  The syntax is ':option=value', appended to the
  end of the filename.

  The following parameters and aliases are accepted as part of the augmented
  filename:
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
                 take on any of the csv QUOTE_ constants (without the QUOTE_
                 prefix).

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      dialect: csv module dialect name or dialect object
  @type       dialect: str or csv.Dialect
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @return            : sequence of rows containing the columns requested
  @rtype             : generator

  >>> from StringIO import StringIO
  >>> o=StringIO()
  >>> w=table_writer(o,dialect='csv')
  >>> w.writerow(['1','2','3'])
  >>> w.writerow(['a','b','c'])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  1,2,3
  a,b,c

  >>> o=StringIO()
  >>> w=table_writer(o,delimiter='|')
  >>> w.writerow(['1','2','3'])
  >>> w.writerow(['a','b','c'])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  1|2|3
  a|b|c

  >>> o=StringIO()
  >>> w=table_writer(o,delimiter='|',columns=0)
  >>> w.writerow(['H1','H2','H3'])
  >>> w.writerow(['1','2','3'])
  >>> w.writerow(['a','b','c'])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  H1
  1
  a

  >>> o=StringIO()
  >>> w=table_writer(o,columns='1-2')
  >>> w.writerow(['H1','H2','H3'])
  >>> w.writerow(['1','2','3'])
  >>> w.writerow(['a','b','c'])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  H1	H2
  1	2
  a	b

  >>> o=StringIO()
  >>> w=table_writer(o,columns='2-H3')
  >>> w.writerows([['H1','H2','H3'],
  ...              ['1','2','3'],
  ...              ['a','b','c']])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  H2  H3
  2   3
  b   c

  >>> o=StringIO()
  >>> w=table_writer(o,drop='H2')
  >>> w.writerows([['H1','H2','H3'],
  ...              ['1','2','3'],
  ...              ['a','b','c']])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  H1  H3
  1   3
  a   c

  >>> o=StringIO()
  >>> w=table_writer(o,columns='1-3',drop='H2-H3')
  >>> w.writerows([['H1','H2','H3'],
  ...              ['1','2','3'],
  ...              ['a','b','c']])
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  H1
  1
  a
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  format  = get_arg(args, ['format'], guess_format(name, TABLE_FORMATS)) or 'tsv'
  columns = get_arg(args, ['c','cols','columns'])
  drop    = get_arg(args, ['d','drop'])

  format  = format.lower()
  if format in ('xls','excel'):
    writer = ExcelWriter(name, extra_args=args)
  elif format in ('delimited','tsv','csv'):
    writer = delimited_table_writer(name, extra_args=args)
  else:
    raise NotImplementedError("File format '%s' is not supported" % format)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  # Use a wrapper writer object to handle column selectors
  if columns is not None or drop:
    writer = TableWriter(writer, columns, drop)

  return writer


def table_reader_delimited(filename, extra_args=None, **kwargs):
  '''
  Return a configured delimited table reader
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  dialect = get_csv_dialect(args, guess_format(name, ['csv']) or 'tsv')
  hyin    = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name) if name!='-' or hyin is None else hyin

  return csv.reader(lfile,**dialect)


def delimited_table_writer(filename, extra_args=None, **kwargs):
  '''
  Return a configured delimited table writer
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  dialect = get_csv_dialect(args, guess_format(name, ['csv']) or 'tsv')
  hyout   = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name,'wb') if name!='-' or hyout is None else hyout

  return csv.writer(lfile,**dialect)


XLS_NULLDATE = (0,0,0)

def _xlate_xls_row_object(book,values,types):
  '''
  Translate a sequence of native Excel values and types into Python objects
  '''
  import xlrd
  from   datetime  import date,time,datetime

  row = []
  t = set(types)
  t -= set([xlrd.XL_CELL_NUMBER,xlrd.XL_CELL_DATE,xlrd.XL_CELL_BOOLEAN,xlrd.XL_CELL_ERROR])
  for value,typ in izip(values,types):
    if typ == xlrd.XL_CELL_NUMBER:
      ivalue = int(value)
      if ivalue == value:
        value = ivalue
    elif typ == xlrd.XL_CELL_DATE:
      if value[:3] == XLS_NULLDATE:
        value = time(*value[3:])
      elif value[3:] == XLS_NULLDATE:
        value = date(*value[:3])
      else:
        value = datetime(*value)
    elif typ == xlrd.XL_CELL_BOOLEAN:
      value = bool(value)
    elif typ == xlrd.XL_CELL_ERROR:
      value = xlrd.error_text_from_code[value]
    elif isinstance(value,unicode):
      value = value.encode('utf8')

    row.append(value)

  return row


def _xlate_xls_row_str(book,values,types):
  '''
  Translate a sequence of native Excel values and types into strings
  '''
  import xlrd

  row = []
  t = set(types)
  t -= set([xlrd.XL_CELL_NUMBER,xlrd.XL_CELL_DATE,xlrd.XL_CELL_BOOLEAN,xlrd.XL_CELL_ERROR])
  for value,typ in izip(values,types):
    if typ == xlrd.XL_CELL_NUMBER:
      ivalue = int(value)
      if ivalue == value:
        value = '%d' % value
      else:
        value = '%f' % value
    elif typ == xlrd.XL_CELL_DATE:
      value = xlrd.xldate_as_tuple(value,book.datemode)
      if value[:3] == XLS_NULLDATE:
        value = '%02d:%02d:%02d' % value[3:]
      elif value[3:] == XLS_NULLDATE:
        value = '%04d/%02d/%02d' % value[:3]
      else:
        value = '%04d/%02d/%02d %02d:%02d:%02d' % value
    elif typ == xlrd.XL_CELL_BOOLEAN:
      value = ['0','1'][value]
    elif typ == xlrd.XL_CELL_ERROR:
      value = xlrd.error_text_from_code[value]
    elif isinstance(value,unicode):
      value = value.encode('utf8')

    row.append(value)

  return row


def table_reader_excel(filename,strdata=True,extra_args=None,**kwargs):
  '''
  Load rows from a Microsoft Excel (XLS) file using the xlrd module

  Supports Excel versions: 2003, 2002, XP, 2000, 97, 95, 5.0, 4.0, 3.0, but not
  Excel 2007 XML (XLSX).
  '''
  try:
    import xlrd
  except ImportError:
    raise ValueError('Missing xlrd module to read Microsoft Excel file')

  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name  = parse_augmented_filename(filename,args)
  sheet = get_arg(args, ['sheet'],0)
  hyin  = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if hyin is not None and name=='-':
    raise IOError('Cannot read Excel file from stdin')

  if compressed_filename(name):
    raise IOError('Cannot read compressed Excel file')

  book = xlrd.open_workbook(name)

  try:
    sheet = book.sheet_by_name(sheet)
  except xlrd.XLRDError:
    try:
      sheet = book.sheet_by_index(tryint1(sheet or 0))
    except (xlrd.XLRDError,TypeError,IndexError):
      raise ValueError('Cannot open Excel sheet %s:%s' % (namefile(name),sheet))

  def _table_reader_excel(book,sheet,strdata):
    if strdata:
      rowfunc = _xlate_xls_row_str
    else:
      rowfunc = _xlate_xls_row_object

    for i in xrange(sheet.nrows):
      values = sheet.row_values(i)
      types  = sheet.row_types(i)
      yield rowfunc(book,values,types)

  return _table_reader_excel(book,sheet,strdata)


class ExcelWriter(object):
  '''
  Write selected columns to a lower-level tabular data writer object

  Supports Excel versions: 2003, 2002, XP, 2000, 97, 95, 5.0, 4.0, 3.0, but not
  Excel 2007 XML (XLSX).
  '''
  def __init__(self, filename, extra_args=None, **kwargs):
    '''
    @param filename: file name or file object
    @type  filename: str or file object
    @param    sheet: sheet name
    @type     sheet: str
    '''
    try:
      import xlwt
    except ImportError:
      raise ValueError('Missing xlwt module to write Microsoft Excel file')

    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    name  = parse_augmented_filename(filename,args)

    sheet = get_arg(args, ['sheet'])
    hyout = get_arg(args, ['hyphen'])

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    self.filename = name if name!='-' or hyout is None else hyout
    self.book     = xlwt.Workbook()
    self.sheet    = self.book.add_sheet(sheet or 'Data')
    self.rownum   = 0

  def writerows(self, rows):
    '''
    Write a sequence of rows as rows in an Excel file

    @param rows: rows to write to Excel
    @type  rows: sequence of sequences of str
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    rownum = self.rownum
    write  = self.sheet.write

    for row in rows:
      for i,value in enumerate(row):
        write(rownum, i, value)
      rownum += 1

    self.row = rownum

  def writerow(self, row):
    '''
    Write a sequence of strings to a row in an Excel file

    @param row: row to write to Excel
    @type  row: sequences of str
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    rownum = self.rownum
    write  = self.sheet.write

    for i,value in enumerate(row):
      write(rownum, i, value)

    self.rownum += 1

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    self.book,book = None,self.book

    book.save(self.filename)

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()

  def __del__(self):
    '''
    Finish saving output when the writer is destroyed, if not already saved.
    '''
    if getattr(self,'book',None) is not None:
      self.close()


def cook_table(table, options):
  '''
  Create categorical variables, subset and sort table
  '''
  vars = ['categorical','columnexpr','includevar','excludevar','filterexpr','sort','uniq']

  if not any(getattr(options,var) for var in vars):
    return table

  table = iter(table)

  try:
    header = table.next()
  except StopIteration:
    return []

  if getattr(options,'categorical',None):
    header,table = create_categorical_variables(header,table,options.categorical)

  if getattr(options, 'columnexpr', None):
    header,table = column_exprs(header,table,options.filterexpr)

  if getattr(options,'includevar',None) or getattr(options,'excludevar',None):
    header,table = subset_variables(header,table,options.includevar,options.excludevar)

  if getattr(options, 'filterexpr', None):
    header,table = filter_expr(header,table,options.filterexpr)

  if getattr(options,'sort',None):
    header,table = sort_table(header,table,options.sort)

  if getattr(options,'uniq',None):
    header,table = uniq_table(header,table)

  return chain([header],table)


def _key_func(indices):
  # Itemgetter returns a single element for one key
  if len(indices) == 1:
    index = indices[0]
    def get_keys(item):
      return tryfloat(item[index])

  # Itemgetter returns a tuple for a compound key
  else:
    get_inds = itemgetter(*indices)
    def get_keys(item):
      return tuple(map(tryfloat, get_inds(item)))

  return get_keys


def sort_table(header, table, keys):
  '''
  Sort a table based on one or more keys

  >>> header =  ['a','b','c']
  >>> data   = [['3','M','-9' ],
  ...           ['2','F','9e9' ],
  ...           ['1','?','abc']]

  >>> h,d = sort_table(header,data,'a')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', '?', 'abc']
  ['2', 'F', '9e9']
  ['3', 'M', '-9']

  >>> h,d = sort_table(header,data,'*')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', '?', 'abc']
  ['2', 'F', '9e9']
  ['3', 'M', '-9']

  >>> h,d = sort_table(header,data,['b','a'])
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', '?', 'abc']
  ['2', 'F', '9e9']
  ['3', 'M', '-9']

  >>> h,d = sort_table(header,data,['c','*'])
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['1', '?', 'abc']

  >>> h,d = sort_table(header,data,'c,a')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['1', '?', 'abc']
  '''
  if is_str(keys):
    keys = [keys]

  indices  = []
  for key in keys:
    if key == '*':
      indices.extend( range(len(header)) )
    else:
      indices.extend(resolve_column_headers(header, key))

  if not indices:
    return header,table

  return header,sorted(table,key=_key_func(indices))


def uniq_table(header, table, keys=None):
  '''
  Generator to produce the unique first occurrence of each item in a table.
  Ordering is stable, since result elements will always appear in the order
  they first first appear in the input sequence.

  >>> header =  ['a','b','c']
  >>> data   = [['3','M','-9' ],
  ...           ['2','F','9e9' ],
  ...           ['3.0','M','-9' ],
  ...           ['3.0','F','-9' ],
  ...           ['2','F','9.0e9' ],
  ...           ['1','?','abc']]

  >>> h,d = uniq_table(header,data)
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['3.0', 'F', '-9']
  ['1', '?', 'abc']

  >>> h,d = uniq_table(header,data,'a')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['1', '?', 'abc']

  >>> h,d = uniq_table(header,data,'*')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['3.0', 'F', '-9']
  ['1', '?', 'abc']

  >>> h,d = uniq_table(header,data,['b','a'])
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['3.0', 'F', '-9']
  ['1', '?', 'abc']

  >>> h,d = uniq_table(header,data,'c,a')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', 'M', '-9']
  ['2', 'F', '9e9']
  ['1', '?', 'abc']
  '''
  if not keys:
    keys = range(len(header))
  elif is_str(keys):
    keys = [keys]

  indices  = []
  for key in keys:
    if key == '*':
      indices.extend( range(len(header)) )
    else:
      indices.extend(resolve_column_headers(header, key))

  if not indices:
    return header,table

  return header,unique(table,key=_key_func(indices))


def subset_variable(header,data,variable,include=None,exclude=None):
  '''
  Subset rows of a table based on inclusion and exclusion criteria for a single variable

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = subset_variable(header,data,'a',include='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = subset_variable(header,data,'b',exclude='?')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variable(header,data,'c',include='')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  '''
  index = resolve_column_header_atom(header,variable)

  if include is not None and exclude is not None:
    include = as_set(include)-as_set(exclude)
    exclude = None

  if include is not None:
    def _subset():
      for row in data:
        if row[index] in include:
          yield row

  elif exclude is not None:
    def _subset():
      for row in data:
        if row[index] not in exclude:
          yield row
  else:
    return header,data

  return header,_subset()


def subset_variables(header,data,include=None,exclude=None):
  '''
  Subset rows of a table based on inclusion and exclusion criteria for one or more variables

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = subset_variables(header,data,include='a=1,2')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variables(header,data,exclude='b=?')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variables(header,data,include='c')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']

  >>> h,d = subset_variables(header,data,include='c=')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variables(header,data,exclude='c')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variables(header,data,exclude='c=')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']
  '''
  if is_str(include):
    include = [include]
  else:
    include = include or []

  if is_str(exclude):
    exclude = [exclude]
  else:
    exclude = exclude or []

  for invar in include:
    if '=' not in invar:
      header,data = subset_variable(header,data,invar,exclude='')
    else:
      var,values = invar.split('=',1)
      values = set(v.strip() for v in values.split(','))
      header,data = subset_variable(header,data,var,include=values)

  for exvar in exclude:
    if '=' not in exvar:
      header,data = subset_variable(header,data,exvar,include='')
    else:
      var,values = exvar.split('=',1)
      values = set(v.strip() for v in values.split(','))
      header,data = subset_variable(header,data,var,exclude=values)

  return header,data


def column_exprs(header,data,exprs,use_globals=None):
  '''
  Create a new column based on a Python expression.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = column_exprs(header,data,"d=1 if fields[a] in ('1','2') else 0")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1']
  ['2', 'F', '', '1']
  ['3', '?', '1', '0']

  >>> h,d = column_exprs(header,data,["d=int(fields[a])**2","e=1 if fields[a] in ('1','2') else 0"])
  >>> h
  ['a', 'b', 'c', 'd', 'e']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '1']
  ['2', 'F', '', '4', '1']
  ['3', '?', '1', '9', '0']
  '''
  if is_str(exprs):
    exprs = [exprs]
  else:
    exprs = exprs or []

  for cexpr in exprs:
    var,expr = cexpr.split('=',1)
    header,data = column_expr(header,data,var,expr,use_globals)

  return header,data


def _make_expr_env(header,use_globals=None):
  use_globals = use_globals or {'__builtins__':__builtins__}
  indices     = use_globals['indices'] = {}
  counts      = tally(header)

  for i,h in enumerate(header):
    if counts[h] == 1:
      indices[h] = i

      spaces = ' ' in h
      digit  = h[0] in '0123456789'

      # Try to fix up h
      if spaces or digit:
        if spaces:
          h = h.replace(' ','_')
        if digit:
          h = '_' + h

        # If fixups make h ambiguous, do not set
        if h in counts:
          continue

      use_globals[h] = i

  return use_globals


def column_expr(header,data,column,expr,use_globals=None):
  '''
  Create a new column based on a Python expression.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = column_expr(header,data,'d',"1 if fields[a] in ('1','2') else 0")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1']
  ['2', 'F', '', '1']
  ['3', '?', '1', '0']

  >>> h,d = column_expr(header,data,'d',"(int(fields[a])+2)**2")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '9']
  ['2', 'F', '', '16']
  ['3', '?', '1', '25']
  '''
  code        = compile(expr,'<expression: %s>' % expr, 'eval')
  use_globals = _make_expr_env(header,use_globals)

  def _column_expr(data, code, use_globals):
    for row in data:
      use_globals['fields'] = row
      result = eval(code, use_globals)
      yield row + [str(result)]

  header = header + [column]
  return header,_column_expr(data, code, use_globals)


def filter_expr(header,data,expr,use_globals=None):
  '''
  Subset rows of a table based on a Python expression that evaluates to True
  or False.  Only rows that evaluate to True are retained.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = filter_expr(header,data,"fields[a] in ('1','2')")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"int(fields[a]) <= 2 and fields[b]=='M'")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = filter_expr(header,data,["int(fields[a]) <= 2","fields[b]=='M'"])
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = filter_expr(header,data,"fields[b]!='?'")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"fields[c]")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']

  >>> h,d = filter_expr(header,data,"not fields[c]")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"fields[2]==''")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"fields[c]")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']
  '''
  if not is_str(expr):
    expr      = ' and '.join('(%s)' % e for e in expr)

  code        = compile(expr,'<expression: %s>' % expr, 'eval')
  use_globals = _make_expr_env(header,use_globals)

  def _filter_expr(data, code, use_globals):
    for row in data:
      use_globals['fields'] = row
      if eval(code, use_globals):
        yield row

  return header,_filter_expr(data, code, use_globals)


def create_categorical_variable(header,data,variable,prefix=None,ref=None,include=None,exclude=None,
                              yes='1',no='0',missing=''):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = create_categorical_variable(header,data,'a')
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0']
  ['2', 'F', '', '0', '1', '0']
  ['3', '?', '1', '0', '0', '1']

  >>> h,d = create_categorical_variable(header,data,'a',ref=['1'])
  >>> h
  ['a', 'b', 'c', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '0', '0']
  ['2', 'F', '', '1', '0']
  ['3', '?', '1', '0', '1']

  >>> h,d = create_categorical_variable(header,data,'b',prefix='',exclude=['?'],yes='Y',no='N')
  >>> h
  ['a', 'b', 'c', 'F', 'M']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', 'N', 'Y']
  ['2', 'F', '', 'Y', 'N']
  ['3', '?', '1', '', '']

  >>> h,d = create_categorical_variable(header,data,'c',ref='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  ['3', '?', '1']

  >>> h,d = create_categorical_variable(header,data,'c',exclude='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  ['3', '?', '1']
  '''
  index = resolve_column_header_atom(header,variable)
  if include is not None and exclude is not None:
    include = as_set(include)-as_set(exclude)
    exclude = None

  if not isinstance(data,(tuple,list)):
    data = list(data)

  values = set(row[index] for row in data if len(row)>=index)
  values.discard('')

  if include is not None:
    values &= as_set(include)
  if exclude is not None:
    values -= as_set(exclude)

  if prefix is None:
    prefix = '%s_' % variable

  if ref is None:
    ref = set()
  else:
    ref = as_set(ref)

  values = sorted(values-ref)

  if not values:
    return header,data

  header = header + [ '%s%s' % (prefix,value) for value in sorted(values) ]
  values = dict( (v,i) for i,v in enumerate(values) )

  def _make_category():
    n = len(values)
    missingrow = [missing]*n
    for row in data:
      if index>=len(row):
        yield row+missingrow
        continue

      val = row[index]
      cats = [no]*n
      if val in values:
        cats[ values[val] ] = yes
      elif val not in ref:
        cats = missingrow

      yield row+cats

  return header,_make_category()


def create_categorical_variables(header,phenos,categorical):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = create_categorical_variables(header,data,['a'])
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0']
  ['2', 'F', '', '0', '1', '0']
  ['3', '?', '1', '0', '0', '1']

  >>> h,d = create_categorical_variables(header,data,['a','b:ref=M:prefix=:exclude=?','c:exclude=1:yes=Y:no=N'])
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3', 'F']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0', '0']
  ['2', 'F', '', '0', '1', '0', '1']
  ['3', '?', '1', '0', '0', '1', '']
  '''
  allowed_args = set(['prefix','ref','include','exclude','missing','yes','no'])

  for cat in categorical:
    opts = {}
    var  = parse_augmented_name(cat,opts)

    illegal = set(opts) - allowed_args
    if illegal:
      raise ValueError('Illegal argument(s) to categorical: %s' % ','.join(sorted(illegal)))

    for arg,val in opts.items():
      if arg in ('ref','include','exclude'):
        opts[arg] = set(v.strip() for v in val.split(','))

    header,phenos = create_categorical_variable(header,phenos,var,**opts)

  return header,phenos


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
