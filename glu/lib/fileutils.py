# -*- coding: utf-8 -*-
'''
File:          fileutils.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-02-21

Abstract:      file related utility functions

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import os
import csv

from   itertools     import islice,chain
from   operator      import itemgetter

from   glu.lib.utils import peekfirst


__all__ = ['autofile','namefile','hyphen',
           'guess_format','related_file','guess_related_file',
           'load_list','load_map','load_table','table_writer']


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
  if not isstr(filename):
    return filename

  if args is None:
    args = {}

  if parse_augmented_filename(filename,args) == '-':
    return defaultfile

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

  if not isstr(filename):
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


_unescape_literal_chars = [('\\r','\r'),('\\t','\t'),('\\n','\n'),('\\\\','\\')]

def unescape_literal(s):
  '''
  Evaluate C/Python character literals \r, \t, \n, and \\ into the
  respective single-character ASCII representation
  '''
  for e,c in _unescape_literal_chars:
    s = s.replace(e,c)
  return s


def tryint(s):
  '''
  Try to coerce an arbitrary object to an integer, otherwise return the
  original value

  @param s: arbitrary item
  @type  s: object
  @return : integer coerced value or the orginal value
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
  @return : integer coerced value or the orginal value
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
  ValueError: Invalid zero string index used where 1-based index is required
  '''
  try:
    ss = int(s)
  except (ValueError,TypeError):
    return s

  if isinstance(s, basestring):
    if ss==0:
      raise ValueError('Invalid zero string index used where 1-based index is required')
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


def load_list(filename,extra_args=None,**kwargs):
  '''
  Load list of values from a file or literal list.

  By default, a file is interpreted as a tab-delimited file with only the
  first field considered.  Blank lines are skipped, whitespace is stripped
  from the beginning and end of each field selected.  These settings can be
  modified via function arguments or by appending additional options to the
  filename.  File loading is delegated to the load_table function, so refer
  to load_table for details on supported arguments.

  The column number or name may be specified as per the load_table function.
  However, for backward compatibility the 'index' parameter is also
  supported as a synonym to 'columns'.

  If the filename is a string that begins with ':', it is interpreted as a
  literal list of comma separated values, ignoring all other parameters
  (including skip!).  No escaping or quoting is allowed.

  The following parameters and aliases are accepted as part of the augmented
  filename (all but the first from load_table):

  [by load_list]
       index, i: column index or header name

  [by load_table]
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
  >>> load_list(f.name + ':skip=1:index=1')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=\\t:index=H1')
  ['v11', 'v21']
  >>> load_list(f.name + ':dialect=excel-tab:skip=1:index=2')
  ['v12', 'v32']
  >>> load_list(f.name + ':index=H2')
  ['v12', 'v32']
  >>> load_list(f.name + ':c=H2')
  ['v12', 'v32']
  >>> load_list(f.name + ':c=2:skip=1')
  ['v12', 'v32']
  >>> load_list(f.name + ':i=2:c=1:skip=1')
  Traceback (most recent call last):
     ...
  ValueError: Invalid specification of both index and columns
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write('H1,H2\\nv11,v12\\nv21\\n,v32')
  >>> f.flush()
  >>> load_list(f.name + ':delimiter=,:skip=1')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=,:skip=1:index=1')
  ['v11', 'v21']
  >>> load_list(f.name + ':delimiter=,:index=H1')
  ['v11', 'v21']
  >>> load_list(f.name + ':dialect=excel:skip=1:index=2')
  ['v12', 'v32']
  >>> load_list(f.name + ':dialect=excel:index=H2')
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

  items = load_table(name, columns=[index], extra_args=args)

  # Needed to trigger arg processing since generators do not start until the
  # first item is requested
  items = chain([items.next()],items)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  return [ intern(i[0]) for i in items if i and i[0] ]


# Anonymous object
_nothing = object()


def load_map(filename,unique=True,extra_args=None,**kwargs):
  '''
  Creates a dictionary representing a list or mapping from a text file.

  By default, valid files should be composed of standard ASCII lines of text
  with tab delimited fields.  Only the first and second fields of each line
  are considered. Whitespace is stripped from the beginning and end of every
  field considered.  If the skip parameter is used to ignore a certain
  number of lines (e.g., headers) at the beginning of the file.  A default
  parameter may be specified to assign values to keys with empty or
  non-existant value fields.  Otherwise, the value will be set equal to the
  key.  File loading is delegated to the load_table function, so refer
  to load_table for details on supported arguments.

  The key and value column number or name may be specified as per the
  load_table function using the 'columns' argument.  The first and second
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

  [by load_map]
      key_index, i: index of field to select or name of header for keys
    value_index, v: index of field to select or name of header for values
     default,def,d: default value for keys with no or empty value

  [by load_table]
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
    items = (kv.split('=') for kv in filename[1:].split(',') if kv)
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

    items = load_table(name,columns=columns,extra_args=args)

    # Needed to trigger arg processing since generators do not start until the
    # first item is requested
    items = chain([items.next()],items)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  # Parse the data file
  def _load_map_generator():
    for row in items:
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

  if isstr(column):
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


def resolve_column_headers(header,columns):
  '''
  @param       header: header line of the input file
  @type        header: sequence of strs
  @param      columns: indices, names, or ranges of columns to select, comma
                       delimited
  @type       columns: list of strings, integers, or 2-tuples for ranges
  @return            : resolved header line
  @rtype             : sequence of strs
  '''
  if isinstance(columns,int):
    columns = [columns]
  elif isinstance(columns,str):
    columns = columns.split(',')

  indices = []
  for column in columns:
    col = resolve_column_header(header,column)
    if isinstance(col,tuple):
      indices.extend( xrange(col[0],col[1]+1) )
    else:
      indices.append(col)
  return indices


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
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[(0,'c2')])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=[('1','c2')])
  [['loc1', ''], ['loc2', ''], ['loc1', 'loc1'], ['loc2', '']]
  >>> test(StringIO("c1\\tc2\\nloc1\\nloc2\\nloc1\\tloc1\\nloc2"),columns=['c2','c1'])
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
  dialect = get_csv_dialect(args,'tsv')
  skip    = int(get_arg(args, ['skip','s'], 0))
  columns = get_arg(args, ['c','cols','columns'])
  hyin    = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name) if name!='-' else hyin
  lines = csv.reader(lfile,**dialect)

  if skip:
    lines = islice(lines,skip,None)

  # All columns are to be returned
  if not columns:
    try:
      row = lines.next()
    except StopIteration:
      return

    # Process row 1 (may or may not be a header)
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
  try:
    indices = resolve_column_headers(None,columns)
  except ValueError:
    # Need headers
    try:
      header = map(str.strip,lines.next())
    except StopIteration:
      return

    indices = resolve_column_headers(header,columns)

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


class TableWriter(object):
  '''
  Write selected columns to a lower-level tabular data writer object
  '''
  def __init__(self, writer, columns):
    self.writer  = writer
    self.columns = columns

    # Try to resolve column headings without a header row
    try:
      self.indices = resolve_column_headers(None,columns)
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

      self.indices = indices = resolve_column_headers(row,self.columns)

    for row in rows:
      m = len(row)
      row = [ (row[j] if j<m else '') for j in indices ]
      self.writer.writerow(row)

  def writerow(self, row):
    indices = self.indices

    # First row and indices require header information
    if indices is None:
      self.indices = indices = resolve_column_headers(row,self.columns)

    m = len(row)
    row = [ (row[j] if j<m else '') for j in indices ]
    self.writer.writerow(row)


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
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  columns = get_arg(args, ['c','cols','columns'])
  dialect = get_csv_dialect(args,'tsv')
  hyout   = get_arg(args,['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if name == '-':
    outfile = hyout
  else:
    outfile = autofile(name,'w')

  writer = csv.writer(outfile, **dialect)

  if columns is not None:
    writer = TableWriter(writer, columns)

  return writer


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
