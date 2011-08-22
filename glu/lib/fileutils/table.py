# -*- coding: utf-8 -*-

__abstract__  = 'file related utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os

from   itertools                import islice

from   glu.lib.utils            import is_str,deprecated_by
from   glu.lib.fileutils.auto   import guess_format, namefile
from   glu.lib.fileutils.parser import parse_augmented_filename, tryint1, get_arg


__all__ = ['list_reader', 'map_reader', 'table_reader', 'table_writer',
           'resolve_column_headers', 'resolve_column_header', 'table_columns']


TABLE_FORMATS = set(['xls','xlsx','csv','dta','stata','db','sqlite'])


def _literal_list(filename):
  return is_str(filename) and not os.path.exists(filename) and filename.startswith(':')


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


def resolve_column_headers(header,include=None,exclude=None):
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

  >>> resolve_column_headers(['a','b','c'],'*')
  [0, 1, 2]
  >>> resolve_column_headers(['a','b','c'],['c','*'],'b')
  [2, 0]
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
      # Wildcards add all un-used indices at that position once all headers are resolved
      if column=='*' and '*' not in header:
        indices.append('*')
        continue

      col = resolve_column_header(header,column)
      if isinstance(col,tuple):
        indices.extend( xrange(col[0],col[1]+1) )
      else:
        indices.append(col)

  # Add all un-used indices at the position of the first wildcard
  if '*' in indices:
    iset    = set(indices)
    iset.discard('*')
    star    = [ i for i in range(len(header)) if i not in iset ]
    index   = indices.index('*')
    indices = indices[:index]+star+[ i for i in indices[index+1:] if i != '*' ]

  # Remove indices for any desired exclusions
  if exclude:
    exclude_indices = set(resolve_column_headers(header,exclude))
    indices = [ i for i in indices if i not in exclude_indices ]

  return indices


def parse_rename(rename):
  if not rename:
    return None

  if isinstance(rename,str):
    rename = [ r.split('=') for r in rename.split(',') ]

  rdict = {}
  if isinstance(rename, (list,tuple)):
    for r in rename:
      if len(r) != 2:
        raise ValueError('Invalid column rename: %s' % str(r))
      rdict[r[0]] = r[1]
    rename = rdict

  return rename


def table_columns(rows, columns, drop=None, rename=None, header=None, want_header=False):
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
  rows   = iter(rows)
  rename = parse_rename(rename)

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

        if rename:
          row = [ rename.get(h,h) for h in row ]

        yield row

      else:
        header = map(str.strip,header)
        n      = len(header)
        if want_header:
          if rename:
            header = [ rename.get(h,h) for h in header ]
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

  def _table_reader_columns(header):
    if header is not None:
      header = [ header[j] for j in indices ]
      if rename:
        header = [ rename.get(h,h) for h in header ]
      yield header

    # Build result rows
    for row in rows:
      m = len(row)
      result  = [ (row[j].strip() if j<m else '') for j in indices ]
      yield result

  return _table_reader_columns(header)


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
  rename  = get_arg(args, ['r','rename'])
  header  = get_arg(args, ['h','header'])

  if header is not None and isinstance(header,str):
    header = map(str.strip,header.split(','))

  format = format.lower()
  if format in ('xls','excel'):
    from glu.lib.fileutils.formats.excel import table_reader_excel
    rows = table_reader_excel(name, extra_args=args)
  elif format in ('xlsx','excel2007','excel2010'):
    from glu.lib.fileutils.formats.xlsx import table_reader_xlsx
    rows = table_reader_xlsx(name, extra_args=args)
  elif format in ('delimited','tsv','csv'):
    from glu.lib.fileutils.formats.delimited import table_reader_delimited
    rows = table_reader_delimited(name, extra_args=args)
  elif format in ('db','sqlite'):
    from glu.lib.fileutils.formats.sqlite import table_reader_sqlite
    rows = table_reader_sqlite(name, extra_args=args)
  elif format in ('dta','stata'):
    from glu.lib.fileutils.formats.stata import table_reader_stata
    rows = table_reader_stata(name, extra_args=args)
  else:
    raise NotImplementedError("Reading file format '%s' is not supported" % format)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if skip:
    rows = islice(rows,skip,None)

  return table_columns(rows,columns,drop=drop,rename=rename,header=header,want_header=want_header)


class TableWriter(object):
  '''
  Write selected columns to a lower-level tabular data writer object
  '''
  def __init__(self, writer, columns, drop, rename):
    self.writer  = writer
    self.columns = columns
    self.drop    = drop
    self.header  = None
    self.rename  = parse_rename(rename)
    self.indices = None

  def _write_header(self, header):
    self.indices = indices = resolve_column_headers(header,self.columns,self.drop)

    m      = len(header)
    header = [ (header[j] if j<m else '') for j in indices ]

    if self.rename:
      rename = self.rename
      header = [ rename.get(h,h) for h in header ]

    self.writer.writerow(header)

  def writerows(self, rows):
    indices = self.indices

    if indices is None:
      rows = iter(rows)
      try:
        header = rows.next()
      except IndexError:
        return

      self._write_header(header)
      indices = self.indices

    for row in rows:
      m = len(row)
      row = [ (row[j] if j<m else '') for j in indices ]
      self.writer.writerow(row)

  def writerow(self, row):
    indices = self.indices

    # First row and indices require header information
    if indices is None:
      self._write_header(row)
    else:
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
  rename  = get_arg(args, ['r','rename'])

  format  = format.lower()
  if format in ('xls','excel'):
    from glu.lib.fileutils.formats.excel import ExcelWriter
    writer = ExcelWriter(name, extra_args=args)
  elif format in ('xlsx','excel2007','excel2010'):
    from glu.lib.fileutils.formats.xlsx import XLSXWriter
    writer = XLSXWriter(name, extra_args=args)
  elif format in ('sqlite','db'):
    from glu.lib.fileutils.formats.sqlite import SQLiteWriter
    writer = SQLiteWriter(name, extra_args=args)
  elif format in ('delimited','tsv','csv'):
    from glu.lib.fileutils.formats.delimited import delimited_table_writer
    writer = delimited_table_writer(name, extra_args=args)
  else:
    raise NotImplementedError("Writing to file format '%s' is not supported" % format)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  # Use a wrapper writer object to handle column selectors and header renaming
  # FIXME: Header renaming without column selectors could be much more
  #        efficient as a special case
  if columns is not None or drop or rename:
    writer = TableWriter(writer, columns, drop, rename)

  return writer


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
