# -*- coding: utf-8 -*-

__abstract__  = 'file related utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   operator                 import itemgetter
from   itertools                import chain

from   glu.lib.utils            import is_str,as_set,unique,Counter
from   glu.lib.fileutils.parser import parse_augmented_name, tryfloat
from   glu.lib.fileutils.table  import resolve_column_headers, resolve_column_header_atom


__all__ = ['cook_table', 'sort_table', 'uniq_table', 'subset_variables',
           'create_categorical_variables', 'column_exprs', 'filter_expr',
           'query_table','table_options']


def table_options(parser):
  parser.add_argument('-c', '--categorical', metavar='VAR', action='append',
                      help='Create indicator variables based on values of VAR')
  parser.add_argument('--columnexpr', metavar='VAR=EXPR', action='append',
                      help='Add a new column VAR with the value determined by expression EXPR')
  parser.add_argument('--includevar', metavar='VAR=VAL', action='append',
                      help='Include only records with variable VAR equal to VAL')
  parser.add_argument('--excludevar', metavar='VAR=VAL', action='append',
                      help='Exclude all records with variable VAR equal to VAL')
  parser.add_argument('--filterexpr', metavar='EXPR', action='append',
                      help='Filter all records where EXPR is not true')
  parser.add_argument('--query', metavar='SQL', action='append',
                      help='Execute in-memory SQL query on all records')
  parser.add_argument('-s', '--sort', metavar='VAR', action='append',
                      help='Sort rows based on values in column VAR')
  parser.add_argument('-u', '--uniq', action='store_true',
                      help='Produce only unique rows by collapsing consecutive duplicate rows')
  return parser


def cook_table(table, options):
  '''
  Create categorical variables, subset and sort table
  '''
  vars = ['categorical','columnexpr','includevar','excludevar','filterexpr','sort','uniq','query']

  if not any(getattr(options,var,None) for var in vars):
    return table

  table = iter(table)

  try:
    header = table.next()
  except StopIteration:
    return []

  if getattr(options,'categorical',None):
    header,table = create_categorical_variables(header,table,options.categorical)

  if getattr(options, 'columnexpr', None):
    header,table = column_exprs(header,table,options.columnexpr)

  if getattr(options, 'query'):
    header,table = query_table_list(header,table,options.query)

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


def column_exprs(header,data,exprs,env=None):
  '''
  Create a new column based on a Python expression.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = column_exprs(header,data,"d=1 if a in ('1','2') else 0")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1']
  ['2', 'F', '', '1']
  ['3', '?', '1', '0']

  >>> h,d = column_exprs(header,data,["d=int(a)**2","e=1 if a in ('1','2') else 0"])
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
    header,data = column_expr(header,data,var,expr,env)

  return header,data


def _make_expr_env(code,header,env=None):
  env     = env or {'__builtins__':__builtins__}
  indices = env['indices'] = {}
  counts  = Counter(header)
  headers = {}

  for i,h in enumerate(header):
    if not h:
      continue

    if counts[h] == 1:
      indices[h] = i

      illegal = ' ' in h or '-' in h
      digit   = h[0] in '0123456789'

      # Try to fix up h
      if illegal or digit:
        if illegal:
          h = h.replace(' ','_').replace('-','_')
        if digit:
          h = '_' + h

        # If fixups make h ambiguous, do not set
        if h in counts:
          continue

      headers[h] = i

  fields = [ (name,headers[name]) for name in code.co_names if name in headers ]

  return env,fields


def column_expr(header,data,column,expr,env=None):
  '''
  Create a new column based on a Python expression.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = column_expr(header,data,'d',"1 if a in ('1','2') else 0")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1']
  ['2', 'F', '', '1']
  ['3', '?', '1', '0']

  >>> h,d = column_expr(header,data,'d',"(int(a)+2)**2")
  >>> h
  ['a', 'b', 'c', 'd']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '9']
  ['2', 'F', '', '16']
  ['3', '?', '1', '25']
  '''
  code       = compile(expr,'<expression: %s>' % expr, 'eval')
  env,fields = _make_expr_env(code,header,env)

  def _column_expr(data, code, fields, env):
    for row in data:
      # Create environment bindings
      env['fields'] = row
      for name,index in fields:
        env[name] = row[index]

      result = eval(code, env)

      yield row + [str(result)]

  header = header + [column]
  return header,_column_expr(data, code, fields, env)


def filter_expr(header,data,expr,env=None):
  '''
  Subset rows of a table based on a Python expression that evaluates to True
  or False.  Only rows that evaluate to True are retained.

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = filter_expr(header,data,"a in ('1','2')")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"int(a) <= 2 and b=='M'")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = filter_expr(header,data,["int(a) <= 2","b=='M'"])
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = filter_expr(header,data,"b!='?'")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = filter_expr(header,data,"c")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']

  >>> h,d = filter_expr(header,data,"not c")
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

  >>> h,d = filter_expr(header,data,"c")
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']
  '''
  if not is_str(expr):
    expr     = ' and '.join('(%s)' % e for e in expr)

  code       = compile(expr,'<expression: %s>' % expr, 'eval')
  env,fields = _make_expr_env(code,header,env)

  def _filter_expr(data, code, fields, env):
    for row in data:
      # Create environment bindings
      env['fields'] = row
      for name,index in fields:
        env[name] = row[index]

      if eval(code, env):
        yield row

  return header,_filter_expr(data, code, fields, env)


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


def create_categorical_variables(header,data,categorical):
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

    header,data = create_categorical_variable(header,data,var,**opts)

  return header,data


def query_table(header,data,query):
  '''
  Execute an SQL query on a table using an in-memory temporary SQLite table

  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = query_table(header,data,'SELECT * FROM data')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  ['3', '?', '1']

  >>> h,d = query_table(header,data,'SELECT * FROM data WHERE a="3"')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']

  >>> h,d = query_table(header,data,'SELECT b,c FROM data WHERE a="3"')
  >>> h
  ['b', 'c']
  >>> for row in d:
  ...   print row
  ['?', '1']
  '''
  import sqlite3
  from   glu.lib.fileutils.formats.sqlite import SQLiteWriter, table_reader_sqlite

  con = sqlite3.connect(':memory:')

  writer = SQLiteWriter('',con=con,table='data')
  writer.writerow(header)
  writer.writerows(data)

  reader = iter(table_reader_sqlite('',con=con,query=query))
  return reader.next(),reader


def query_table_list(header,data,queries):
  '''
  Execute multiple sequential SQL table queries.

  The current implementation isn't ideal, since in-memory tables are
  re-created for each query.
  '''
  for query in queries:
    header,data = query_table(header,data,query)

  return header,data


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
