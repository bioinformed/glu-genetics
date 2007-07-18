# -*- coding: utf-8 -*-
'''
File:          rowops.py

Author:        Kevin Jacobs (jacobs@theopalgroup.com)

Created:       October 17, 2002

Purpose:       Operators over database rows

Compatibility: Python 2.2 and above

Requires:

Revision:      $Id: rowops.py 202 2006-06-13 18:17:19Z jacobske $

Copyright (c) 2001,2002,2003 The OPAL Group.  All rights reserved.
'''

import sys

if __name__ == '__main__':
  sys.path.append('../../..')


__all__ = ['unique', 'extract_values', 'project', 'projection', 'extend',
           'hash_unique', 'hash_list', 'order_by', 'hash_aggregate', 'hash_sum',
           'aggregate_by', 'sum_by', 'uberize', 'extract_names']


import operator
import db_row
import dbexceptions


try:
  Nothing
except NameError:
  Nothing = object()


def iskey(v):
  return isinstance(v, (str,int)) or callable(v)


def extract_names(fields, descr=None):
  used = {}
  names = []
  i = 1
  for field in fields:
    name = field
    if isinstance(field, int) and descr:
      name = descr.__fieldnames__[field]
    elif callable(field) and getattr(field,'__doc__',None):
      name = field.__doc__

    while name in used:
      name = '_field%d' % i
      i += 1

    used[name] = True
    names.append(name)

  return names


class FieldNamer(object):
  __slots__ = ('i','names')
  def __init__(self):
    self.i = 1
    self.names = {}

  def name_field(self, descr, field):
    if isinstance(field, str):
      name = descr[field].name
    elif isinstance(field, int):
      name = descr[descr.__fieldnames__[field]].name
    elif callable(field) and getattr(field,'__doc__',None):
      name = field.__doc__
    else:
      name = self.get_anon()

    self.names[name] = True
    return name

  def name_extra(self, descr, field):
    if isinstance(field, str):
      name = field
    elif callable(field) and getattr(field,'__doc__',None):
      name = field.__doc__
    else:
      name = self.get_anon()

    self.names[name] = True
    return name

  def get_anon(self):
    while True:
      name = '_field%d' % self.i
      self.i += 1
      if name not in self.names:
        return name


def extract_descrs(descr, fields, extra_fields=None):
  if fields is None:
    fields = project(descr, 'name')

  namer = FieldNamer()
  descrs = [ namer.name_field(descr, field) for field in fields ]

  if extra_fields is not None:
    if iskey(extra_fields):
      extra_fields = [extra_fields]
    descrs += [ namer.name_extra(descr, field) for field in extra_fields ]

  return descrs


def get_row_class(rows):
  row_class = None
  if isinstance(rows, db_row.RowList):
    row_class = rows.row_class
  elif rows and isinstance(rows[0], db_row.Row):
    row_class = type(rows[0])
  return row_class


def project_row_class(rows, fields, extra_fields=None):

  row_class = get_row_class(rows)

  if not row_class or iskey(fields):
    return tuple

  descrs = extract_descrs(row_class.field_descriptors, fields, extra_fields)

  if issubclass(row_class, db_row.IRow):
    row_class = db_row.IMetaRow(descrs,driver=row_class.driver)
  else:
    row_class = db_row.MetaRow(descrs,driver=row_classs.driver)

  return row_class


def promote_rows(original, derived, fields):
  row_class = project_row_class(original, fields)
  if not row_class:
    return derived
  return db_row.RowList([ row_class(row) for row in derived ], row_class)


def unique(rows, reorder=True):
  '''unique(rows) -> [sorted distinct/unique elements of rows]'''

  if not rows:
    return rows

  if reorder:
    rows.sort()
  out = rows[0:1]
  for i in range(1,len(rows)):
    if rows[i] != out[-1]:
      out.append(rows[i])

  return out


def extract_values(row, fields, row_class=tuple):
  '''extract_values(row, fields)

    where row is a sequence or mapping object,
      and fields is an integer, a string, a callable object that takes
          a row as a parameter, or a tuple consisting of a mix of tuples,
          integers and strings.

    returns a values or row (or a value computed from the row), possibly
    nested within sets of tuples.

    This function performs a recursive extraction of fields from row
    based on a template given by the fields argument.

    e.g.: extract_values(['a','b','c'], 0)         -> 'a'
          extract_values(['a','b','c'], (1,0))     -> ('b','a')
          extract_values(['a','b','c'], (2,(1,0))) -> ('c',('b','a'))
          extract_values({1:'a',2:'b'},  1)        -> 'a'
          extract_values({1:'a',2:'b'}, (2,1))     -> ('b','a')
          extract_values({1:'a',2:'b'}, (1,(1,1))) -> ('a',('a','a'))
          extract_values([1,2,3],    ((0,max,min)) -> (1,3,1)
  '''

  # Handle non-recursive case
  if isinstance(fields, (str,int)):
    return row[fields]
  elif callable(fields):
    return fields(row)
  else:
    return row_class([ extract_values(row, field) for field in fields ])


def projection(rows, fields):
  '''project(rows, fields)

     where rows are a list of sequences or mappings,
           fields is a (possibly recursive) sequence of keys, indices that
           index the elements of rows, or callable objects (see
           extract_values).

     This function returns a new row list by computing the specified new
     fields from each row.
  '''
  row_class = project_row_class(rows, fields)
  return [ extract_values(row, fields, row_class) for row in rows ]


project = projection


def extend(rows, *args, **kwargs):
  '''extend(rows, field1=value1, field2=value2, ...)

     where rows are a list of sequences or mappings,
           fields is a (possibly recursive) sequence of keys, indices that
           index the elements of rows, or callable objects (see
           extract_values).

     This function returns a new row list built by extending each row by
     adding the fields specified.
  '''
  names  = []
  fields = []

  for field in args:
    if callable(field) and field.__doc__:
      names.append(field.__doc__)
    else:
      names.append(field)
    fields.append(field)

  names  += kwargs.keys()
  fields += kwargs.values()

  row_class = project_row_class(rows, None, names)
  return [ row_class(row+extract_values(row, fields)) for row in rows ]

#######################################################################################

def hash_unique(rows, key_fields, value_fields=None):
  '''hash_unique(rows, key_fields, [value_fields])

     where rows are a list of sequences or mappings,
           key_fields is a (possibly recursive) sequence of keys, indices that
           index the elements of rows, or callable objects (see
           extract_values) (see extract_values), and (optionally)
           value_fields to extract from each row (again, see extract_values).

     This function builds an arbitrary dimensional hash table based on the
     values extracted from key_fields, mapping to the original rows or
     values extracted from them (if value_fields are specified).  The last
     level of the hash is required to be unique or else a ValueError
     exception will be thrown.
  '''

  if iskey(key_fields):
    key_fields = [key_fields]

  if value_fields is not None:
    row_class = project_row_class(rows, value_fields)

  root = {}
  inner_fields = key_fields[:-1]
  for row in rows:
    results = root
    for key_field in inner_fields:
      key = extract_values(row, key_field)
      results = results.setdefault(key, {})
    key = extract_values(row, key_fields[-1])
    if key in results:
      raise ValueError, 'Non-unique key (%s) extracted from row %s' % \
                        (extract_values(row, key_fields), row)
    if value_fields is not None:
      value = extract_values(row, value_fields, row_class)
    else:
      value = row
    results[key] = value

  return root


def hash_list(rows, key_fields, value_fields=None):
  '''hash_list(rows, key_fields, [value_fields])

     where rows are a list of sequences or mappings, key_fields is a
           (possibly recursive) sequence of keys, indices that index the
           elements of rows, or callable objects (see extract_values), and
           (optionally) value_fields to extract from each row (again see
           extract_values).

     This function builds an arbitrary dimensional hash table based on the
     values extracted from key_fields, mapping to a list of the original
     rows or values extracted from them (if value_fields are specified).
     The last level of the hash table is not required to be unique and a
     list of rows of values extracted from rows will always be stored.
  '''

  if iskey(key_fields):
    key_fields = [key_fields]

  if value_fields is not None:
    row_class = project_row_class(rows, value_fields)

  root = {}
  inner_fields = key_fields[:-1]
  for row in rows:
    results = root
    for key_field in inner_fields:
      key = extract_values(row, key_field)
      results = results.setdefault(key, {})
    key = extract_values(row, key_fields[-1])

    if value_fields is not None:
      value = extract_values(row, value_fields, row_class)
    else:
      value = row
    results.setdefault(key,[]).append(value)

  return root

#######################################################################################

def order_by(rows, keys):
  '''order_by(rows, key_fields)

     where rows are a list of sequences or mappings, keys is a (possibly
           recursive) sequence of keys, indices that index the elements of
           rows, or callable objects (see extract_values).

     This function sorts rows by the fields extracted by the keys specified.
     An old technique called a Schwartzian Transform is used to make this
     function very fast.

     NOTE 1: This function DOES result in a stable ordering.
     NOTE 2: This function currently supports only ascending sort orders.
     NOTE 3: Extracting computed keys by passing in callable keys is
             efficient, in the sense that each computed key is calculated
             exactly once.  Some implementations would naively recompute keys
             for each pairwise comparison during the sort, making
             the use of computed keys prohibitively slow.
  '''

  rows = [ (extract_values(row,keys),i,row) for i,row in enumerate(rows) ]
  rows.sort()
  return [ x[-1] for x in rows ]

#######################################################################################

def reduce_hash_list(function, hash):
  '''reduce_hash_list(function, hash)

     where function is a binary aggregating function,
           and hash is an arbitrary dimension list hash.

     This function performs data reduction on the lists at the lowest level
     of the hash (see hash_list and builtin reduce for more information).

     NOTE: This function IS destructive and will overwrite the supplied data.
  '''
  stack = [hash]
  while stack:
    top = stack.pop()
    for key,value in top.iteritems():
      if isinstance(value, dict):
        stack.append(value)
      else:
        top[key] = reduce(function, value)
  return hash


def hash_aggregate(function, rows, keys, value_fields=None):
  '''hash_aggregate(function, rows, keys, [value_fields])

     where function is a binary aggregating function,
           rows are a list of sequences or mappings,
           key_fields is a (possibly recursive) sequence of keys or
           indices that index the elements of rows (see extract_values),
           and (optionally) value_fields to extract from each row (see
           extract_values).

     This function builds an arbitrary dimensional hash table based on the
     values extracted from key_fields, mapping to a value that is determined by the
     reduction of all elements with the selected keys from the original
     rows or values extracted from them (if value_fields are specified).

     e.g.: a = [('a',1), ('b',1), ('a',1), ('c',1), ('a',1), ('b',1), ('a',1)]
           hash_aggregate(operator.add, a, 0, 1) -> {'a': 4, 'c': 1, 'b': 2}

           b = [('a','a','a',1),('b','b','b',1),('a','a','a',1),('c','c','c',1),
                ('a','a','a',1),('b','b','b',1),('a','a','a',1)]
           hash_aggregate(operator.add, b, (0,1,2), -1)
           -> {'a': {'a': {'a': 4}},
               'c': {'c': {'c': 1}},
               'b': {'b': {'b': 2}} }
  '''
  hash = hash_list(rows, keys, value_fields)
  return reduce_hash_list(function, hash)


def hash_sum(rows, keys, value_fields=None):
  '''hash_sum(rows, keys, [value_fields])

     where rows are a list of sequences or mappings,
           key_fields is a (possibly recursive) sequence of keys or
           indices that index the elements of rows (see extract_values),
           and (optionally) value_fields to extract from each row (see
           extract_values).

     This function builds an arbitrary dimensional hash table based on the
     values extracted from key_fields, mapping to the sum of all elements
     with the selected keys from the original rows or values extracted from
     them (if value_fields are specified).

     NOTE: hash_sum(x,y) === hash_aggregate(operator.add, x, y)

     e.g.: a = [('a',1), ('b',1), ('a',1), ('c',1), ('a',1), ('b',1), ('a',1)]
           hash_sum(a, 0, 1) -> {'a': 4, 'c': 1, 'b': 2}

           b = [('a','a','a',1),('b','b','b',1),('a','a','a',1),('c','c','c',1),
                ('a','a','a',1),('b','b','b',1),('a','a','a',1)]
           hash_sum(b, (0,1,2), -1)
           -> {'a': {'a': {'a': 4}},
               'c': {'c': {'c': 1}},
               'b': {'b': {'b': 2}} }
  '''
  return hash_aggregate(operator.add, rows, keys, value_fields)

#######################################################################################

def _vector_op(functions):
  def _vector_op(lvalues,rvalues):
    return [ f(l,r) for f,l,r in zip(functions,lvalues,rvalues) ]
  return _vector_op


def aggregate_by(function, rows, keys, value_fields, name='aggregate'):
  '''aggregate_by(function, rows, keys, value_fields, [name])

     where function is a binary aggregating function,
           rows are a list of sequences or mappings,
           key_fields is a (possibly recursive) sequence of keys or
           indices that index the elements of rows (see extract_values),
           and a value_field to extract from each row (see extract_values).

     This function builds a new row list based on the values extracted from
     key_fields and the the reduction of all corresponding value fields
     extracted from the original rows.

     e.g.: a = [('a',1), ('b',1), ('a',1), ('c',1), ('a',1), ('b',1), ('a',1)]
           hash_aggregate(operator.add, a, 0, 1) -> [('a', 4), ('b', 2), ('c', 1)]

           b = [('a','a','a',1),('b','b','b',1),('a','a','a',1),('c','c','c',1),
                ('a','a','a',1),('b','b','b',1),('a','a','a',1)]
           hash_aggregate(operator.add, b, (0,1,2), -1)
           -> [('a', 'a', 'a', 4), ('c', 'c', 'c', 1), ('b', 'b', 'b', 2)]
  '''
  if iskey(keys):
    keys = [keys]

  if iskey(value_fields):
    multivalued  = False
    names = extract_names([value_fields])
  else:
    multivalued = True
    names = extract_names(value_fields)

  row_class = project_row_class(rows, keys, value_fields)
  hash = hash_aggregate(function, rows, [keys], value_fields)

  if multivalued:
    return [ row_class( group+tuple(value) ) for group,value in hash.iteritems() ]
  else:
    return [ row_class( group+(value,) ) for group,value in hash.iteritems() ]


def sum_by(rows, keys, value_fields, name='sum'):
  '''sum_by(rows, keys, value_fields, [name])

     where function is a binary aggregating function,
           rows are a list of sequences or mappings,
           key_fields is a (possibly recursive) sequence of keys or
           indices that index the elements of rows (see extract_values),
           and a value_field to extract from each row (see extract_values).

     This function builds a new row list based on the values extracted from
     key_fields and the the reduction of all corresponding value fields
     extracted from the original rows.

     e.g.: a = [('a',1), ('b',1), ('a',1), ('c',1), ('a',1), ('b',1), ('a',1)]
           hash_aggregate(operator.add, a, 0, 1) -> [('a', 4), ('b', 2), ('c', 1)]

           b = [('a','a','a',1),('b','b','b',1),('a','a','a',1),('c','c','c',1),
                ('a','a','a',1),('b','b','b',1),('a','a','a',1)]
           hash_aggregate(operator.add, b, (0,1,2), -1)
           -> [('a', 'a', 'a', 4), ('c', 'c', 'c', 1), ('b', 'b', 'b', 2)]
  '''
  if iskey(value_fields):
    function = operator.add
  else:
    function = _vector_op([operator.add]*len(value_fields))

  return aggregate_by(function, rows, keys, value_fields, name)


#######################################################################################

def build_merge_tuple(left_rows, right_rows, on_fields):
  left_row_class  = get_row_class(left_rows)
  right_row_class = get_row_class(right_rows)

  if left_row_class and right_row_class:

    identity_fields = [ f1 for f1,f2 in on_fields if isinstance(f1,str) and f1==f2 ]

    left_descrs     = left_row_class.field_descriptors
    left_driver     = left_row_class.driver

    right_descrs    = right_row_class.field_descriptors
    right_driver    = right_row_class.driver

    driver = (left_driver == right_driver and left_driver) or None

    left_descrs     = [ db_row.FieldDescriptor(d) for d in left_descrs  ]
    right_descrs    = [ db_row.FieldDescriptor(d) for d in right_descrs
                          if d['name'] not in identity_fields ]

    left_project    = project(left_descrs,  'name')
    right_project   = project(right_descrs, 'name')

    for f in left_descrs:
      if f['name'] in right_project:
        f['name'] = 'left_%s' % f['name']

    for f in right_descrs:
      if f['name'] in left_project:
        f['name'] = 'right_%s' % f['name']

    descrs = left_descrs + right_descrs

    if issubclass(left_row_class, db_row.IRow) or issubclass(right_row_class, db_row.IRow):
      row_class = db_row.IMetaRow(descrs,driver=driver)
    else:
      row_class = db_row.MetaRow(descrs,driver=driver)

    empty = (None,) * len(right_descrs)
    def merge_func(l,r, right_project=right_project, row_class=row_class, empty=empty):
      if r is not None:
        return row_class(l+extract_values(r, right_project))
      else:
        return row_class(l+empty)

  elif isinstance(left_rows[0], dict) and isinstance(right_rows[0], dict):
    def merge_func(l,r):
      m=dict(l)
      if r:
        m.update(r)
      return m
  else:
    empty = (None,) * len(right_rows[0])
    merge_func = lambda l,r,empty=empty: l + (r or empty)

  return merge_func


def build_join_hashes(left_rows, right_rows, on_fields):
  if on_fields:
    left_fields,right_fields = zip(*on_fields)
    if len(left_fields) == len(right_fields) == 1:
      left_fields  = left_fields[0]
      right_fields = right_fields[0]
    left  = hash_list(left_rows,  left_fields)
    right = hash_list(right_rows, right_fields)
  else:
    left  = { None : left_rows  }
    right = { None : right_rows }

  return left,right


def expand_join_condition(on_fields):
  if on_fields is None:
    return []
  elif iskey(on_fields):
    return [(on_fields,on_fields)]
  else:
    onf = []
    for f in on_fields:
      if iskey(f):
        onf.append( (f,f) )
      else:
        onf.append(f)
    return onf


def equi_join(left_rows, right_rows, on_fields, left_join=False):

  if not left_rows or not right_rows:
    return []

  on_fields  = expand_join_condition(on_fields)
  left,right = build_join_hashes(left_rows, right_rows, on_fields)
  merge_func = build_merge_tuple(left_rows, right_rows, on_fields)

  rows = []
  for l_key,l_rows in left.iteritems():
    if l_key is None or (isinstance(l_key, (list,tuple)) and None in l_key):
       continue
    elif l_key in right:
      r_rows = right[l_key]
    elif left_join:
      r_rows = [None]
    else:
      continue
    for l_row in l_rows:
      for r_row in r_rows:
        rows.append( merge_func(l_row,r_row))

  return rows

#######################################################################################

def uberize(rows, hashby=None, orderby=None, project=None, single=False, default=Nothing):
  # Handle default
  if not rows and default is not Nothing:
    return default

  # Order rows
  if orderby is not None:
    rows = order_by(rows, orderby)

  # Add grouping or projections
  if hashby is not None:
    # Enforce single-flag semantics
    if not single:
      rows = hash_list(rows, hashby, project)
    else:
      rows = hash_unique(rows, hashby, project)

  # otherwise, just apply projections
  elif project is not None:
    rows = projection(rows, project)

  # Enforce single-flag semantics if not grouping
  if single and not hashby:
    if not rows:
      raise dbexceptions.TooFewRowsError
    elif len(rows) > 1:
      raise dbexceptions.TooManyRowsError
    rows = rows[0]

  return rows



def main():
  a = [{'a':1, 'b':2},{'a':3,'b':4}]
  b = [{'a':1, 'c':5},{'a':1,'c':6},{'a':2,'c':7}]

  print "equi_join(%s, %s, 'a')\n  = %s" % (str(a).replace(' ',''),
                                            str(b).replace(' ',''),
                                            str(equi_join(a,b, 'a')).replace(' ',''))
  print

  print "equi_join(%s, %s, [('b','a'])\n  = %s" % (str(a).replace(' ',''),
                                            str(b).replace(' ',''),
                                            str(equi_join(a,b, [('b','a')])).replace(' ',''))
  print

  a = [(1,2),(3,4)]
  b = [(1,5),(1,6),(2,7)]

  print "equi_join(%s, %s, 0)\n  = %s" % (str(a).replace(' ',''),
                                          str(b).replace(' ',''),
                                          str(equi_join(a,b, 0)).replace(' ',''))
  print

  print "equi_join(%s, %s, [(1,0)])\n  = %s" % (str(a).replace(' ',''),
                                                str(b).replace(' ',''),
                                                str(equi_join(a,b, [(1,0)])).replace(' ',''))
  print

  Row1 = db_row.IMetaRow( ('a','b') )
  Row2 = db_row.IMetaRow( ('a','c') )
  a = map(Row1, a)
  b = map(Row2, b)

  print "equi_join(%s, %s, 'a')\n  = %s" % (str(a).replace(' ',''),
                                            str(b).replace(' ',''),
                                            str(equi_join(a,b, 'a')).replace(' ',''))
  print

  print "equi_join(%s, %s, [('b','a'])\n  = %s" % (str(a).replace(' ',''),
                                            str(b).replace(' ',''),
                                            str(equi_join(a,b, [('b','a')])).replace(' ',''))
  print [ d.dict() for d in equi_join(a,b, [('b','a')]) ]
  print

  print "left_join(%s, %s, [('b','a'])\n  = %s" % (str(a).replace(' ',''),
                                            str(b).replace(' ',''),
                                            str(equi_join(a,b, [('b','a')],left_join=1)).replace(' ',''))
  print [ d.dict() for d in equi_join(a,b, [('b','a')],left_join=1) ]
  print


def main2():
  import db_row

  r = [(1,2,3)]*5

  print extend(r)
  print extend(r, 0)
  print extend(r, new_a=0)

  R = db_row.IMetaRow( ['a','b','c'] )
  r = [R([1,2,3])]*5

  print extend(r)
  print extend(r, 0)
  print extend(r, new_a=0)

def main3():
  import db_row

  R = db_row.IMetaRow( ['a','b','c'] )
  r = [R([1,2,3])]*3

  print r
  print sum_by(r, 'a', 'b')
  print sum_by(r, ['a','b'], 'c')
  print sum_by(r, 'a', ['b','c'])

if __name__ == '__main__':
  main()
  main2()
  main3()
