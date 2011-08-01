# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Utility to merge rows with matching values from two files'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   collections         import defaultdict

from   glu.lib.utils       import pick,unique
from   glu.lib.fileutils   import table_reader,tryint1,table_writer,resolve_column_headers, \
                                  cook_table,table_options


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('table1', help='First table to join')
  parser.add_argument('table2', help='Second table to join')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')

  parser.add_argument('-1', '--key1',
                    help='Table 1 key column numbers or names.  Comma separated, with '
                         'defaulting to common headers')
  parser.add_argument('-2', '--key2',
                    help='Table 2 key column numbers or names')
  parser.add_argument('--prefix1', metavar='P1',
                    help='prefix to prepend to each non-key header name in table1')
  parser.add_argument('--prefix2', metavar='P2',
                    help='prefix to prepend to each non-key header name in table2')
  parser.add_argument('-U', '--uniquejoin', action='store_true',
                    help='Require that rows with valid keys from table2 are unique')
  parser.add_argument('-j', '--join', default='left', choices=['left','inner'],
                    help="Join type: 'left' (default) or 'inner'")

  table_options(parser)

  return parser


def _parse_header(table,key):
  table = iter(table)

  try:
    header = table.next()
  except StopIteration:
    return [],[],None

  if isinstance(key, str):
    key = [ tryint1(k) for k in key.split(',') ]
  elif isinstance(key, int):
    key = [ key ]

  if key is not None:
    key = tuple(unique(resolve_column_headers(header,key)))

  return header,table,key


def left_join(table1,table2,key1=None,key2=None,unique=False,inner=False,
              null='',prefix1=None,prefix2=None):
  '''
  A generator that performs a relational join on two tables streams based on
  key equality (see http://en.wikipedia.org/wiki/Join_(SQL) ).  Both input
  tables are required to be an iterable sequence of a header and subsequent
  data rows, all of equal length.

  If no keys are specified, an equi-join is assumed on all matching headers.
  Otherwise, keys may be specified as:

  1) A sequence of header names, integer column indices (zero-based), string
     column indices (zero-based)

  or

  2) a single comma delimited string of header names or column indices
     (zero-based)

  If inner is False, then a relational LEFT OUTER JOIN is performed such
  that all rows in table1 will appear at least once in the results.  table1
  rows without any matching table2 rows will be emitted with all null values
  for table2 fields.

  If inner is True, then a relational INNER JOIN is performed such that only
  rows in table1 with one or more matching rows in table2 will be emitted.

  Row keys containing null values, as determined by the null parameter, will
  not match any other key.  All keys must be hashable since the underlying
  implementation is a hash join (see
  http://en.wikipedia.org/wiki/Hash_join).

  The results contain all fields from table 1 and all non-key fields from
  table2 and begins with a header and each subsequent result data row.

  @param table1: Left input table
  @type  table1: Sequence of header and data rows, all of equal row length
  @param table2: Right input table
  @type  table2: Sequence of header and data rows, all of equal row length
  @param   key1: Key selector for table1, if None an equi-join is performed
  @type    key1: Sequence of mixed int or str, or str, or None
  @param   key2: Key selector for table1, if None an equi-join is performed
  @type    key2: Sequence of mixed int or str, or str, or None
  @param unique: If True, raise an exception if table2 keys are not unique
  @type  unique: bool
  @param  inner: If True perform a INNER JOIN, otherwise a LEFT OUTER JOIN
  @type   inner: bool
  @param   null: Null data value
  @type    null: object
  @return      : join results with header and data rows
  @rtype       : generator

  >>> table1 = [[ 'A', 'B'],
  ...           [   1,   2],
  ...           [   1,   3],
  ...           [   2,   1],
  ...           [None,None]]
  >>> table2 = [[ 'B', 'C'],
  ...           [   3,   1],
  ...           [None,   1],
  ...           [   2,   0]]

  Left outer equi-join:

  >>> for row in left_join(table1,table2,null=None):
  ...   print row
  ['A', 'B', 'C']
  [1, 2, 0]
  [1, 3, 1]
  [2, 1, None]
  [None, None, None]

  Left outer join with specified keys:

  >>> for row in left_join(table1,table2,key1=1,key2=0,null=None):
  ...   print row
  ['A', 'B', 'C']
  [1, 2, 0]
  [1, 3, 1]
  [2, 1, None]
  [None, None, None]

  Inner equi-join:

  >>> for row in left_join(table1,table2,null=None,inner=True):
  ...   print row
  ['A', 'B', 'C']
  [1, 2, 0]
  [1, 3, 1]
  '''
  header1,table1,index1 = _parse_header(table1, key1)
  header2,table2,index2 = _parse_header(table2, key2)

  if index1 is None:
    index1 = index2
  if index2 is None:
    index2 = index1

  if index1 is None and index2 is None:
    equikeys = set(header1) & set(header2)
    index1 = tuple(i for i,h in enumerate(header1) if h in equikeys)
    index2 = tuple(i for i,h in enumerate(header2) if h in equikeys)

  if not index1 or not index2:
    raise ValueError('No compatible keys found')

  if len(index1) != len(index2):
    raise ValueError('Keys must be of same length')

  index2r   = [ i for i in range(len(header2)) if i not in index2 ]
  header2   = pick(header2,index2r)
  len1,len2 = len(header1),len(header2)
  max1,max2 = max(index1),max(index2)

  if prefix1:
    header1 = [ prefix1+h if i not in index1 else h for i,h in enumerate(header1) ]

  if prefix2:
    header2 = [ prefix2+h for h in header2 ]

  # Compute the left hash-join

  # Step 1: Extract keys and hash all rows in table2
  rowmap = defaultdict(list)
  for row in table2:
    # Skip the row if it is too short to have all keys
    if len(row) <= max2:
      continue

    key = tuple(pick(row,index2))

    # Keys are null values treated as NULL
    if any(k==null for k in key):
      continue

    if unique and row in rowmap[key]:
      raise ValueError('Non-unique key %s' % (key,))

    # Extract non-key fields and normalize length
    row  = pick(row,index2r)
    row += [null]*(len2-len(row))
    row  = row[:len2]

    rowmap[key].append(row)

  # Step 2: Perform the join over all rows in table 1

  # Yield the new header
  yield header1 + header2

  # This default value for non-key matches is what enforces the LEFT OUTER
  # or INNER JOIN semantics.
  blank2 = [[null]*len2] if not inner else []

  for row in table1:
    key  = tuple(pick(row,index1)) if len(row) > max1 else None

    # Normalize length
    row += [null]*(len1-len(row))
    row  = row[:len1]

    # Combine with
    for row2 in rowmap.get(key,blank2):
      yield row+row2


def main():
  parser  = option_parser()
  options = parser.parse_args()

  table1  = table_reader(options.table1,hyphen=sys.stdin,want_header=True)
  table2  = table_reader(options.table2,hyphen=sys.stdin,want_header=True)

  table   = left_join(table1,table2,key1=options.key1,key2=options.key2,
                      unique=options.uniquejoin,inner=(join_type=='inner'),
                      prefix1=options.prefix1,prefix2=options.prefix2)

  out     = table_writer(options.output,hyphen=sys.stdout)

  table   = cook_table(table,options)

  out.writerows(table)


if __name__=='__main__':
  import doctest
  doctest.testmod()
  main()
