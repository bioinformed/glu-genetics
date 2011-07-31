# -*- coding: utf-8 -*-

__abstract__  = 'Functions to perform cross-tabulations of data'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   operator  import itemgetter
from   itertools import groupby, chain


def _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc):
  for row in data:
    rowkey = rowkeyfunc(row)
    colkey = colkeyfunc(row)
    value  = valuefunc(row)
    yield rowkey,colkey,value


def _infer_columns(data):
  '''
  Infer columns from a data stream by looking ahead in the stream to detect
  the first row key change.  Once found, the column keys and look-ahead
  buffer chained together with the remaining data generator are chained
  together and returned.

  This method of detecting columns can accommodate repeated column keys, but
  not missing column keys.  The first set of data records must have at least
  one instance of each column key.
  '''
  seen       = set()
  columns    = []
  look_ahead = []

  # Obtain the first row
  row = next(data,None)

  if row is None:
    return columns,data

  # Record the details
  look_ahead.append(row)
  rowkey = row[0]
  colkey = row[1]
  seen.add(colkey)
  columns.append(colkey)

  # Process the remaining rows for the first rowkey
  for row in data:
    look_ahead.append(row)

    if row[0]!=rowkey:
      break

    # Record first instance of each column key
    colkey = row[1]
    if colkey not in seen:
      seen.add(colkey)
      columns.append(colkey)

  # Return columns and entire data sequence
  return columns,chain(look_ahead,data)


def rowsby(data, columns=None, rowkeyfunc=None, colkeyfunc=None, valuefunc=None, aggregatefunc=None):
  '''
  Build genomatrix from genotriples when the necessary columns are given.
  Genotypes from the same row and column are grouped together.  Empty list
  will be used for missing genotypes.

  The columns that is passed in will be returned first followed by a
  generator of the resulting table.

  @param          data: genotriple stream
  @type           data: sequence
  @param       columns: sequence of column names
  @type        columns: sequence of strs
  @param    rowkeyfunc: function to identify the rowkey
  @type     rowkeyfunc: function
  @param    colkeyfunc: function to identify the colkey
  @type     colkeyfunc: function
  @param     valuefunc: function to identify the value
  @type      valuefunc: function
  @param aggregatefunc: function to aggregate columns
  @type  aggregatefunc: function
  @return             : genomatrix
  @rtype              : genomatrix generator

  >>> columns = ['l1','l2']
  >>> data = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
  ...         ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
  ...         ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
  >>> cols, rows = rowsby(data, columns, itemgetter(0), itemgetter(1), itemgetter(2))
  >>> print cols
  ['l1', 'l2']
  >>> for row in rows:
  ...   print row
  ('s1', [[('G', 'G')], [('A', 'A')]])
  ('s2', [[('G', 'T')], [('T', 'T')]])
  ('s3', [[('G', 'G')], [('A', 'A')]])
  >>> cols, rows = rowsby(data)
  >>> print cols
  ['l1', 'l2']
  >>> for row in rows:
  ...   print row
  ('s1', [[('G', 'G')], [('A', 'A')]])
  ('s2', [[('G', 'T')], [('T', 'T')]])
  ('s3', [[('G', 'G')], [('A', 'A')]])
  >>> data = [('l1','s1','AA'),('l1','s1',   0),('l1','s2','AB'),('l2','s1','AA'),
  ...         ('l2','s1','AA'),('l3','s1','BB'),('l3','s1','BB'),('l3','s1','AB')]
  >>> columns = ['s1','s2']
  >>> cols, rows = rowsby(data, columns, itemgetter(0), itemgetter(1), itemgetter(2))
  >>> print cols
  ['s1', 's2']
  >>> for row in rows:
  ...   print row
  ('l1', [['AA', 0], ['AB']])
  ('l2', [['AA', 'AA'], []])
  ('l3', [['BB', 'BB', 'AB'], []])
  '''
  rowkeyfunc = rowkeyfunc or itemgetter(0)
  colkeyfunc = colkeyfunc or itemgetter(1)
  valuefunc  = valuefunc  or itemgetter(2)

  data = _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc)

  if columns is None:
    columns,data = _infer_columns(data)

  # Build column key index
  columns = list(columns)
  colkeys = dict( (colkey,j) for j,colkey in enumerate(columns))

  if len(colkeys) != len(columns):
    raise ValueError('Column values must be unique')

  def _rowsby():
    # Build and yield result rows
    get0 = itemgetter(0)
    for rowkey,rowdata in groupby(data, get0):
      row = [ [] for i in xrange(len(colkeys)) ]

      for rowkey,colkey,value in rowdata:
        j = colkeys[colkey]
        row[j].append(value)

      if aggregatefunc:
        for colkey,j in colkeys.iteritems():
          row[j] = aggregatefunc(rowkey, colkey, row[j])

      yield rowkey,row

  return columns,_rowsby()


def xtab(data, rowkeyfunc=None, colkeyfunc=None, valuefunc=None,
               aggregatefunc=None, roworder=None, colorder=None):
  '''
  Build a table of rows by columns from a sequence that can be interpreted
  as a triple of (row key, column key, value).  Values with the same row and
  column keys are grouped together.  Cells with no values default to None.

  The list of colkeys and rowkeys will be returned first followed by a
  generator of the table rows.  Row and column keys appear in the optional
  roworder and colorder parameters.  All other keys will be in same order
  they are first observed in the data sequence.

  @param          data: genotriple stream
  @type           data: sequence
  @param    rowkeyfunc: function to identify the rowkey
  @type     rowkeyfunc: function
  @param    colkeyfunc: function to identify the colkey
  @type     colkeyfunc: function
  @param     valuefunc: function to identify the value
  @type      valuefunc: function
  @param aggregatefunc: function to aggregate columns
  @type  aggregatefunc: function
  @param      roworder: ordered row keys
  @type       roworder: sequence
  @param      colorder: ordered column keys
  @type       colorder: sequence
  @return             : genomatrix
  @rtype              : genomatrix generator

  >>> data = [('s1','l1', ('G', 'G')),('s1','l2', ('A', 'A')),
  ...         ('s2','l1', ('G', 'T')),('s2','l2', ('T', 'T')),
  ...         ('s3','l1', ('G', 'G')),('s3','l2', ('A', 'A'))]
  >>> cols, rows, vals = xtab(data, itemgetter(0), itemgetter(1), itemgetter(2))
  >>> print cols
  ['l1', 'l2']
  >>> print rows
  ['s1', 's2', 's3']
  >>> for val in vals:
  ...   print val
  [[('G', 'G')], [('A', 'A')]]
  [[('G', 'T')], [('T', 'T')]]
  [[('G', 'G')], [('A', 'A')]]
  >>> data = [('l1','s1','AA'),('l1','s1',   0),('l1','s2','AB'),('l2','s1','AA'),
  ...         ('l2','s1','AA'),('l3','s1','BB'),('l3','s1','BB'),('l3','s1','AB')]
  >>> cols, rows, vals = xtab(data, itemgetter(0), itemgetter(1), itemgetter(2))
  >>> print cols
  ['s1', 's2']
  >>> print rows
  ['l1', 'l2', 'l3']
  >>> for val in vals:
  ...   print val
  [[0, 'AA'], ['AB']]
  [['AA', 'AA'], None]
  [['AB', 'BB', 'BB'], None]
  '''
  get0 = itemgetter(0)
  get1 = itemgetter(1)
  get2 = itemgetter(2)

  rowkeyfunc = rowkeyfunc or get0
  colkeyfunc = colkeyfunc or get1
  valuefunc  = valuefunc  or get2

  rowkeys  = dict( (rowkey,i) for i,rowkey in enumerate(roworder or []) )
  colkeys  = dict( (colkey,i) for i,colkey in enumerate(colorder or []) )
  datalist = []

  # Pass 1: Build row, column, and data list
  for row in data:
    rowkey = rowkeyfunc(row)
    colkey = colkeyfunc(row)
    value  = valuefunc(row)
    i=rowkeys.setdefault(rowkey, len(rowkeys))
    j=colkeys.setdefault(colkey, len(colkeys))
    datalist.append( (i,j,value) )

  # Invert and sort the row and column keys
  rowkeys = sorted(rowkeys.iteritems(), key=get1)
  colkeys = sorted(colkeys.iteritems(), key=get1)

  datalist.sort()

  # Output column metadata
  columns = map(get0, colkeys)
  rows    = map(get0, rowkeys)

  def _xtab():
    # Pass 2: Build and yield result rows
    for i,rowdata in groupby(datalist, get0):
      row = [None]*len(colkeys)

      for j,vs in groupby(rowdata, get1):
        row[j] = map(get2, vs)

      if aggregatefunc:
        for colkey,j in colkeys:
          row[j] = aggregatefunc(rows[i], colkey, row[j])

      yield row

  return columns,rows,_xtab()


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
