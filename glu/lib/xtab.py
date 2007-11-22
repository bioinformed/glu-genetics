# -*- coding: utf-8 -*-
'''
File:          xtab.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
from   operator  import itemgetter
from   itertools import groupby,chain


def _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc):
  for row in data:
    rowkey = rowkeyfunc(row)
    colkey = colkeyfunc(row)
    value  = valuefunc(row)
    yield rowkey,colkey,value


def rowsby(data, columns, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=None):
  '''
  Build genomatrix from genotriples when the necessary columns are given. Genotypes from the same row and column are grouped   together. Empty list will be used for missing genotypes.

  The columns that is passed in will be returned first followed by A generator of the genomatrix which is built

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
  get0 = itemgetter(0)
  get1 = itemgetter(1)
  get2 = itemgetter(2)

  columns = list(columns)

  data = _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc)

  # Build column key index
  colkeys = dict( (colkey,j) for j,colkey in enumerate(columns))

  if len(colkeys) != len(columns):
    raise ValueError, 'Column values must be unique'

  def _rowsby():
    # Build and yield result rows
    for rowkey,rowdata in groupby(data, get0):
      row = [ [] for i in range(len(colkeys)) ]

      for rowkey,colkey,value in rowdata:
        j = colkeys[colkey]
        row[j].append(value)

      if aggregatefunc:
        for colkey,j in colkeys.iteritems():
          row[j] = aggregatefunc(rowkey, colkey, row[j])

      yield rowkey,row

  return columns,_rowsby()


def xtab(data, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=None):
  '''
  Build genomatrix from genotriples. Genotypes from the same row and column are grouped together.
  Missing genotypes default to None.

  The list of colkeys and rowkeys will be returned first followed by A generator of the genomatrix which is built

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

  rowkeys  = {}
  colkeys  = {}
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

