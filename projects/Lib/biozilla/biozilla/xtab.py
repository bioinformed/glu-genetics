import sys
from   operator  import itemgetter
from   itertools import groupby,chain,izip


def _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc):
  for row in data:
    rowkey = rowkeyfunc(row)
    colkey = colkeyfunc(row)
    value  = valuefunc(row)
    yield rowkey,colkey,value


def rowsby(data, columns, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=None):
  get0 = itemgetter(0)
  get1 = itemgetter(1)
  get2 = itemgetter(2)

  columns = list(columns)

  data = _rowstream(data, rowkeyfunc, colkeyfunc, valuefunc)

  # Build column key index
  colkeys = dict( (colkey,j) for j,colkey in enumerate(columns))

  if len(colkeys) != len(columns):
    raise ValueError, 'Column values must be unique'

  # Output column metadata
  yield columns

  # Build and yield result rows
  for rowkey,rowdata in groupby(data, get0):
    row = [ [] for i in range(len(colkeys)) ]

    for rowkey,colkey,value in rowdata:
      j = colkeys[colkey]
      row[j].append(value)

    if aggregatefunc:
      for colkey,j in colkeys.iteritems():
        row[j] = aggregatefunc(rowkey, colkey, row[j] or None)

    yield rowkey,row


def xtab(data, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=None):
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


def xtab_list(data, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=None):
  columns,rows,results = xtab(data, rowkeyfunc, colkeyfunc, valuefunc, aggregatefunc=aggregatefunc)
  return chain([columns],izip(rows,results))
