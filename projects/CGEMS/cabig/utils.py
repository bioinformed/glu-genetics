import csv
from   itertools      import islice

from   biozilla.utils import autofile


#FIXME: Consider to enhance the similar function in biozilla.genodata
def load_map(file_or_name,keycol,valcol,skip=0,dialect='excel-tab'):
  '''
  @param keycol: the column index for the key and it is 0 based
  @type  keycol: int
  @param valcol: the column index for the value and it is 0 based
  @type  valcol: int
  '''
  rows = csv.reader(autofile(file_or_name), dialect=dialect)

  if skip:
    rows = islice(rows,skip,None)

  map = {}
  for row in rows:
    if row and len(row) > max(keycol,valcol):
      map[row[keycol]] = row[valcol]

  return map


#FIXME: Consider to enhance the similar function in biozilla.genodata
def load_list(file_or_name,col,skip=0,dialect='excel-tab'):
  '''
  @param col: the column index for the list entries and it is 0 based
  @type  col: int
  '''
  rows = csv.reader(autofile(file_or_name), dialect=dialect)

  if skip:
    rows = islice(rows,skip,None)

  return set(row[col] for row in rows if row and len(row)>col)


def load_rows(filename,skip=1,delimit='excel'):
  rows = csv.reader(autofile(filename),dialect=delimit)
  return islice(rows,skip,None)

