from db_row import *

try:
  Nothing
except NameError:
  Nothing = object()


def check_expected_rows(rowcount, expectedRows, minRows, maxRows):
  if rowcount < 0:
    return

  if expectedRows is not None and rowcount != expectedRows:
    if rowcount > expectedRows:
      raise dbexceptions.TooManyRowsError
    else:
      raise dbexceptions.TooFewRowsError

  if minRows is not None and rowcount < minRows:
    raise dbexceptions.TooFewRowsError
  if maxRows is not None and rowcount > maxRows:
    raise dbexceptions.TooManyRowsError


def query(cursor, sql, hashby=None, orderby=None, project=None, single=False,
                       default=Nothing, expectedRows=None, minRows=None, maxRows=None):

  if single and not hashby:
    expectedRows = 1

  cursor.execute(sql)
  rows = cursor.fetchall()
  rowcount = len(rows)

  if rowcount == 0 and default is not Nothing:
    return default
  else:
    check_expected_rows(rowcount, expectedRows, minRows, maxRows)
  rowclass = db_row.IMetaRow(cursor.description)
  rows = db_row.RowList([ rowclass(row) for row in rows ], row_class=rowclass)

  return uberize(rows, hashby, orderby, project, single)
