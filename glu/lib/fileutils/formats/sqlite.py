# -*- coding: utf-8 -*-

__abstract__  = 'fileutils sqlite format support'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sqlite3

from   itertools                import chain

from   glu.lib.fileutils.auto   import compressed_filename
from   glu.lib.fileutils.parser import parse_augmented_filename, get_arg


__all__ = ['table_reader_sqlite','SQLiteWriter']


def table_reader_sqlite(filename, con=None, extra_args=None, **kwargs):
  '''
  Return a configured delimited table reader
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name  = parse_augmented_filename(filename,args)
  table = get_arg(args, ['table','t'])
  query = get_arg(args, ['query','q'])
  hyin  = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if hyin is not None and name=='-':
    raise IOError('Cannot read Excel file from stdin')

  if compressed_filename(name):
    raise IOError('Cannot read compressed Excel file')

  if not (bool(table) ^ bool(query)):
    raise IOError('Either a table name or query must be specified')

  if table:
    query = 'SELECT * from %s' % table

  if not con:
    con   = sqlite3.connect(name)

  cur     = con.execute(query)
  header  = [ d[0] for d in cur.description ]
  data    = ( map(str,row) for row in cur )

  return chain([header],data)


class SQLiteWriter(object):
  def __init__(self, filename, con=None, extra_args=None, **kwargs):
    '''
    Return a configured SQLite table writer
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    name  = parse_augmented_filename(filename,args)
    table = get_arg(args, ['table','t'])
    hyout = get_arg(args, ['hyphen'])

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if hyout is not None and name=='-':
      raise IOError('Cannot write sqlite file to stdout')

    if compressed_filename(name):
      raise IOError('Cannot read compressed sqlite file')

    if not table:
      raise ValueError('Invalid sqlite table name specified')

    self.con    = con or sqlite3.connect(name)
    self.owncon = con is None
    self.table  = table
    self.header = None

  def create_table(self, header):
    if self.con is None:
      raise IOError('Writer object already closed')

    # FIXME: Escape SQL identifiers
    self.header = ["'%s'" % h for h in header]
    sql = '''CREATE TABLE IF NOT EXISTS '%s'(%s)''' % (self.table, ','.join(self.header))
    self.con.execute(sql)

  def writerows(self, rows):
    '''
    Write a sequence of rows as a table in a SQLite database

    @param rows: rows to write
    @type  rows: sequence of sequences of str
    '''
    if self.con is None:
      raise IOError('Writer object already closed')

    rows = iter(rows)
    if self.header is None:
      self.create_table(rows.next())

    header = self.header
    n      = len(header)

    nones  = [None]*n
    def _extend(row):
      m = len(row)
      if n == m:
        return row

      if m < n:
        row  = list(row)
        row += none[n-m]
      else:
        row  = row[:n]

      return row

    sql = '''INSERT INTO '%s'(%s) VALUES (%s)''' % (self.table, ','.join(header), ','.join('?'*n))
    self.con.executemany(sql, (_extend(row) for row in rows))

  def writerow(self, row):
    '''
    Write a row to a table in a SQLite database

    @param row: row to write
    @type  row: sequences of str
    '''
    if self.con is None:
      raise IOError('Writer object already closed')

    if self.header is None:
      self.create_table(row)
      return

    header = self.header
    n      = len(header)
    m      = len(row)

    if m < n:
      row  = list(row)
      row += [None]*(n-m)
    elif m > n:
      row  = row[:n]

    sql = '''INSERT INTO '%s'(%s) VALUES (%s)''' % (self.table, ','.join(header), ','.join('?'*n))
    self.con.execute(sql, row)

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.con is None:
      raise IOError('Writer object already closed')

    self.con,con = None,self.con

    con.commit()
    if self.owncon:
      con.close()

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()

  def __del__(self):
    '''
    Commit results when the writer is destroyed, if not already closed.
    '''
    if getattr(self,'con',None) is not None:
      self.close()


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
