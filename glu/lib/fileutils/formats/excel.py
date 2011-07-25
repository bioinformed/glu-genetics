# -*- coding: utf-8 -*-

__abstract__  = 'fileutils Excel file support'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   itertools                import izip

from   glu.lib.fileutils.auto   import compressed_filename, namefile
from   glu.lib.fileutils.parser import parse_augmented_filename, tryint1, get_arg


XLS_NULLDATE = (0,0,0)


def _xlate_xls_row_object(book,values,types):
  '''
  Translate a sequence of native Excel values and types into Python objects
  '''
  import xlrd
  from   datetime  import date,time,datetime

  row = []
  t = set(types)
  t -= set([xlrd.XL_CELL_NUMBER,xlrd.XL_CELL_DATE,xlrd.XL_CELL_BOOLEAN,xlrd.XL_CELL_ERROR])
  for value,typ in izip(values,types):
    if typ == xlrd.XL_CELL_NUMBER:
      ivalue = int(value)
      if ivalue == value:
        value = ivalue
    elif typ == xlrd.XL_CELL_DATE:
      if value[:3] == XLS_NULLDATE:
        value = time(*value[3:])
      elif value[3:] == XLS_NULLDATE:
        value = date(*value[:3])
      else:
        value = datetime(*value)
    elif typ == xlrd.XL_CELL_BOOLEAN:
      value = bool(value)
    elif typ == xlrd.XL_CELL_ERROR:
      value = xlrd.error_text_from_code[value]
    elif isinstance(value,unicode):
      value = value.encode('utf8')

    row.append(value)

  return row


def _xlate_xls_row_str(book,values,types):
  '''
  Translate a sequence of native Excel values and types into strings
  '''
  import xlrd

  row = []
  t = set(types)
  t -= set([xlrd.XL_CELL_NUMBER,xlrd.XL_CELL_DATE,xlrd.XL_CELL_BOOLEAN,xlrd.XL_CELL_ERROR])
  for value,typ in izip(values,types):
    if typ == xlrd.XL_CELL_NUMBER:
      ivalue = int(value)
      if ivalue == value:
        value = '%d' % value
      else:
        value = '%f' % value
    elif typ == xlrd.XL_CELL_DATE:
      value = xlrd.xldate_as_tuple(value,book.datemode)
      if value[:3] == XLS_NULLDATE:
        value = '%02d:%02d:%02d' % value[3:]
      elif value[3:] == XLS_NULLDATE:
        value = '%04d/%02d/%02d' % value[:3]
      else:
        value = '%04d/%02d/%02d %02d:%02d:%02d' % value
    elif typ == xlrd.XL_CELL_BOOLEAN:
      value = ['0','1'][value]
    elif typ == xlrd.XL_CELL_ERROR:
      value = xlrd.error_text_from_code[value]
    elif isinstance(value,unicode):
      value = value.encode('utf8')

    row.append(value)

  return row


def table_reader_excel(filename,strdata=True,extra_args=None,**kwargs):
  '''
  Load rows from a Microsoft Excel (XLS) file using the xlrd module

  Supports Excel versions: 2003, 2002, XP, 2000, 97, 95, 5.0, 4.0, 3.0, but not
  Excel 2007 XML (XLSX).
  '''
  try:
    import xlrd
  except ImportError:
    raise ValueError('Missing xlrd module to read Microsoft Excel file')

  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name  = parse_augmented_filename(filename,args)
  sheet = get_arg(args, ['sheet'],0)
  hyin  = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if hyin is not None and name=='-':
    raise IOError('Cannot read Excel file from stdin')

  if compressed_filename(name):
    raise IOError('Cannot read compressed Excel file')

  book = xlrd.open_workbook(name,on_demand=True)

  try:
    sheet = book.sheet_by_name(sheet)
  except xlrd.XLRDError:
    try:
      sheet = book.sheet_by_index(tryint1(sheet or 0))
    except (xlrd.XLRDError,TypeError,IndexError):
      raise ValueError('Cannot open Excel sheet %s:%s' % (namefile(name),sheet))

  def _table_reader_excel(book,sheet,strdata):
    if strdata:
      rowfunc = _xlate_xls_row_str
    else:
      rowfunc = _xlate_xls_row_object

    for i in xrange(sheet.nrows):
      values = sheet.row_values(i)
      types  = sheet.row_types(i)
      yield rowfunc(book,values,types)

  return _table_reader_excel(book,sheet,strdata)


class ExcelWriter(object):
  '''
  Write selected columns to a lower-level tabular data writer object

  Supports Excel versions: 2003, 2002, XP, 2000, 97, 95, 5.0, 4.0, 3.0, but not
  Excel 2007 XML (XLSX).
  '''
  def __init__(self, filename, extra_args=None, **kwargs):
    '''
    @param filename: file name or file object
    @type  filename: str or file object
    @param    sheet: sheet name
    @type     sheet: str
    '''
    try:
      import xlwt
    except ImportError:
      raise ValueError('Missing xlwt module to write Microsoft Excel file')

    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    name  = parse_augmented_filename(filename,args)

    sheet = get_arg(args, ['sheet'])
    hyout = get_arg(args, ['hyphen'])

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    self.filename = name if name!='-' or hyout is None else hyout
    self.book     = xlwt.Workbook()
    self.sheet    = self.book.add_sheet(sheet or 'Data')
    self.rownum   = 0

  def writerows(self, rows):
    '''
    Write a sequence of rows as rows in an Excel file

    @param rows: rows to write to Excel
    @type  rows: sequence of sequences of str
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    rownum = self.rownum
    write  = self.sheet.write

    for row in rows:
      for i,value in enumerate(row):
        write(rownum, i, value)
      rownum += 1

    self.row = rownum

  def writerow(self, row):
    '''
    Write a sequence of strings to a row in an Excel file

    @param row: row to write to Excel
    @type  row: sequences of str
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    rownum = self.rownum
    write  = self.sheet.write

    for i,value in enumerate(row):
      write(rownum, i, value)

    self.rownum += 1

  def close(self):
    '''
    Close the writer

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.book is None:
      raise IOError('Writer object already closed')

    self.book,book = None,self.book

    book.save(self.filename)

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
    Finish saving output when the writer is destroyed, if not already saved.
    '''
    if getattr(self,'book',None) is not None:
      self.close()


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
