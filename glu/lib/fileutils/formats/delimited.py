# -*- coding: utf-8 -*-

__abstract__  = 'fileutils delmited text file support'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import csv

from   glu.lib.fileutils.auto   import guess_format, autofile
from   glu.lib.fileutils.parser import parse_augmented_filename, trybool, get_arg


__all__ = ['get_csv_dialect','table_reader_delimited','delimited_table_writer']


# Create more standard aliases for Python CSV module dialects.  'excel' is
# now available as 'csv' and 'excel-tab' is now available as 'tsv'.
csv.register_dialect('csv', csv.get_dialect('excel'))
csv.register_dialect('tsv', csv.get_dialect('excel-tab'))

DIALECT_KWARGS = ['dialect','delimiter','doublequote','escapechar','lineterminator',
                  'quotechar','quoting','skipinitialspace','strict']

_unescape_literal_chars = [('\\r','\r'),('\\t','\t'),('\\n','\n'),('\\\\','\\')]


def unescape_literal(s):
  '''
  Evaluate C/Python character literals \r, \t, \n, and \\ into the
  respective single-character ASCII representation
  '''
  for e,c in _unescape_literal_chars:
    s = s.replace(e,c)
  return s


def get_csv_dialect(args,default_dialect='tsv'):
  '''
  Retrieve the standardized csv argument list

  @param             args: csv arguments
  @type              args: dict
  @param  default_dialect: set the default value for dialect to 'tsv'
  @type   default_dialect: str
  @return                : list of tuples of (option,value)
  @rtype                 : list of tuples

  >>> args = dict(dialect='csv',quoting='none',skipinitialspace='y',doublequote='false',unrelated='baz',lineterminator='\\n',
  ...             delimiter='\\t')
  >>> sorted(get_csv_dialect(args).iteritems())
  [('delimiter', '\\t'), ('dialect', 'csv'), ('doublequote', False), ('lineterminator', '\\n'), ('quoting', 3), ('skipinitialspace', True)]
  >>> sorted(args.iteritems())
  [('unrelated', 'baz')]
  '''
  dargs = {'dialect':default_dialect}
  for darg in DIALECT_KWARGS:
    if darg in args:
      dargs[darg] = unescape_literal(args.pop(darg))

  if 'quoting' in dargs:
    q = 'QUOTE_' + dargs['quoting'].upper()
    dargs['quoting'] = getattr(csv, q, q)

  # Parse truth values
  for tf in ('skipinitialspace','doublequote'):
    if tf in dargs:
      dargs[tf] = trybool(dargs[tf])

  return dargs


def table_reader_delimited(filename, extra_args=None, **kwargs):
  '''
  Return a configured delimited table reader
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  dialect = get_csv_dialect(args, guess_format(name, ['csv']) or 'tsv')
  hyin    = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name) if name!='-' or hyin is None else hyin

  return csv.reader(lfile,**dialect)


def delimited_table_writer(filename, extra_args=None, **kwargs):
  '''
  Return a configured delimited table writer
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  name    = parse_augmented_filename(filename,args)
  dialect = get_csv_dialect(args, guess_format(name, ['csv']) or 'tsv')
  hyout   = get_arg(args, ['hyphen'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  lfile = autofile(name,'wb') if name!='-' or hyout is None else hyout

  return csv.writer(lfile,**dialect)


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
