# -*- coding: utf-8 -*-
'''
File:          sections.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-07-13

Abstract:      Set of tools and parsers to read, write, and manipulate analysis result files.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import  sys
import  csv
import  time
import  os

from    itertools           import groupby
from    fileutils           import autofile


class SectionInUseError(RuntimeError): pass
class SectionNotInUseError(RuntimeError): pass
class NonUniqueSectionHeadingError(RuntimeError): pass


class GroupKey(object):
  '''
  A grouping section for reading sections from a file.
  '''
  def __init__(self):
    self.key='data'
    self.section_num=0

  def __call__(self, row):
    if (row and row[0].startswith('[')
            and row[0].endswith(']')):
      self.key = row[0][1:-1]
      self.section_num += 1

    return self.section_num,self.key


class SectionWriter(object):
  '''
  The sole purpose of SectionWriter is to govern section writing to a file
  to ensure that two separate pieces of code do not attempt to write to the
  same section.
  '''
  def __init__(self, filename, dialect='tsv'):
    self.writer = csv.writer(autofile(filename, 'w'),dialect=dialect)
    self.header = None

  def checkout(self,header):
    if self.header is not None:
      raise SectionInUseError, "Section '%s' is already checked out" % self.header

    self.header = header
    self.writer.writerow(['[%s]' % header])
    return self.writer

  def checkin(self,header):
    if self.header is None:
      raise SectionNotInUseError, "Invalid check-in attempted on section '%s'" % header
    self.header = None


def save_section(swriter, heading, rows):
  '''
  Writes a section to file.

  @param swriter: A SectionWriter object for checking out a section
  @type swriter: SectionWriter
  @param heading: A section heading to checkout from the SectionWriter
  @type heading: string
  @param rows: A list of data for the specified section heading
  @type rows: sequence

  >>> import StringIO
  >>> sio = StringIO.StringIO()
  >>> sw = SectionWriter(sio)
  >>> rows = [['analysis','completion'], ['version','66.69']]
  >>> save_section(sw,'head',rows)
  >>> sio.getvalue().split('\\r\\n')
  ['[head]', 'analysis\\tcompletion', 'version\\t66.69', '']
  '''
  writer = swriter.checkout(heading)
  try:
    writer.writerows(rows)
  finally:
    swriter.checkin(heading)


def index_sections(sections):
  '''
  Creates a dictionary of materialized sections with the key being the
  heading of each section.

  @param sections: A list of tuples containing the section
  heading and section
  @type sections: sequence

  >>> sections = [('head',      [['analysis','completion'], ['version','66.69']]),
  ...             ('samples',  [['sb12345'], ['sb12345']]),
  ...             ('loci',     [['a-00023'], ['a-00082']])]
  >>> sections = filter_sections(sections,[],exclude=True)
  >>> index_sections(sections)
  {'head': [[['analysis', 'completion'], ['version', '66.69']]], 'loci': [[['a-00023'], ['a-00082']]], 'samples': [[['sb12345'], ['sb12345']]]}
  '''
  idx = {}
  for heading,section in sections:
    idx.setdefault(heading,[]).append(list(section))
  return idx


def materialize_sections(sections):
  '''
  Materializes each section within a list sections

  @param sections: A list of heading and section tuples
  @type sections: sequence

  >>> sections = [('head',      [['analysis','completion'], ['version','66.69']]),
  ...             ('samples',  [['sb12345'], ['sb12345']]),
  ...             ('loci',     [['a-00023'], ['a-00082']])]
  >>> list(materialize_sections(sections))
  [('head', [['analysis', 'completion'], ['version', '66.69']]), ('samples', [['sb12345'], ['sb12345']]), ('loci', [['a-00023'], ['a-00082']])]
  '''
  return ((heading,list(section)) for heading,section in sections)


def filter_sections(sections,sectionset,exclude=False):
  '''
  Includes (default) or excludes sections indicated by a list of sections.

  @param sections: A list of tuples, each being a section head and its section
  @type sections: sequence
  @param sectionset: A list or set of section headings
  @type sectionset: sequence
  @param exclude: A flag to exclude the provide sections instead of including
  @type exclude: boolean

  >>> sections = [('head',     [['analysis','completion'], ['version','66.69']]),
  ...             ('samples',  [['sb12345'], ['sb12345']]),
  ...             ('loci',     [['a-00023'], ['a-00082']])]
  >>> sectionset = ['head','loci']
  >>> sections = filter_sections(sections,sectionset)
  >>> list(materialize_sections(sections))
  [('head', [['analysis', 'completion'], ['version', '66.69']]), ('loci', [['a-00023'], ['a-00082']])]
  '''
  if exclude:
    for heading,rows in sections:
      if heading not in sectionset:
        yield heading,rows
  else:
    for heading,rows in sections:
      if heading in sectionset:
        yield heading,rows


def read_sections(data):
  '''
  A generator for section headings and their content.  Contents are valid
  only until the next section is read.

  @param data: a sequence of rows, generally from a file
  @type data: sequence of sequences
  @return: An iterable of section heading and section content

  >>> data = [['[head]'],    ['analysis','completion'], ['version','66.69'],
  ...         ['[samples]'], ['sb12345'], ['sb12345'],
  ...         ['[samples]'], ['sb12315'], ['sb12311'],
  ...         ['[loci]'],    ['a-00023'], ['a-00082']]
  >>> for k,g in read_sections(data):
  ...   print k,list(g)
  head [['analysis', 'completion'], ['version', '66.69']]
  samples [['sb12345'], ['sb12345']]
  samples [['sb12315'], ['sb12311']]
  loci [['a-00023'], ['a-00082']]
  '''
  grouper = GroupKey()
  for (section_num,heading),row in groupby(data, grouper):
    row.next()
    yield heading,row


def save_metadata_section(swriter, **kwargs):
  '''
  A general method for writing a section of meta data to a file.
  Default section data are:
    generated:  current time provided by localtime()
    cwd:        the current working directory from os.getcwd()
    uname:      platform identity from os.uname()
    user:       users login name from os.getlogin()
    argv:       command line arguments used from sys.argv

  @param swriter: A SectionWriter object for checking out a section
  @type swriter: SectionWriter

  >>> import StringIO
  >>> sio = StringIO.StringIO()
  >>> sw = SectionWriter(sio)
  >>> save_metadata_section(sw, analysis='completion', project='GR-0460')
  >>> sorted(sio.getvalue().split('\\r\\n')) # doctest: +ELLIPSIS
  ['', '[header]', 'analysis\\tcompletion', "argv\\t['...']", 'cwd\\t...', 'generated\\t...', 'project\\tGR-0460', 'uname\\t...', 'user\\t...']
  '''
  header = kwargs.pop('section','header')

  kwargs.setdefault('generated', time.asctime())
  kwargs.setdefault('cwd',       os.getcwd())
  kwargs.setdefault('uname',     ' '.join(os.uname()))
  kwargs.setdefault('user',      os.getlogin())
  kwargs.setdefault('argv',      str(sys.argv))

  save_section(swriter, header, kwargs.iteritems())


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
