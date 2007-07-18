# -*- coding: utf-8 -*-
'''
File:          regionparser.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)

Created:       Oct 25, 2006

Abstract:      This library module can be called to parse the region section file

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv
import unittest

from StringIO             import StringIO
from itertools            import chain,izip

from sections             import SectionWriter,read_sections,save_section
from fileutils            import autofile


class Regions(object):
  '''
  Store the name and a set of samples or loci for each region.
  '''
  def __init__(self,regions=[]):
    self.regions = {}
    for name,samples,loci in regions:
      self.add_region(name,samples,loci)

  def add_region(self,name,samples,loci):
    self.regions[name] = samples,loci

  def __contains__(self,key):
    sample,locus = key
    for samples,loci in self.regions.itervalues():
      if sample in samples and locus in loci:
        return True
    return False

  def __len__(self):
    return len(self.regions)

  def __iter__(self):
    return ((name,samples,loci) for name,(samples,loci) in self.regions.iteritems())

  def __getitem__(self,name):
    return self.regions[name]


def save_regions(file_or_name,regions):
  '''
  Save the regions into a file using SectionWriter.
  '''
  sw = SectionWriter(file_or_name)
  for name,samples,loci in regions:
    save_section(sw,'region', [[name]])
    save_section(sw,'samples',[[sample] for sample in samples])
    save_section(sw,'loci',   [[locus]  for locus  in loci])


def build_transition_table():
  '''
  Build an FSA look up table for parsing the region section file
  The region section file format should have two invariants:
  1. The section header [region] is optional for the first region but mandatory for the following regions
  2. In each region section, there should be one[only] sample sub-section and one[only] loci sub-section in any order
  '''
  conditions  = ['samples','region','loci']
  transitions = [['start',     ['samples',    'region',     'loci']],
                ['samples',    ['error',      'error',      'sampleloci']],
                ['loci',       ['locisample', 'error',      'error']],
                ['region',     ['samples',    'error',      'loci']],
                ['sampleloci', ['error',      'region',     'error']],
                ['locisample', ['error',      'region',     'error']]]

  dfa = {}
  for state,trans in transitions:
      dfa.update( ((state,cond),next) for cond,next in zip(conditions,trans) )

  return dfa


class RegionParserError(RuntimeError): pass


def parse_header(rows):
  '''
  Parse the region header.
  It should be one line and one column to simply give the name of the region.
  '''
  rows = list(rows)
  if len(rows) != 1 or len(rows[0]) != 1:
    raise RegionParserError, 'Format Error: only one row and one column is allowed for specifying the region name'
  return rows[0][0]


def parse_list(rows):
  '''
  Parse either the samples or loci subsection.
  It sould be one sample or loci per row.
  '''
  aset = set()
  rows = list(rows)
  for row in rows:
    if len(row) != 1:
      raise RegionParserError, 'Format Error: only one sample or loci per row'
    aset.add(row[0])
  return aset


def build_error_msg(dfa,state,heading):
  expect = [cond for (s,cond),new_s in dfa.iteritems() if state == s and new_s != 'error']
  if not expect:
    expect = ['samples','loci','region']
  return 'Parsing Error: expect %s but got %s' % (' or '.join(expect),heading)


def new_state(dfa,state,heading):
  if heading=='data':
    raise RegionParserError, 'Parsing Error: The section file has no header'

  new_state = dfa.get( (state,heading),'error' )

  if new_state == 'error':
    raise RegionParserError, build_error_msg(dfa,state,heading)

  return new_state


def load_regions(filename):
  '''
  Load the regions from a file.
  '''
  data = csv.reader(autofile(filename),dialect='tsv')
  state = 'start'
  endstates = ['start','sampleloci','locisample']
  rname = 'default'
  dfa = build_transition_table()
  for heading,rows in read_sections(data):
    state = new_state(dfa,state,heading)

    if heading == 'loci':
      loci = parse_list(rows)
    elif heading == 'samples':
      samples = parse_list(rows)
    elif heading == 'region':
      rname = parse_header(rows)

    if state in endstates:
      yield rname,samples,loci

  if state not in endstates:
    raise RegionParserError, build_error_msg(dfa,state,heading)


class testRegion(unittest.TestCase):
  def assertRaisesMsg(self, exception, callable, *args, **kwargs):
    exc_msg = kwargs["exc_msg"]
    try:
      list(callable(*args))
    except RegionParserError, exc:
      self.failIf(exc.args[0] != exc_msg, "An unexpected exception message is seen. Expected=%s, Actual=%s"
                                           % (exc_msg, exc.args[0]))
    except:
      exc_info = sys.exc_info()
      self.fail("An unexpected exception type: %s is raised. Expected=%s, actual=%s"
                 % (exc_info, 'RegionParserError', exc_info[0]))
    else:
      self.fail("no exception is raised")


  def test_load_regions(self):
    def eq(region1,region2):
      return sorted(region1) == sorted(region2)

    #test valid cases
    cases   = ["",
               "[region]\nregion1\n[samples]\ns1\ns2\n[loci]\nl1\nl2\n",
               "[samples]\ns1\ns2\n[loci]\nl1\nl2\n",
               "[samples]\ns1\ns2\n[loci]\nl1\n[region]\nregion2\n[loci]\nl1\nl2\n[samples]\ns1\n",
               "[region]\nregion1\n[samples]\ns1\n[loci]\nl1\nl2\n[region]\nregion2\n[loci]\nl1\n[samples]\ns1\ns2\n"]

    results = [Regions(),
               Regions([('region1',set(['s1','s2']),set(['l1','l2']))]),
               Regions([('default',set(['s1','s2']),set(['l1','l2']))]),
               Regions([('default',set(['s1','s2']),set(['l1'])),('region2',set(['s1']),set(['l1','l2']))]),
               Regions([('region1',set(['s1']),set(['l1','l2'])),('region2',set(['s1','s2']),set(['l1']))])]

    for case,result in izip(cases,results):
      regions = Regions()
      for region in load_regions(StringIO(case)):
        regions.add_region(*region)
      self.assertTrue(eq(regions,result))

    #test invalid cases
    cases      = ["[samples]\ns1\n[samples]\ns1\n",
                  "[samples]\ns1\n[region]\nregion2\n",
                  "[loci]\nl1\n[loci]\nl1\n",
                  "[loci]\nl1\n[region]\nregion2\n",
                  "[region]\nregion1\n[region]\nregion2\n",
                  "[samples]\ns1\n[loci]\nl1\n[samples]\ns1\n",
                  "[samples]\ns1\n[loci]\nl1\n[loci]\nl1\n",
                  "[loci]\nl1\n[samples]\ns1\n[samples]\ns1\n",
                  "[loci]\nl1\n[samples]\ns1\n[loci]\nl1\n",
                  "abc\n",
                  "[abc]\n"]

    exceptions = ['Parsing Error: expect loci but got samples',
                  'Parsing Error: expect loci but got region',
                  'Parsing Error: expect samples but got loci',
                  'Parsing Error: expect samples but got region',
                  'Parsing Error: expect samples or loci but got region',
                  'Parsing Error: expect region but got samples',
                  'Parsing Error: expect region but got loci',
                  'Parsing Error: expect region but got samples',
                  'Parsing Error: expect region but got loci',
                  'Parsing Error: The section file has no header',
                  'Parsing Error: expect region or loci or samples but got abc']

    for case,exception in izip(cases,exceptions):
      self.assertRaisesMsg(RegionParserError, load_regions, StringIO(case), **{'exc_msg':exception})


if __name__ == '__main__':
  unittest.main()
