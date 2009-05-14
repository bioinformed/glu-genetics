# -*- coding: utf-8 -*-

__abstract__  = 'parser for region files'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys
import csv
import unittest

from   StringIO             import StringIO
from   itertools            import izip

from   glu.lib.sections     import SectionWriter,read_sections,save_section
from   glu.lib.fileutils    import autofile


class Regions(object):
  '''
  Python implementation of a region representation object
  which stores the name and a set of samples or loci for each region.
  '''
  def __init__(self,regions=[]):
    '''
    Regions object constructor

    @param  regions: region that was passed in
    @type   regions: regions object
    '''
    self.regions = {}
    for name,samples,loci in regions:
      self.add_region(name,samples,loci)

  def add_region(self,name,samples,loci):
    '''
    Add new region into the current regions object with the supplied information

    @param     name: region name
    @type      name: str
    @param  samples: list of sample ids
    @type   samples: list
    @param     loci: list of locus ids
    @type      loci: list
    '''
    self.regions[name] = samples,loci

  def __contains__(self,key):
    '''
    Return whether or not the sample-locus pair that was passed in exists in the current regions object

    @param    key: sampel id and locus id
    @type     key: tuple
    @return      : flag indicating if the sample-locus exists in the current regions object
    @rtype       : bool
    '''
    sample,locus = key
    for samples,loci in self.regions.itervalues():
      if sample in samples and locus in loci:
        return True
    return False

  def __len__(self):
    '''
    Return the count of regions in the current regions object
    '''
    return len(self.regions)

  def __iter__(self):
    '''
    Iterate through the current regions object for each ssample and locus in each region
    '''
    return ((name,samples,loci) for name,(samples,loci) in self.regions.iteritems())

  def __getitem__(self,name):
    '''
    Return all samples and loci in the specified region

    @param   name: region name
    @type    name: str
    @return      : sample and locus sets for the specified region
    @rtype       : list of set
    '''
    return self.regions[name]


def save_regions(file_or_name,regions):
  '''
  Save the regions into a file using SectionWriter.

  @param file_or_name: file name or file object
  @type  file_or_name: str or file object
  @param      regions: regions object that was passed in
  @type       regions: Regions object
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
    raise RegionParserError('Format Error: only one row and one column is allowed for specifying the region name')
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
      raise RegionParserError('Format Error: only one sample or loci per row')
    aset.add(row[0])
  return aset


def build_error_msg(dfa,state,heading):
  expect = [cond for (s,cond),new_s in dfa.iteritems() if state == s and new_s != 'error']
  if not expect:
    expect = ['samples','loci','region']
  return 'Parsing Error: expect %s but got %s' % (' or '.join(expect),heading)


def new_state(dfa,state,heading):
  '''
  Return a valid state, otherwise raise the error

  @param      dfa: FSA look up table
  @type       dfa: dict
  @param    state: current FSA state
  @type     state: str
  @param  heading: section name
  @type   heading: str
  @return        : next FSA state
  @rtype         : str

  >>> dfa = build_transition_table()
  >>> state = 'start'
  >>> heading = 'samples'
  >>> new_state(dfa,state,heading)
  'samples'
  >>> state = 'locisample'
  >>> heading = 'region'
  >>> new_state(dfa,state,heading)
  'region'
  '''
  if heading=='data':
    raise RegionParserError('Parsing Error: The section file has no header')

  new_state = dfa.get( (state,heading),'error' )

  if new_state == 'error':
    raise RegionParserError(build_error_msg(dfa,state,heading))

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
    raise RegionParserError(build_error_msg(dfa,state,heading))


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


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
  unittest.main()
