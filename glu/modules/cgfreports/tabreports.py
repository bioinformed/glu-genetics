# -*- coding: utf-8 -*-
'''
File:          tabreports.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-09-15

Abstract:      A report module for generating tabular reports.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import  sys
import  csv

from    itertools        import islice,chain
from    operator         import itemgetter

from    glu.lib.genodata import load_map_file
from    glu.lib.utils    import autofile,percent
from    glu.lib.sections import read_sections,filter_sections,index_sections, materialize_sections


COMPLETION_SECTIONS=['summary', 'data']


def retrieve_section(sections, heading, section_type=None):
  section = sections.get(heading)

  if not section_type and len(section) >1:
    raise NoneUniqueSectionHeadingError, 'More than one section with the heading: ' + heading

  if not section_type:
    return section[0]

  for items in section:
    for item in items:
      if section_type in set(item):
        return items


def completion_percent(data, summary):
  data      = islice(data,2,None)
  summary   = dict(summary)
  total     = summary['total']

  for loc,num in data:
    yield (loc,percent(int(num),int(total)))


def tabulate_sample(completion, concordance, sampledef):
  header,sampledef = sampledef[0],sampledef[1:]

  comp = conc = None
  if completion:
    loc_summ  = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Loci')
    sam_data  = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Sample')
    comp      = dict(completion_percent(sam_data, loc_summ))

  if concordance:
    raise NotImplementedError,'Tabular reports do not yet support concordance data'

  data = [header + ['COMPLETION','CONCORDANCE']]

  for row in sampledef:
    # use either sample or ind until we standardize input and unique ids
    (sample,ind),rest = row[:2],row[2:]

    if comp:
      row.append(comp.get(sample) or comp.get(ind) or '')
    if conc:
      row.append(conc.get(sample) or conc.get(ind) or '')

    data.append(row)

  return data


def tabulate_locus(completion, concordance, hwpdata, assaydef):
  header,assaydef = assaydef[0],assaydef[1:]

  comp = conc = hwp = None
  if completion:
    sam_summ = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Sample')
    loc_data = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Loci')

    comp = dict(completion_percent(loc_data, sam_summ))

  if concordance:
    raise NotImplementedError,'Tabular reports do not yet support concordance data'

  if hwpdata:
    raise NotImplementedError,'Tabular reports do not yet support hwp data'

  data = [header+['COMPLETION','CONCORDANCE','HWP_P_VALUE','X2_ON_GENOTYPES','P_VALUE_GENOTYPES','X2_ON_ALLELES','P_VALUE_ALLELES','OR']]

  for row in assaydef:
    # use either assay or dbsnp until we standardize input and unique ids
    (assay,externalassay,dbsnp),rest = row[:3],row[3:]

    if comp:
      row.append(comp.get(assay) or comp.get(dbsnp) or '')
    if conc:
      row.append(conc.get(assay) or conc.get(dbsnp) or '')
    if hwp:
      row.append(hwp.get(assay) or hwp.get(dbsnp) or '')

    data.append(row)

  return data


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-l', '--locusreport', dest='locusreport', metavar='N', type='str',
                    help='File name to be given the locus tabulated report')
  parser.add_option('-s', '--samplereport', dest='samplereport', metavar='N', type='str',
                    help='File name to be given the sample tabulated report')

  parser.add_option('-c', '--completion', dest='completion', metavar='FILE',
                    help='Completion results (completion.dat) data file')
  parser.add_option('-d', '--dupcheck', dest='dupcheck', metavar='FILE',
                    help='Dupcheck results (dupcheck.dat) data file')
  parser.add_option('-n', '--concordance', dest='concordance', metavar='FILE',
                    help='Concordance results (concordance.dat) data file')
  parser.add_option('-w', '--hwp', dest='hwp', metavar='FILE',
                    help='HWP results (hwp.dat) data file')

  parser.add_option('-p', '--projectdef', dest='projectdef', metavar='FILE',
                    help='Project definition file (project.def)')
  parser.add_option('-a', '--assaydef', dest='assaydef', metavar='FILE',
                    help='Assay definition file (assay.def)')
  parser.add_option('-m', '--sampledef', dest='sampledef', metavar='FILE',
                    help='Sample definition file (sample.def)')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  assaydef = sampledef = completion = dupcheck = concordance = hwp = None
  if options.completion:
    completion    =  list(csv.reader(autofile(options.completion), dialect='excel-tab'))
    completion    =  index_sections(filter_sections(read_sections(completion),COMPLETION_SECTIONS))
  if options.dupcheck:
    dupcheck      =  list(csv.reader(autofile(options.dupcheck), dialect='excel-tab'))
  if options.concordance:
    concordance   =  list(csv.reader(autofile(options.concordance), dialect='excel-tab'))
  if options.hwp:
    hwp           =  list(csv.reader(autofile(options.hwp), dialect='excel-tab'))

  assaydef = sampledef = None
  if options.assaydef:
    assaydef   = csv.reader(autofile(options.assaydef), dialect='excel-tab')
    sections   = index_sections(filter_sections(read_sections(assaydef),'assays'))
    assaydef   = retrieve_section(sections, 'assays')
  if options.sampledef:
    sampledef   = csv.reader(autofile(options.sampledef), dialect='excel-tab')
    sections    = index_sections(filter_sections(read_sections(sampledef), 'samples'))
    sampledef   = retrieve_section(sections, 'samples')

  if options.locusreport:
    locusreport   = csv.writer(autofile(options.locusreport, 'w'),dialect='excel-tab')
    tab           = tabulate_locus(completion, concordance, hwp, assaydef)
    locusreport.writerows(tab)
  if options.samplereport:
    samplereport  = csv.writer(autofile(options.samplereport, 'w'),dialect='excel-tab')
    tab           = tabulate_sample(completion, concordance, sampledef)
    samplereport.writerows(tab)


if __name__ == '__main__':
  main()