# -*- coding: utf-8 -*-
'''
File:          tabreports.py

Authors:       Brian Staats (staatsb@mail.nih.gov), Nick Xiao (xiaon@mail.nih.gov)

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
from    glu.lib.utils    import autofile,percent,pick
from    glu.lib.sections import read_sections,filter_sections,index_sections, materialize_sections
from    reports          import retrieve_section

COMPLETION_SECTIONS=['summary', 'data']


def completion_percent(data, summary):
  data      = islice(data,2,None)
  summary   = dict(summary)
  total     = summary['total']

  for loc,num in data:
    yield (loc,percent(int(num),int(total)))


def tabulate_sample1(completion, sampledef):
  header,sampledef = sampledef[0],sampledef[1:]

  comp = None
  if completion:
    loc_summ  = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Loci')
    sam_data  = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Sample')
    comp      = dict(completion_percent(sam_data, loc_summ))

  data = [header[1:] + ['COMPLETION']]

  data_dict = {}

  for row in sampledef:
    (sample,ind),rest = row[:2],row[2:]

    if ind in data_dict:
      if rest[6] != "TRUE":
        continue

    if rest[6] == "QC":
      continue

    if comp:
      comp_rate = '%8.2f'% comp[ind] if ind in comp else ''
      row.append(comp_rate)

    data_dict[ind] = row[1:]

  for row in data_dict.values():
    data.append(row)

  return data


def tabulate_sample(completion, sampledef):
  comp = None
  if completion:
    loc_summ  = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Loci')
    sam_data  = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Sample')
    comp      = dict(completion_percent(sam_data, loc_summ))

  header = ['SID','ASSAYABLE','IDENTIFILERREMARKS']
  def_header,sampledef = sampledef[0],sampledef[1:]
  pos = [def_header.index(h) for h in header]
  pos_assayable = header.index('ASSAYABLE')
  pos_sid = header.index('SID')

  data_dict = {}

  for srow in sampledef:
    srow = srow + [''] if len(srow) != pos[-1] + 1 else srow
    row = pick(srow, pos)

    if row[pos_assayable] == 'QC':
      continue

    sid   = row[pos_sid]
    if sid in data_dict:
      if data_dict.get(sid)[pos_assayable] == 'TRUE':
        continue
      elif row[pos_assayable] == "TRUE":
        pass
      elif data_dict.get(sid)[pos_assayable] == row[pos_assayable]:
        continue
      else:
        row[pos_assayable] == 'Multiple Failure'

    if comp:
      if sid not in comp and row[pos_assayable] == "TRUE":
        raise ValueError, 'Stopped: Sample %s with TRUE status but NO completion data' % sid
      comp_rate = '%8.2f'% comp[sid] if sid in comp else ''
      row.append(comp_rate)

    data_dict[sid] = row

  data = [header + ['COMPLETION']]
  for row in data_dict.values():
    data.append(row)

  return data


def tabulate_locus(completion, assaydef, lindex="dbsnp", rstrand='gene'):
  sam_summ = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Sample')
  loc_data = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Loci')
  comp = dict(completion_percent(loc_data, sam_summ))

  header_map = {'Strand':'Gene Strand', 'dbsnpid':'dbSNP ID', 'Chr':'Chromosome', 'Alias':'CGF Alias'}
  header_include = ('PanelID','InternalAssayid','ExternalAssayID','dbsnpid','Alias','Contig Pos',
                    'Contig ID','Gene','Strand','Chr','Genomic Systematic','Protein Systematic')

  header,assaydef = assaydef[0],assaydef[1:]

  pos = [header.index(h) for h in header_include]
  new_header = [header_map.get(h,h) for h in header_include]
  new_header.insert(5, 'Reporting Strand')
  data = [new_header+['COMPLETION']]

  for lrow in assaydef:
    assay,externalassay,dbsnp, alias= lrow[1:5]
    row = pick(lrow,pos)

    if not rstrand:
      rstrand = 'gene'

    if rstrand == 'gene':
      row.insert(new_header.index('Reporting Strand'), row[new_header.index('Gene Strand')-1])
    elif rstrand == 'contig':
      row.insert(new_header.index('Reporting Strand'), 'F')
    elif rstrand == 'custom':
      row.insert(new_header.index('Reporting Strand'), ' ')
    else:
      raise ValueError, 'Reporting strand specified is not valid: %s ' %rstrand

    if lindex == "dbsnp":
      if not dbsnp:
        print >> sys.stderr, 'Warning: No dbSNP id for the assay: %s ' % assay
      comp_rate = '%8.2f'% comp[dbsnp] if dbsnp in comp else ''
    elif lindex == "alias":
      comp_rate = '%8.2f'% comp[alias] if alias in comp else ''
    elif lindex == "external":
      comp_rate = '%8.2f'% comp[externalassay] if externalassay in comp else ''
    else:
      raise ValueError, 'Assay/locus identifier specified is not valid: %s ' %lindex

    row.append(comp_rate)
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
  parser.add_option('-a', '--assaydef', dest='assaydef', metavar='FILE',
                    help='Assay definition file (assay.def)')
  parser.add_option('-m', '--sampledef', dest='sampledef', metavar='FILE',
                    help='Sample definition file (sample.def)')
  parser.add_option('-i', '--locusindex', dest='locusindex', metavar='N', type='str', default='dbsnp',
                    help='Locus ID used for in data reporting (dbsnp/alias/external)')
  parser.add_option('-r', '--reportingstrand', dest='rstrand', metavar='STRAND', type='str', default='gene',
                    help='The strand that the genotyping data is reported on gene(default),or contig')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not any([options.completion, options.assaydef, options.sampledef]):
    print >> sys.stderr, 'At least one of the input files are missing!'
    parser.print_help()
    return

  if not all([options.locusreport, options.samplereport]):
    print >> sys.stderr, 'At least one of the output files need to be specified!'
    parser.print_help()
    return

  assaydef = sampledef = completion = None

  completion    =  list(csv.reader(autofile(options.completion), dialect='excel-tab'))
  completion    =  index_sections(filter_sections(read_sections(completion),COMPLETION_SECTIONS))

  assaydef   = list(csv.reader(autofile(options.assaydef), dialect='excel-tab'))
  sampledef   = list(csv.reader(autofile(options.sampledef), dialect='excel-tab'))

  if options.locusreport:
    locusid = options.locusindex
    locusreport   = csv.writer(autofile(options.locusreport, 'w'),dialect='excel-tab')
    tab           = tabulate_locus(completion, assaydef, locusid, options.rstrand)
    locusreport.writerows(tab)
  if options.samplereport:
    samplereport  = csv.writer(autofile(options.samplereport, 'w'),dialect='excel-tab')
    tab           = tabulate_sample(completion, sampledef)
    samplereport.writerows(tab)


if __name__ == '__main__':
  main()
