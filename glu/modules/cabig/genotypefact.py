# -*- coding: utf-8 -*-
'''
File:          genotypefact.py

Inputs:        scan1_merged.loci:    check for QC
               sample.def:           illumina sample id to specimen de-identified id map
               specimen.csv:         filter using specimen id by inclusion
               sample.map:           a map betweent the id that Kevin constructed and illumina sample id
               genotype_fact.csv.gz: the output data file for GENOTYPE_FACT table
               snpid.map:            a map between the rs number and snp annotation id

Requires:      glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import csv

from   itertools         import chain,islice

from   glu.lib.fileutils import autofile,load_map,load_list
from   utils             import load_rows



HEADER = ['ALLELE1','ALLELE2','QUALITY_SCORE','NORMAL_X','NORMAL_Y','RAW_X','RAW_Y','ASSAY_ID',
          'SPECIMEN_ID','SNPANNO_ID','STUDY_NAME','SCAN_ID','STATUS']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genotriplefile'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-o', '--outfile',      dest='outfile',      default='-', metavar='FILE',
                    help="the name of the output file, '-' for standard out")
  parser.add_option('-s', '--samdef',       dest='samdef',       default='-', metavar='FILE',
                    help="the name of the sample definition file, '-' for standard in")
  parser.add_option('-p', '--specimenfile', dest='specimenfile', default='-', metavar='FILE',
                    help="the name of the containing list of specimens, '-' for standard in")
  parser.add_option('-q', '--qcsnps',       dest='qcsnps',                    metavar='FILE')
  parser.add_option('-r', '--snpidmap',     dest='snpidmap',                  metavar='FILE')

  return parser


def process(genotypes,options):
  qcsnps = set(row[0] for row in load_rows(options.qcsnps))
  samspec = load_map(options.samdef)
  snpidmap = load_map(options.snpidmap)
  specimens = load_list(options.specimenfile)

  for genotype in islice(genotypes,0,None):
    snp = genotype[0]
    illsid = genotype[1]

    status = 'QC-'
    if snp in qcsnps:
      status = 'QC+'

    sbnumber = illsid.split('_')[0]
    specimen_id = samspec[sbnumber]
    if specimen_id not in specimens:
      continue

    snpanno_id = snpidmap[snp]
    study_name = 'CGEMS Prostate Cancer WGAS Phase 1'
    scan_id = 1
    results = genotype[2:]
    if results[0] == '-':
      results[0] = None
    if results[1] == '-':
      results[1] = None
    if results[3] == 'NaN':
      results[3] = None
    if results[4] == 'NaN':
      results[4] = None
    # assay_id and snpanno_id are always equal in the way we generated snp_dim and snp_assay
    assay_id = snpanno_id
    results.extend([assay_id,specimen_id,snpanno_id,study_name,scan_id,status])

    yield results


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  out = csv.writer(autofile(options.outfile,'w'))
  out.writerow(HEADER)

  genotypes = csv.reader(autofile(args[0]),dialect='excel-tab')
  results = process(genotypes,options)
  for row in results:
    out.writerow(row)


if __name__=='__main__':
  main()
