# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      concordance check using a config file as input

Requires:      Python 2.4 or higher

Revision:
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys,os
import csv

from   biozilla.utils         import autofile,hyphen
from   biozilla.genoarray     import snp_acgt,get_genorepr
from   biozilla.genodata      import load_map
from   biozilla.genomerge     import VoteMerger
from   app_utils              import *
from   operator               import itemgetter

def remove_QC(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name =path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  sampledef = data_dir + '/' + proj_name + '_Sample.txt'
  fi = csv.reader(autofile(sampledef,'r'),dialect='excel-tab')
  header = fi.next()
  sb  = header.index('LIMS')
  pos = header.index('ASSAYABLE')
  status = {}
  for r in fi:
    if  len(r) > pos:
      status[r[sb]] = r[pos]

  sdat   = work_dir + '/' + proj_name + '.sdat'
  n_sdat = work_dir + '/' + proj_name + '_n.sdat'
  fi = csv.reader(autofile(sdat,'r'),dialect='excel-tab')
  fo = csv.writer(autofile(n_sdat,'w'),dialect='excel-tab')
  fo.writerow(fi.next())
  for r in fi:
    new_id = mapper(r[0])
    if status[new_id] == 'TRUE':
      fo.writerow(r)


def process(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name =path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  # clean old files
  FILES = [proj_name+'_completion.out',proj_name + '.sdat']
  clean_dir(res_dir,FILES)

  pf = path_dir['platform']
  if pf == 'gg':
    genorepr = 'snp_acgt'
    prep_gg_data(path_dir)
    input_file = work_dir + '/' + proj_name + '.sdat'
  elif pf in ('taqman','snplex'):
    genorepr = 'generic'
    prep_taq_data(path_dir,outgenorepr=genorepr)
    remove_QC(path_dir)
    input_file = work_dir + '/' + proj_name + '_n.sdat'
  else:
    raise ValueError,'Platform should be "gg","taqman",or "snplex" !'

  # run completion
  output_file= res_dir  + '/' + proj_name + '_completion.out'

  completion= os.path.abspath(path_dir['script_dir']) + '/completion2.py'
  f = os.popen('python %s -f sdat %s -o %s -r %s' \
            % (completion,input_file,output_file,genorepr) ).close()
  if f:
    raise ValueError,'Failed in completion!'


def option_parser():
  import optparse

  usage = 'usage: %prog configFile'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    print >> sys.stderr, 'Stopped: require a config file!'
    return

  path_dir = load_map(args[0])
  path_dir = prep_dir(path_dir)
  process(path_dir)


if __name__ == '__main__':
  main()
