# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      dupcheck for web-app

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


def run_identifiler_dupcheck(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  prep_ident_data(path_dir)

  id_data        = data_dir + '/' +proj_name+'_Identifiler.trip'
  sample_def     = data_dir + '/' +proj_name+'_Sample.txt'
  new_sample_def = work_dir + '/' + 'new_sample_def'
  expdup_id      = work_dir + '/' + 'expdup_id'
  _ids = os.path.abspath(path_dir['script_dir2']) + '/ids.py'
  f = os.popen('python %s -v -s 1 -S 0 %s %s --sampledef=%s ' % (_ids,id_data,new_sample_def,sample_def)).close()
  if f:
    raise ValueError,'Failed in creating new sample definition file!'
  f = os.popen( 'python %s --cols=1,2 %s %s ' % (_ids,new_sample_def,expdup_id) ).close()
  if f:
    raise ValueError,'Failed in creating expected dup for identifiler file!'

  # run dupcheck on identifiler
  input_file = work_dir + '/' + proj_name + '_Identifiler.sdat'
  output_file= res_dir  + '/' + proj_name + '_dupcheck_Identifiler.out'
  _dupcheck2 = os.path.abspath(path_dir['script_dir']) + '/dupcheck2.py'

  dups       = work_dir + '/' + 'dups.txt'
  f = os.popen('python %s -f sdat %s  -T %s -m %s -o %s -g generic -e %s -d %s --phenofile=%s --phenos=%s ' \
            % (_dupcheck2,input_file,path_dir['Ident_threshold'],path_dir['Ident_mincount'], \
               output_file,expdup_id,dups,new_sample_def,'sex') ).close()
  if f:
    raise ValueError,'Failed in creating Identifiler dupcheck.out file!'


def dupcheck_withExp(path_dir,genorepr):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  _dupcheck2 = os.path.abspath(path_dir['script_dir']) + '/dupcheck2.py'
  if 'yes' in path_dir['Geno_dupcheck_lite']:
    _dupcheck2 = os.path.abspath(path_dir['script_dir2']) + '/dupcheck2_lite.py'

  _ids = os.path.abspath(path_dir['script_dir2']) + '/ids.py'

  duplicates =  work_dir + '/duplicates.txt'
  input_file =  work_dir + '/' + proj_name + '.sdat'

  if path_dir['Geno_expdup_type']=='ident_dups':
    output_file= res_dir  + '/' + proj_name + '_dupcheck_identExp.out'
    dups       = work_dir + '/' + 'dups.txt'
    if not os.path.exists(dups):
      if 'Ident_threshold' in path_dir and 'Ident_mincount' in path_dir:
        run_identifiler_dupcheck(path_dir)
      else:
        print >> sys.stderr,'Stopped: identifiler dupcheck output NOT found! '
        return
    f = os.popen('python %s -s 0 -S 0 -v %s --dupfile=%s %s ' % (_ids,input_file,dups,duplicates) ).close()
    if f:
      raise ValueError,'Failed in creating expected dups (from identifiler data) for genotype data!'
  elif path_dir['Geno_expdup_type']=='pi_dups':
    output_file = res_dir  + '/' + proj_name + '_dupcheck_PIExp.out'
    expdup_pi   = work_dir + '/expdup_pi'
    sample_def  = data_dir + '/' + proj_name + '_Sample.txt'

    new_sample_def = work_dir + '/new_sample_def'
    f = os.popen('python %s -v -s 1 -S 0 %s %s --sampledef=%s ' % (_ids,input_file,new_sample_def,sample_def) ).close()
    if f:
      raise ValueError,'Failed in creating new sample definition file!'
    f = os.popen( 'python %s --cols=1,4 %s %s ' % (_ids,new_sample_def,duplicates) ).close()
    if f:
      raise ValueError,'Failed in creating PI expected dups for genotype dupcheck !'
  else:
    raise NotImplementedError,'The expected duplicate type: %s is not supported yet' % path_dir['Geno_expdup_type']


  f = os.popen('python %s -f sdat %s -o %s -T %s -m %s -e %s -g %s '
            % (_dupcheck2,input_file,output_file,path_dir['Geno_threshold'],path_dir['Geno_mincount'], \
               duplicates,genorepr) ).close()
  if f:
    raise ValueError,'Failed in creating dupcheck.out file!'


def process(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  FILES = [proj_name+'_dupcheck_PIExp.out',proj_name+'_dupcheck_identExp.out',proj_name+'_dupcheck.out']
  FILES += ['duplicates.txt','expdup_pi']
  clean_dir(res_dir,FILES)

  pf = path_dir['platform']
  to_process = 'yes' in path_dir['preprocess'].lower()
  if pf == 'gg':
    if to_process:
      prep_gg_data(path_dir)
    genorepr = 'snp_acgt'
  elif pf in ('taqman','snplex'):
    genorepr = 'generic'
    if to_process:
      prep_taq_data(path_dir,genorepr)
  else:
    raise ValueError,'Platform should be "gg","taqman",or "snplex" !'

  input_file = work_dir + '/' + proj_name + '.sdat'

  _dupcheck2 = os.path.abspath(path_dir['script_dir']) + '/dupcheck2.py'

  #-- case I: dupcheck without expected dups
  if 'no' in path_dir['Geno_with_expdup'].lower():
    output_file= res_dir  + '/' + proj_name + '_dupcheck.out'
    f = os.popen('python %s -f sdat %s  -T %s -m %s -o %s -g %s' \
        % (_dupcheck2,input_file,path_dir['Geno_threshold'],path_dir['Geno_mincount'],output_file,genorepr) ).close()
    if f:
      raise ValueError,'Failed in creating dupcheck.out file (without expected dups)!'
  else:
    #-- case II: dupcheck with expected dups
    dupcheck_withExp(path_dir,genorepr)


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
