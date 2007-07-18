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

def process(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name = path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  # clean old files
  FILES = [proj_name+'_dupcheck_Identifiler.out',proj_name+'_Identifiler.sdat']
  FILES += ['new_sample_def','expdup_id']
  clean_dir(res_dir,FILES)

  prep_ident_data(path_dir)

  input_file = work_dir + '/' + proj_name + '_Identifiler.sdat'
  output_file= res_dir  + '/' + proj_name + '_dupcheck_Identifiler.out'
  _dupcheck2 = os.path.abspath(path_dir['script_dir']) + '/dupcheck2.py'
  dups       = work_dir + '/' + 'dups.txt'

  # dupcheck on identifiler without expected dups
  if 'no' in path_dir['Ident_with_expdup'].lower():
    f = os.popen('python %s -f sdat %s  -T %d -m %d -o %s -g generic -d %s ' \
              % (_dupcheck2,input_file,int(path_dir['Ident_threshold']),int(path_dir['Ident_mincount']), \
               output_file,dups) ).close()
    if f:
      raise ValueError,'Failed in creating dupcheck.out file (without expected dups)!'

    return None

  # create exp_dups for identifiler data
  id_data        = data_dir + '/' +proj_name+'_Identifiler.trip'
  sample_def     = data_dir + '/' +proj_name+'_Sample.txt'
  new_sample_def = work_dir + '/' + 'new_sample_def'
  expdup_id      = work_dir + '/' + 'expdup_id'
  _ids = os.path.abspath(path_dir['script_dir2']) + '/ids.py'
  pheno_file   = ''
  pheno_var    = 'sid'

  is_pre_geno = path_dir['pre_geno_qc'].lower()=='yes'
  if is_pre_geno:
    if 'pheno_file' in path_dir:
      pheno_file = os.path.abspath(path_dir['pheno_file'])

  if os.path.exists(pheno_file) and is_pre_geno:
    pheno_var  = 'pid'
    output_file= res_dir  + '/' + proj_name + '_dupcheck_Identifiler_PID.out'
    f = os.popen('python %s -v -s 1 -S 0 %s %s --sampledef=%s --phenofile=%s' \
              % (_ids,id_data,new_sample_def,sample_def,pheno_file)).close()
    if f:
      raise ValueError,'Failed in creating new sample definition file!'

    # gender check
    _transform2 = os.path.abspath(path_dir['script_dir']) + '/transform2.py'
    amel_include= os.path.abspath(path_dir['meta_dir'])   + '/include_AMEL'
    amel_sdat   = work_dir + '/' + proj_name + '_AMEL.sdat'
    f = os.popen( 'python %s --cols=1,2 %s %s ' % (_ids,new_sample_def,expdup_id) ).close()
    f = os.popen('python %s -f sdat -g generic %s -F sdat -o %s -u %s -c' \
               %(_transform2,input_file,amel_sdat,amel_include) ).close()
    if f:
      raise ValueError,'Failed in creating AMEL sdat file!'

    #--
    output_gender = res_dir  + '/' + proj_name + '_dupcheck_gender.out'
    dup_lite      = os.path.abspath(path_dir['script_dir2']) + '/dupcheck2_lite.py'
    f = os.popen('python %s -f sdat %s  -T 100 -m 1 -o %s -g generic -e %s  --phenofile=%s --phenos=%s ' \
              % (dup_lite,amel_sdat, output_gender,expdup_id,new_sample_def,pheno_var) ).close()
    if f:
      raise ValueError,'Failed in creating dupcheck_gender.out file!'

    #--
    output_gender = res_dir  + '/' + proj_name + '_diff_gender.out'
    gender_check(amel_sdat,new_sample_def,output_gender)
  else:
    print >> sys.stderr,'\n Note: Run SID-based dupcheck !\n'
    output_file= res_dir  + '/' + proj_name + '_dupcheck_Identifiler_SID.out'
    f = os.popen('python %s -v -s 1 -S 0 %s %s --sampledef=%s ' \
              % (_ids,id_data,new_sample_def,sample_def)).close()
    if f:
      raise ValueError,'Failed in creating new sample definition file!'

  f = os.popen( 'python %s --cols=1,2 %s %s ' % (_ids,new_sample_def,expdup_id) ).close()
  if f:
    raise ValueError,'Failed in creating expected dup for identifiler file!'

  # run dupcheck on identifiler with expected dups
  if 'yes' in path_dir['Ident_dupcheck_lite'].lower():
    _dupcheck2 = os.path.abspath(path_dir['script_dir2']) + '/dupcheck2_lite.py'

  f = os.popen('python %s -f sdat %s  -T %s -m %s -o %s -g generic -e %s -d %s --phenofile=%s --phenos=%s ' \
            % (_dupcheck2,input_file,path_dir['Ident_threshold'],path_dir['Ident_mincount'], \
               output_file,expdup_id,dups,new_sample_def,pheno_var) ).close()
  if f:
    raise ValueError,'Failed in creating dupcheck.out file!'


def gender_check(amel_sdat,new_sdef,output_file):
  iid2pheno   = {}
  f      = csv.reader(autofile(new_sdef,'r'),dialect='excel-tab')
  header = f.next()
  header = [t.lower() for t in header]
  pos    = [header.index(t) for t in ('pid','sex')]
  for row in f:
    iid2pheno[row[0]] = pick(row,pos)

  fout   = csv.writer(autofile(output_file,'w'),dialect='excel-tab')
  fin    = csv.reader(autofile(amel_sdat,'r'),dialect='excel-tab')
  fin.next()
  fout.writerow(['SB#','SEX_IDENT','PID','SEX_PI'])
  outrows = []
  for row in fin:
    sex_ident = '-'
    if row[1] == 'X/X':
      sex_ident = 'F'
    elif row[1] in ('X/Y','Y/X'):
      sex_ident = 'M'
    if row[0] in iid2pheno:
      sex_pi    = iid2pheno.get(row[0])[1]
      if sex_pi and sex_pi.lower() != sex_ident.lower():
        outrows.append([row[0],sex_ident.upper()] + iid2pheno.get(row[0]))
      if not sex_pi:
        outrows.append([row[0],sex_ident.upper()] + [iid2pheno.get(row[0])[0],'-'])
    else:
      outrows.append([row[0],sex_ident.upper()] + ['-','-'])

  outrows_1  = []
  outrows_2  = []
  for r in outrows:
    if (r[1],r[3]) in (('M','F'),('F','M')):
      outrows_1.append(r)
    else:
      outrows_2.append(r)

  fout.writerows(outrows_1 + outrows_2)


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
