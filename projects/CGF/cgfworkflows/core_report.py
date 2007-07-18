# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      Final report delivery

Requires:      Python 2.4 or higher

Revision:
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys,os
from subprocess import Popen,PIPE
import csv

from   biozilla.utils         import autofile,hyphen
from   biozilla.genomerge     import VoteMerger
from   app_utils              import *

ASSAY_ID_MAP = dict([('dbsnp',4),('alias',5),('external',3)])

def create_report(path_dir,col,strand):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name = path_dir['proj_name']
  final_dir = os.path.abspath(path_dir['final_dir'])
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  sample_def = final_dir + '/' + proj_name + '_Sample.txt'
  assay_def  = data_dir  + '/' + proj_name + '_Assay.txt'
  trip_file  = data_dir  + '/' + proj_name + '.trip'

  _ids = os.path.abspath(path_dir['script_dir2']) + '/ids.py'

  # new sample def (with updated reason annotation)
  _sampledef = os.path.abspath(path_dir['script_dir2']) + '/sampledef_prep.py'
  reason_f   = os.path.abspath(path_dir['meta_dir'])    + '/exclude_reason_list.txt'
  rank_f     = os.path.abspath(path_dir['meta_dir'])    + '/exclude_priority_list.txt'
  n_sample_def = work_dir+ '/new_sample_def'

  f = os.popen('python %s %s %s -m %s -r %s' \
              % (_sampledef,sample_def,n_sample_def,reason_f,rank_f) ).close()
  if f:
    raise ValueError,'Failed in creating new_sample_def file!'

  # sample include
  sample_include  = work_dir+ '/sample_include'
  f  = os.popen('python %s -s 0 -S 0 -v %s --filtercol=9 --filterval=TRUE %s '    \
             % (_ids,n_sample_def,sample_include) ).close()
  if f:
    raise ValueError,'Failed in creating sample_include file!'

  # sample and assay map
  data_file    = work_dir + '/'+ proj_name+'.sdat'
  if col==3:
    data_file  = work_dir+ '/' + proj_name + '.sdat'

  to_col       = ASSAY_ID_MAP[path_dir['report_assayid']]
  long2short   = work_dir+ '/long2short_map'
  sid_map      = work_dir+ '/sid_map'
  locus_map    = work_dir+ '/locus_map'
  f   = os.popen('python %s -f 1 -t 2 %s %s'     % (_ids,n_sample_def,sid_map) ).close()
  f_s = os.popen('python %s -s 0 -S 0 -v %s %s ' % (_ids,data_file,long2short) ).close()
  f_a = os.popen('python %s -f %d -t %d %s --filtercol=8 --filterval=!Failed %s ' \
                % (_ids,col,to_col,assay_def,locus_map) ).close()
  if f or f_s or f_a:
    raise ValueError,'Failed in creating long2short map, or locus_map, or sid map file(s)!'

  # locus includes; long2short map
  _transform = os.path.abspath(path_dir['script_dir']) + '/transform2.py'
  trip_v1    = work_dir+ '/' + proj_name + '_v1.trip'
  mergestat_sample_0 = work_dir+ '/mergestat_SB.txt'
  f = os.popen('python -u %s -f sdat %s -F trip -o %s -u %s -m %s --samplemerge=%s '
            % (_transform,data_file,trip_v1,locus_map,long2short,mergestat_sample_0) ).close()
  if f:
    raise ValueError,'Failed in creating trip_v1(long SB to short SB)!'

  # sample includes; sample and locus remap
  sdat_v2          = work_dir + '/' + proj_name + '_v2.sdat'
  mergestat_sample = res_dir  + '/mergestat_sid.txt'
  mergestat_locus  = res_dir  + '/mergestat_locus.txt'
  f = os.popen('python -u %s -f trip %s -F sdat -o %s -m %s -n %s -r %s --samplemerge=%s --locusmerge=%s'
            % (_transform,trip_v1,sdat_v2,sid_map,sample_include,locus_map,mergestat_sample,mergestat_locus) ).close()
  if f:
    raise ValueError,'Failed in creating sdat_v2(rename samples and loci)!'

  # final sdat
  sdat_final = res_dir + '/' + proj_name + '_validated.sdat'
  f = os.popen('python -u %s -f sdat %s -F sdat -o %s -c' % (_transform,sdat_v2,sdat_final) ).close()
  if f:
    raise ValueError,'Failed in creating final sdat!'

  # final reports
  if 'yes' in path_dir['tabular_report'].lower():
    compl_dat   = work_dir+ '/' + proj_name +'_completion.dat'
    #compl_out   = work_dir+ '/' + proj_name +'_completion.out'
    _completion = os.path.abspath(path_dir['script_dir']) + '/completion.py'
    #f = os.popen('python -u %s -f sdat %s  --tabularoutput %s > %s' \
    #       %(_completion,sdat_v2,compl_dat,compl_out) ).close()
    compl_out   = open(work_dir+ '/' + proj_name +'_completion.out','w')
    p1 = Popen(['python',_completion, '-f', 'sdat', sdat_v2, '--tabularoutput',compl_dat],stdout=PIPE)
    compl_out.write(p1.communicate()[0])
    compl_out.close()

    if f:
      raise ValueError,'Failed in completion (for generating tabular report)!'

    sample_prn  = res_dir + '/' + proj_name + '_sample.txt'
    locus_prn   = res_dir + '/' + proj_name + '_locus.txt'
    snpid       = path_dir['report_assayid']
    _tabreport  = os.path.abspath(path_dir['script_dir']) + '/tabreports2.py'
    f = os.popen('python -u %s -c %s -a %s -m %s -l %s -s %s -i %s -r %s' \
              % (_tabreport,compl_dat,assay_def,n_sample_def,locus_prn,sample_prn,snpid,strand) ).close()
    if f:
      raise ValueError,'Failed in creating sample or locus prn files (in tabular report)!'

    if 'yes' in path_dir['summary_report'].lower():
      summary     = res_dir  + '/' + proj_name + '_summary'
      proj_def    = final_dir+ '/' + proj_name + '_project.def'
      qc          = final_dir+ '/QC_Report.txt'
      _sumreports = os.path.abspath(path_dir['script_dir']) + '/summaryreports.py'
      f = os.popen('python2.4 -u %s -c %s -a %s -m %s -p %s -s %s -t %s' \
                % (_sumreports,compl_dat,assay_def,sample_prn,proj_def,summary,qc) ).close()
      if f:
        raise ValueError,'Failed in creating summary pdf report!'

  if 'yes' in path_dir['analysis_ready_data'].lower():
    #pheno_file = work_dir + '/' + proj_name + '_Pheno.txt'
    #f = os.popen('cp %s %s ' % (path_dir['pheno_file'], pheno_file) ).close()
    #if f:
    #  raise ValueError,'Failed in copy pheno_file!'
    pheno_file = path_dir['pheno_file']
    _pids      = os.path.abspath(path_dir['script_dir2']) + '/pids.py'
    #pid_1
    pheno_2    = work_dir + '/' + proj_name + '_Pheno.txt'
    f = os.popen('python %s  %s  %s --cols=1,2,7' \
              % (_ids,pheno_file,pheno_2) ).close()
    if f:
      raise ValueError,'Failed in subsetting pheno file (for analysis ready data)'
    #pid
    pid_file   = work_dir + '/' + 'pid.txt'
    sid_include= work_dir + '/' + 'sid_include'
    f = os.popen('python %s  %s  %s --sidout=%s' \
              % (_pids,pheno_2,pid_file,sid_include) ).close()
    if f:
      raise ValueError,'Failed in creating pheno file withou QC (for analysis ready data)'
    #print 'after qc','check pip error'
    #rename2
    trip_v3   = work_dir+ '/' + proj_name + '_v3.trip'
    f = os.popen('python -u %s -f trip %s -F trip -o %s -m %s -n %s -r %s'
              % (_transform,trip_v1,trip_v3,sid_map,sample_include,locus_map) ).close()
    if f:
      raise ValueError,'Failed in creating trip_v3(for analysis ready data)!'
    #rename3
    sdat_pid         = res_dir + '/' + proj_name + '_pid.sdat'
    mergestat_sample = res_dir + '/' + proj_name + '_mergestat_pid.txt'
    f = os.popen('python -u %s -f trip %s -F sdat -o %s -m %s -n %s --samplemerge=%s -c'
              % (_transform,trip_v3,sdat_pid,sid_include,sid_include,mergestat_sample) ).close()
    if f:
      raise ValueError,'Failed in creating analysis ready sdata(pid based sdat)!'
    #subset
    proj_info = res_dir + '/' + proj_name + '.info'
    f = os.popen('python %s  %s  %s --filtercol=1 --rowinclude=%s -T' \
              % (_ids,pid_file,proj_info,sdat_pid) ).close()
    if f:
      raise ValueError,'Failed in createing .info file'

  return None


def process(path_dir):
  pf = path_dir['platform']
  if pf in ('taqman','snplex'):
    prep_taq_data(path_dir)
    col = 2
    strand = 'gene'
  elif pf == 'gg':
    prep_gg_data(path_dir)
    col = 3
    strand = 'contig'
  else:
    raise ValueError,'Stopped: Platformat has to be taqman, snplex, or gg'

  create_report(path_dir,col,strand)


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
