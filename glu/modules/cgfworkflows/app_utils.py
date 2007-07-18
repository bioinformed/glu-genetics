# -*- coding: utf-8 -*-
'''
File:          app_utils.py

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      Utility function for web-applications

Requires:      Python 2.5 or higher

Revision:
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys
import csv

from   itertools         import islice,imap,chain,groupby
from   operator          import itemgetter

from   glu.lib.utils     import autofile,peekfirst
from   glu.lib.genoarray import snp_acgt,get_genorepr
from   glu.lib.genodata  import  *
from   glu.lib.genomerge import VoteMerger

def mapper2(astr):
  if '~' in astr:
    tt = astr.split('~')
    jj = '~'.join([tt[0],tt[2]])
    return jj
  if '_' in astr:
    tt = astr.split('_')
    jj = '_'.join([tt[0],tt[2]])
    return jj
  return astr


def mapper(astr):
  if '~' in astr:
    jj = astr.split('~')[0]
    return jj
  if '_' in astr:
    tt = astr.split('_')
    jj = tt[0]
    if len(tt)>3:
      print >> sys.stderr,'Warnings: Multiple underscores in the sample names!'
      jj = '_'.join(tt[:-2])
    return jj
  return astr


def clean_dir(res_dir,files):
  res_dir = os.path.abspath(res_dir)
  for f in files:
    if os.path.exists(res_dir+'/'+f):
      os.remove(res_dir+'/'+f)
    if os.path.exists(res_dir+'/workspace/'+f):
      os.remove(res_dir+'/workspace/'+f)


def get_dir(proj_name,data_type='QC'):
  if proj_name.startswith('GR'):
    return '/'.join([proj_name,'analysis','data'])
  elif proj_name.startswith('GP'):
    p_num = proj_name.split('-')[0][2:]
    return '/'.join(['MR-'+p_num,proj_name,'Analysis',data_type,data_type+'Data'])
  else:
    raise ValueError,'Project name should start with GP or GR; the current name is %s!' % proj_name


def trip2sdat(trip_file,sdat_file,sample_map=None,assay_map=None,outgenorepr='snp_acgt'):
  trip = autofile(trip_file)
  header,trip = peekfirst(trip)
  if 'SAMPLE' in header:
    trip.next()
  ingenorepr  = get_genorepr('generic')
  outgenorepr = get_genorepr(outgenorepr)
  merger = VoteMerger(threshold=1.0,missingrepr=0)
  genos  = load_genostream(trip,'trip',genorepr=ingenorepr).transformed(rename_samples=sample_map,rename_loci=assay_map)
  save_genostream(sdat_file,genos.recoded(outgenorepr),'sdat',merger)


def trip2ldat(trip_file,ldat_file,sample_map=None,assay_map=None,outgenorepr='snp_acgt'):
  trip = autofile(trip_file)
  header,trip = peekfirst(trip)
  if 'SAMPLE' in header:
    trip.next()
  ingenorepr  = get_genorepr('generic')
  outgenorepr = get_genorepr(outgenorepr)
  merger = VoteMerger(threshold=1.0,missingrepr=0)
  genos  = load_genostream(trip,'trip',genorepr=ingenorepr).transformed(rename_samples=sample_map,rename_loci=assay_map)
  save_genostream(ldat_file,genos.recoded(outgenorepr),'ldat',merger)


def sdat2ldat(sdat_f,ldat_f):
  infiles = [sdat_f]
  ingenorepr  = get_genorepr('snp_acgt')
  genos   = transform_files(infiles,'sdat',ingenorepr,'ldat',ingenorepr)
  merger = VoteMerger(threshold=1.0,missingrepr=0)
  save_genostream(ldat_f,genos,'ldat',merger)


def output_map(infile,outfile,cols,nskiprows=0):
  rows = csv.reader(autofile(infile),dialect='excel-tab')
  rows = islice(rows,nskiprows,None)

  out = csv.writer(autofile(outfile, 'w'),dialect='excel-tab')
  for row in rows:
    out.writerow([row[i-1] for i in cols])


def get_OPA_id(assay_file):
  rows = csv.reader(autofile(assay_file),dialect='excel-tab')
  rows.next()
  r = rows.next()
  return r[0]


def create_map(data_dir,proj_name,work_dir,sdat,ids,assay_from,assay_to):
  sampledef= data_dir + '/'+proj_name+'_Sample.txt'
  assaydef = data_dir + '/'+proj_name+'_Assay.txt'

  output_map(assaydef,work_dir + '/locus_map',[assay_from,assay_to])
  output_map(assaydef,work_dir + '/locus_include',[assay_to,13])

  output_map(sampledef,work_dir + '/sample_include',[2])
  sample_map = work_dir + '/sample_map'
  f = os.popen('python %s -s 1 -S 0 -v %s %s --sidfile=%s ' % (ids,sdat,sample_map,sampledef) ).close()
  if f:
    raise ValueError,'Failed in creating sample_map!'


def prep_gg_data(path_dir):
  data_dir = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir  = os.path.abspath(path_dir['result_dir'])
  work_dir = res_dir+'/'+'workspace'

  check_data_dir(data_dir,proj_name,res_dir)

  sampledef = data_dir + '/'+proj_name+'_Sample.txt'
  assaydef = data_dir + '/'+proj_name+'_Assay.txt'

  s_data      = work_dir + '/'+proj_name+'.sdat'
  lbd2sdat    = os.path.abspath(path_dir['script_dir']) + '/lbd2sdat.py'

  opa_id = get_OPA_id(assaydef)
  ab_map = os.path.abspath(path_dir['abmap_dir']) + '/'+ opa_id + '_nt.ab'
  assert os.path.exists(ab_map),'Require ab_map file: %s ' % ab_map

  lbd_file = data_dir +'/'+proj_name+'.lbd'
  assert os.path.exists(lbd_file), 'Require lbd file: %s ' % lbd_file
  f = os.popen('python -u %s %s -o %s -m %s' % (lbd2sdat,lbd_file,s_data,ab_map)).close()
  if f:
    raise ValueError,'Failed in transforming lbd to sdat ' % lbd_file

  # subset
  if 'include_loci' in path_dir and os.path.exists(path_dir['include_loci']):
    _transform2  = os.path.abspath(path_dir['script_dir']) + '/transform2.py'
    os.popen('mv %s %s' % (s_data,s_data+'.orig')).close()
    f = os.popen('python -u %s -f sdat %s -o %s -u %s' \
              % (_transform2,s_data+'.orig',s_data,path_dir['include_loci'])).close()


def prep_taq_data(path_dir,outgenorepr='snp_acgt'):
  data_dir = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir  = os.path.abspath(path_dir['result_dir'])
  work_dir = res_dir+'/'+'workspace'

  check_data_dir(data_dir,proj_name,res_dir)

  s_data      = work_dir + '/'+proj_name+'.sdat'

  ids = os.path.abspath(path_dir['script_dir2']) + '/ids.py'
  trip_file  = data_dir +'/'+proj_name+'.trip'
  assert os.path.exists(trip_file), 'Require trip file: %s ' % trip_file
  sample_map = work_dir + '/remap_sample'
  f = os.popen('python %s -s 0 -S 0 -v -w %s %s ' % (ids,trip_file,sample_map) ).close()
  if f:
    raise ValueError,'Failed in creating sample_map (remove plate# from sample ID)!'

  trip2sdat(trip_file,s_data,outgenorepr=outgenorepr,sample_map=sample_map)

  # subset
  if 'include_loci' in path_dir and os.path.exists(path_dir['include_loci']):
    _transform2  = os.path.abspath(path_dir['script_dir']) + '/transform2.py'
    os.popen('mv %s %s' % (s_data,s_data+'.orig')).close()
    f = os.popen('python -u %s -f sdat -g generic %s -o %s -u %s' \
              % (_transform2,s_data+'.orig',s_data,path_dir['include_loci'])).close()


def prep_ident_data(path_dir):
  data_dir = os.path.abspath(path_dir['data_dir'])
  proj_name= path_dir['proj_name']
  res_dir  = os.path.abspath(path_dir['result_dir'])
  work_dir = res_dir+'/'+'workspace'

  check_data_dir(data_dir,proj_name,res_dir)
  s_data    = work_dir + '/'+proj_name+'_Identifiler.sdat'
  trip_file = data_dir + '/'+proj_name+'_Identifiler.trip'
  assert os.path.exists(trip_file), 'Require trip file: %s ' % trip_file

  trip2sdat(trip_file,s_data,outgenorepr='generic')


def prep_dir(raw_dir):
  new_dir = {}
  for k,v in raw_dir.iteritems():
    new_dir[k] = v.strip()

  if 'result_dir' in new_dir:
    new_dir['result_dir'] = new_dir['result_dir'] + '/' + new_dir['proj_name']

  return new_dir


def check_data_dir(data_dir,proj_name,res_dir):
  sampledef = data_dir + '/' + proj_name + '_Sample.txt'
  assaydef  = data_dir + '/' + proj_name + '_Assay.txt'

  assert os.path.exists(data_dir),  'Directory: %s not found!' % data_dir

  assert os.path.exists(sampledef), 'File: %s not found!' % sampledef
  assert os.path.exists(assaydef),  'File: %s not found!' % assaydef

  res_dir  = os.path.abspath(res_dir)
  up_resdir = os.path.dirname(res_dir)
  assert os.path.exists(up_resdir),   'Directory: %s not found!' % up_resdir
  work_dir  = res_dir + '/workspace'
  report_dir= res_dir + '/reports'

  for d in (work_dir,report_dir):
    if not os.path.exists(d):
      os.makedirs(d)

  return None


def guess_genorepr(sdat_f):
  f = csv.reader(autofile(sdat_f),dialect='excel-tab')
  f.next()
  for row in f:
    genos = [r for r in row if r]
    if len(genos)>1:
      if '/' in genos[1]:
        return 'generic'
      elif ',' in genos[1]:
        return 'snp_marker'

  return 'snp_acgt'


def guess_platform(proj):
  a,b = proj.split('-')
  pf  = 'unknown'
  if b.startswith('GG'):
    pf = 'gg'
  elif b.startswith('UN'):
    pf = 'taqman'
  elif b.startswith('SN'):
    pf = 'snplex'

  return pf
