# -*- coding: utf-8 -*-
'''
File:          core_concordance.py

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      concordance check using a config file as input

Requires:      Python 2.5 or higher

Revision:
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys,os
import csv

from   glu.lib.utils       import autofile,hyphen
from   glu.lib.genoarray   import snp_acgt,get_genorepr
from   glu.lib.genodata    import load_map
from   glu.lib.genomerge   import VoteMerger
from   app_utils           import *
from   scripts.concordance import *


def concordance_check(path_dir):
  data_dir  = os.path.abspath(path_dir['data_dir'])
  proj_name =path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  comp_dir  = os.path.abspath(path_dir['comp_dir'])
  comp_format = path_dir['comp_format']

  # get chromosome info
  sample_include = load_map(work_dir + '/sample_include')
  locus_include  = load_map(work_dir + '/locus_include')

  comp_files = []
  if 'hapmap' in comp_format:
    hapsamples= load_map(os.path.abspath(path_dir['hapsamples']))
    sample_include = load_map(work_dir + '/sample_include')
    locus_include  = load_map(work_dir + '/locus_include')
    rows = csv.reader(autofile(work_dir+ '/locus_include'),dialect='excel-tab')
    rows = islice(rows,1,None)
    chrs = set(row[1] for row in rows)
    pops = set(hapsamples[s] for s in sample_include if s in hapsamples)
    print >> sys.stderr, '\nExtracting hapmap data on %d chromosomes and %d population(s) ...' % ( len(chrs),len(pops) )
    file_s = comp_dir+'/genotypes_chr'
    for a in pops:
      for b in chrs:
        if a=='CEU':
          comp_files.append(file_s+b+'_'+a+'_r21a_nr_fwd.txt')
        else:
          comp_files.append(file_s+b+'_'+a+'_r21a_nr_fwd.txt.gz')

  elif 'ldat' in comp_format:
    import glob
    comp_files = glob.glob(comp_dir+'/*')

  else:
    raise ValueError, 'comp_format has to be either ldat or hapmap!'

  assert comp_files,'Comparison files (e.g. hapmap) not found!'
  # extract comp data (e.g. hapmap or other ldat)
  comp_dat = work_dir + '/comp.ldat'
  transform = GenoTransform.from_kwargs(include_samples=sample_include, include_loci=locus_include)
  merge   = VoteMerger()
  genos   = transform_files(comp_files,comp_format,snp_acgt,'ldat',snp_acgt,transform=transform,mergefunc=merge)
  save_genostream(comp_dat,genos,'ldat')

  # load ref and comp geno data
  locuseq = eqlocus = None
  locuseq = load_map(work_dir + '/locus_map')
  eqlocus = invert_dict(locuseq)

  sampleeq = eqsample = None
  sampleeq = load_map(work_dir + '/sample_map')
  eqsample = invert_dict(sampleeq)

  refgenos,samples,loci = load_reference_genotypes(work_dir+'/'+proj_name+'.ldat','ldat',locuseq,sampleeq,None)
  compgenos = load_comparison_genotypes(comp_dat,'ldat', eqlocus, eqsample,None,None)

  # generate allele remap
  sampleconcord = SampleConcordStat()
  locusconcord  = LocusConcordStat()

  concordance(refgenos,samples,loci,compgenos,eqsample,eqlocus,sampleconcord,locusconcord)

  amap = compute_allele_maps(locusconcord)
  output_allele_maps(amap,work_dir+'/allele_map')

  # concordance check
  compgenos = load_comparison_genotypes(comp_dat,'ldat', eqlocus, eqsample,None,None)
  allelemaps   = load_remap_file(work_dir+'/allele_map')
  genomappings = build_geno_remap(allelemaps)
  compgenos    = remap_genotypes(compgenos, genomappings)

  sampleconcord = SampleConcordStat()
  locusconcord  = LocusConcordStat()

  concordance(refgenos,samples,loci,compgenos,eqsample,eqlocus,sampleconcord,locusconcord)

  output_sample_concordstat(res_dir+'/'+ proj_name + '_concord_by_sample.out',sampleconcord)
  output_locus_concordstat(res_dir+'/'+proj_name + '_concord_by_locus.out',locusconcord,allelemaps)


def process(path_dir):
  data_dir= path_dir['data_dir']
  proj_name =path_dir['proj_name']
  res_dir   = os.path.abspath(path_dir['result_dir'])
  work_dir  = res_dir + '/workspace'

  # clean old files
  FILES = [proj_name+'_concord_by_sample.out',proj_name+'_concord_by_locus.out']
  FILES = FILES + ['locus_map','locus_include','sample_map','sample_include']
  FILES = FILES + [proj_name+'.ldat',proj_name+'.sdat','comp.ldat','allele_map']
  clean_dir(res_dir,FILES)

  ids= os.path.abspath(path_dir['script_dir2']) + '/ids.py'
  s_data = work_dir + '/'+proj_name+'.sdat'
  pf = path_dir['platform']
  if pf == 'gg':
    prep_gg_data(path_dir)
    create_map(data_dir,proj_name,work_dir,s_data,ids,3,4)
  elif pf in ('taqman','snplex'):
    prep_taq_data(path_dir)
    create_map(data_dir,proj_name,work_dir,s_data,ids,2,4)
  else:
    raise ValueError,'Platform should be "gg","taqman",or "snplex" !'

  l_data  = work_dir + '/'+proj_name+'.ldat'

  sdat2ldat(s_data,l_data)
  concordance_check(path_dir)


def option_parser():
  import optparse

  usage = 'usage: %prog configFile'
  parser = optparse.OptionParser(usage=usage)

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    print >> sys.stderr, 'Stopped: require a config file!'
    return

  raw_dir = load_map(args[0])
  path_dir= prep_dir(raw_dir)

  process(path_dir)


if __name__ == '__main__':
  main()
