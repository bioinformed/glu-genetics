# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      dupcheck for pregenotype identifiler data

Requires:      Python 2.4 or higher

Revision:
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys,os
import csv

from core_report   import *


def option_parser():
  import optparse

  usage = 'usage: %prog [options] project_name'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-o', '--outdir', dest='outdir', metavar='string',
                    help='specifying the directory for the project output')
  parser.add_option('-i', '--inputdir', dest='inputdir', metavar='string',
                    help='specifying the directory of input data, e.g. trip file')
  parser.add_option('-f', '--platform', dest='platform', metavar='string',
                    help='specifying the platform of the project, values=gg,snplex,taqman')

  parser.add_option('-d', '--finaldir', dest='finaldir', metavar='string',
                    help='specifying the final directory')
  parser.add_option('-a', '--assayid', dest='assayid', metavar='string',default='dbsnp',
                    help='specifying the reporting assayid, value=dbsnp(default),external,alias')
  parser.add_option('-t', '--tabreport', dest='tabreport', metavar='string', default='yes',
                    help='specifying whether to create tabular report, values = yes(default),no')
  parser.add_option('-s', '--sumreport', dest='sumreport', metavar='string', default='yes',
                    help='specifying whether to create pdf summary report, values = yes(default),no')
  parser.add_option('-r', '--analysisdata', dest='analysisdata', metavar='string', default='no',
                    help='specifying whether to create analysis ready data, values = no(default),yes')

  parser.add_option('-p', '--phenofile', dest='phenofile', metavar='string',default='none',
                    help='specifying location of PI phenotype file')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    print >> sys.stderr, 'Program only accepts 1 argument (a project name)'
    return

  if not options.outdir:
    print >> sys.stderr,'Stopped: need to specify an output directory'
    return

  path_dir = load_map('config_app')
  assert path_dir

  path_dir['proj_name'] = args[0]
  path_dir['data_dir']  = os.path.abspath(options.inputdir)
  path_dir['result_dir']= options.outdir
  path_dir['platform']  = options.platform

  path_dir['final_dir'] = options.finaldir
  path_dir['report_assayid']      = options.assayid
  path_dir['tabular_report']      = options.tabreport
  path_dir['summary_report']      = options.sumreport
  path_dir['analysis_ready_data'] = options.analysisdata
  path_dir['pheno_file']          = options.phenofile

  path_dir = prep_dir(path_dir)
  process(path_dir)


if __name__ == '__main__':
  main()
