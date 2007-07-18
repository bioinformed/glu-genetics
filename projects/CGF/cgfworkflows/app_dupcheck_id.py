# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      specifiy a project name and a data source, then output sdat files

Requires:      Python 2.4 or higher

Revision:
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys,os
import csv

from core_dupcheck_id   import *


def option_parser():
  import optparse

  usage = 'usage: %prog [options] project_name'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=85,
                    help='Threshold for the percentage of identity shared between two individuals (default=85)')
  parser.add_option('-m', '--mincount', dest='mincount', metavar='N', type='int', default=5,
                    help='Minimum number of concordant genotypes required for duplicate checking (default=5)')
  parser.add_option('-o',  '--outdir', dest='outdir', metavar='string',
                    help='specifying the directory for the project output')

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

  proj_name = args[0]
  path_dir['proj_name'] = proj_name
  path_dir['data_dir']  = os.path.abspath(path_dir['cgfdata_dir']) + '/' + get_dir(proj_name)
  path_dir['result_dir']= options.outdir
  path_dir['Ident_threshold'] = str(options.threshold)
  path_dir['Ident_mincount']  = str(options.mincount)

  path_dir['Ident_with_expdup'] = 'yes'
  path_dir['Ident_dupcheck_lite'] = 'no'

  path_dir['pre_geno_qc'] = 'no'

  path_dir = prep_dir(path_dir)
  process(path_dir)


if __name__ == '__main__':
  main()
