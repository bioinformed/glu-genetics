# -*- coding: utf-8 -*-
'''
File:          app_dupcheck.py

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      specifiy a project name and a data source, then output sdat files

Requires:      Python 2.5 or higher

Revision:
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys

from   glu.lib.genodata import load_map
from   app_utils        import guess_platform,get_dir
from   core_dupcheck    import *

def option_parser():
  import optparse

  usage = 'usage: %prog [options] project_name'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-p',  '--platform', dest='platform', metavar='string',
                    help='specifying a platform, required (value=gg,taqman or snplex)')
  parser.add_option('-T', '--threshold', dest='threshold', metavar='N%', type='int', default=90,
                    help='Threshold for the percentage of identity shared between two individuals (default=90)')
  parser.add_option('-m', '--mincount', dest='mincount', metavar='N', type='int', default=20,
                    help='Minimum number of concordant genotypes required to be duplicates (default=20)')
  parser.add_option('-e',  '--expdups', dest='expdups', metavar='string',default=None,
                    help='specifying type of expected duplicates,(value=iddups, pidups,default=None)')
  parser.add_option('-l',  '--lite', dest='lite', action='store_true', default=False,
                    help='Run dupcheck lite (i.e. ignore unexpected duplicates)')

  parser.add_option('-o',  '--outdir', dest='outdir', metavar='string',
                    help='specifying the directory for the project output')
  parser.add_option('-s',  '--preprocess', dest='preprocess', metavar='string', default='yes',
                    help='specifying whether preprocess step is needed, value=yes(default),no')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  if not options.outdir:
    print >> sys.stderr,'Stopped: need to specify an output directory'
    return

  if not options.platform:
    print >> sys.stderr,'Stopped: need to specify a platform'
    return

  path_dir = load_map('./config_app')
  assert path_dir

  proj_name = args[0]
  path_dir['proj_name'] = proj_name
  path_dir['data_dir']  = os.path.abspath(path_dir['cgfdata_dir']) + '/' + get_dir(proj_name)
  path_dir['result_dir']= options.outdir
  path_dir['Geno_threshold'] = str(options.threshold)
  path_dir['Geno_mincount']  = str(options.mincount)
  path_dir['platform']       = options.platform
  path_dir['Geno_dupcheck_lite'] = 'yes' if  options.lite else 'no'

  path_dir['preprocess'] = options.preprocess

  path_dir['Geno_with_expdup']   = 'yes'
  if not options.expdups:
    path_dir['Geno_with_expdup'] = 'no'
  elif options.expdups == 'iddups':
    path_dir['Geno_expdup_type'] = 'ident_dups'
  elif options.expdups == 'pidups':
    path_dir['Geno_expdup_type'] = 'pi_dups'
  else:
    raise ValueError, 'Invalid expected dup types: %s ' % options.expdups

  path_dir = prep_dir(path_dir)

  process(path_dir)


if __name__ == '__main__':
  main()
