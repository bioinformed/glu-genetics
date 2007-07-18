# -*- coding: utf-8 -*-
'''
File:          app_completion.py

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      specifiy a project name and a platform, and the completion results will be output to a directory

Requires:      Python 2.5 or higher

Revision:
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import os
import sys

from   glu.lib.genodata import load_map
from   app_utils        import get_dir,guess_platform
from   core_completion  import *

def option_parser():
  import optparse

  usage = 'usage: %prog project_name [options]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o',  '--outdir', dest='outdir', metavar='string',
                    help='specifying the directory for the project output')

  parser.add_option('-p',  '--platform', dest='platform', metavar='string',
                    help='specifying a platform, required (value=gg,taqman or snplex)')

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

  path_dir['proj_name'] = args[0]
  path_dir['data_dir']  = os.path.abspath(path_dir['cgfdata_dir']) + '/' + get_dir(args[0])
  path_dir['result_dir']= options.outdir
  path_dir['platform']  = options.platform

  path_dir = prep_dir(path_dir)
  process(path_dir)


if __name__ == '__main__':
  main()
