# -*- coding: utf-8 -*-

__abstract__  = 'Run GLU lib doctest suite'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import pkg_resources
pkg_resources.require('glu')

def test(report=True):
  import os
  import doctest

  import glu
  from   glu.lib import (utils,fileutils,association,glm,imerge,sections,sequence,stats,
                         utils,hwp,xtab,union_find,regionparser)
  import utils as local_utils

  # Standardize on the .py source, since one can load the .py and the other the .pyc
  packagefile =  utils.__file__.replace('.pyc','.py')
  localfile   =  local_utils.__file__.replace('.pyc','.py')

  try:
    if not os.path.samefile(packagefile,localfile):
      raise OSError
  except OSError:
    raise ImportError('Your PYTHONPATH/pkg_resources are not set correctly to find this GLU tree')

  for module in association,fileutils,glm,imerge,sections,sequence,stats,utils,hwp,xtab,union_find,regionparser:
    doctest.testmod(module,report=False)

  from glu.lib.genolib.test import test as test_genolib

  test_genolib(False)

  if report:
    doctest.master.summarize()


if __name__ == '__main__':
  test()
