#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# To use:
#       python setup.py install
#

import os
import sys


from ez_setup import use_setuptools
use_setuptools()


# minimum versions required
min_python_version = (2,5)
min_numpy_version  = '1.0.2'
min_scipy_version  = '0.5.2'
min_tables_version = '2.0'
min_ply_version    = '2.5'
min_nose_version   = '0.10.4'

if sys.version_info < min_python_version:
  sys.stderr.write('Python 2.5 or newer to required to install and run GLU!\n')
  sys.exit(1)

import distutils
from   setuptools           import setup, find_packages, Extension
from   numpy.distutils.core import setup as numpy_setup
    
import numpy

def get_version():
  vfile = os.path.join(os.path.dirname(__file__), 'glu', 'VERSION')
  version = file(vfile).readline().strip()
  return version


classifiers = '''\
Development Status :: 5 - Production/Stable
Home-page: http://code.google.com/p/glu-genetics/
Intended Audience :: Developers
Intended Audience :: Information Technology
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Programming Language :: Python :: 2.5
Programming Language :: C
Topic :: Scientific/Engineering
Topic :: Software Development
Topic :: Database
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: Unix
Operating System :: POSIX
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
'''

#entry_points    = { 'console_scripts':['glu = glu.lib.glu_launcher:main'] },


def glmnet_config(parent_package='', top_path=None):
  from numpy.distutils.misc_util import Configuration
  config = Configuration('_glmnet', parent_package, top_path)
  return config.add_extension('glu.lib.glm._glmnet', sources=['glu/lib/glm/glmnet.pyf','glu/lib/glm/GLMnet.f'])
                

def main():
  numpy_setup(name             = 'glu',
        version          = get_version(),
        author           = 'Kevin Jacobs',
        author_email     = 'jacobske@mail.nih.gov',
        maintainer       = 'Kevin Jacobs',
        maintainer_email = 'jacobske@mail.nih.gov',
        platforms        = ['any'],
        description      = 'Genotype Library and Utilities (GLU)',
        long_description = ('Genotype Library and Utilities (GLU): Tools for the management of large '
                            'amounts of SNP genotype data and programs to check its quality and to '
                            'test for association between SNP markers with continuous or discrete '
                            'trait phenotypes.'),
        classifiers      = filter(None, classifiers.split('\n')),
        install_requires = ['numpy>=%s'  % min_numpy_version,
                            'scipy>=%s'  % min_scipy_version,
                            'tables>=%s' % min_tables_version,
                            'ply>=%s'    % min_ply_version  ],
        setup_requires   = ['numpy>=%s'  % min_numpy_version],
        tests_require    = ['nose>=%s'   % min_nose_version ],
        packages         = find_packages(),
        include_package_data = True,
        scripts          = ['bin/glu'],
        zip_safe         = False,
        test_suite       = 'nose.collector',
        ext_modules = [ Extension('glu.lib.genolib.bitarrayc',      sources = ['glu/lib/genolib/bitarrayc.c']),
                        Extension('glu.lib.genolib._genoarray',     sources = ['glu/lib/genolib/_genoarray.c',
                                                                               'glu/lib/genolib/bitarrayc.c',
                                                                               'glu/lib/genolib/_ibs.c',
                                                                               'glu/lib/genolib/_ld.c'],
                                                                    include_dirs = [numpy.get_include()]),
                        Extension('glu.modules.ld.pqueue',    sources = ['glu/modules/ld/pqueue.c']),
                        glmnet_config('glu.lib'),
                      ])


if __name__ == '__main__':
  main()
