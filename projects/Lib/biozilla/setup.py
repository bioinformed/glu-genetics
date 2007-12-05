#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# To use:
#       python setup.py install
#

import distutils, os
from   distutils.core import setup, Extension

import numpy

setup (name = 'biozilla',
       version = '0.01',
       maintainer = 'Kevin Jacobs',
       maintainer_email = 'jacobske@mail.nih.gov',
       description = '',
       packages  = ['db_row','biozilla'],
       ext_modules = [ Extension('db_row/db_rowc',       ['src/fields.c','src/row.c'],
                                                         include_dirs = ['./src']),
                       Extension('biozilla/bitarrayc',   ['biozilla/bitarrayc.c']),
                       Extension('biozilla/_genoarray',  ['biozilla/_genoarray.c','biozilla/bitarrayc.c'],
                                                         include_dirs = [numpy.get_include()]),
                     ])
