# -*- coding: utf-8 -*-
'''
File:          __init__.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU genotype data input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

# FIXME: Format support should ultimately be pluggable with a registration protocol
from   glu.lib.genolib.formats.text       import *
from   glu.lib.genolib.formats.prettybase import *
from   glu.lib.genolib.formats.hapmap     import *
from   glu.lib.genolib.formats.binary     import *
