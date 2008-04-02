# -*- coding: utf-8 -*-
'''
File:          __init__.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU genotype data input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


# FIXME: Format support should ultimately be pluggable with a registration protocol
from   glu.lib.genolib.formats.text       import *
from   glu.lib.genolib.formats.binary     import *
from   glu.lib.genolib.formats.prettybase import *
from   glu.lib.genolib.formats.hapmap     import *
from   glu.lib.genolib.formats.plink      import *
from   glu.lib.genolib.formats.merlin     import *
from   glu.lib.genolib.formats.structure  import *
from   glu.lib.genolib.formats.phase      import *
from   glu.lib.genolib.formats.wtccc      import *
from   glu.lib.genolib.formats.eigensoft  import *
