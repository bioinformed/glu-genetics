# -*- coding: utf-8 -*-
'''
File:          lbd2sdat.py

Authors:       Brian Staats (staatsb@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-15

Abstract:      Parses an Illumina Locus by DNA report file and outputs a
               genotype matrix or triple stream.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from glu.modules.convert.from_lbd import main

if __name__ == '__main__':
  main()
