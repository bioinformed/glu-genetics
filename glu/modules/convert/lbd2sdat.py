# -*- coding: utf-8 -*-
'''
File:          lbd2sdat.py

Authors:       Brian Staats (staatsb@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-15

Abstract:      Parses an Illumina Locus by DNA report file and outputs a
               genotype matrix or triple stream.

Requires:      Python 2.5, glu

Revision:      $Id: lbd2sdat.py 432 2007-12-03 01:17:07Z jacobske $
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

from glu.modules.convert.from_lbd import main

if __name__ == '__main__':
  main()
