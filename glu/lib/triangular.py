# -*- coding: utf-8 -*-
'''
File:          triangular.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

# Routines for computing triangular numbers and their inverses
from math import sqrt


def triangle(n):
  return n*(n+1)//2


def triangle_inv(x):
  return int((sqrt(8*x+1)-1)//2)
