# -*- coding: utf-8 -*-

__abstract__  = 'Compute triangular numbers and their inverse'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

# Routines for computing triangular numbers and their inverses
from math import sqrt


def triangle(n):
  return n*(n+1)//2


def triangle_inv(x):
  return int((sqrt(8*x+1)-1)//2)
