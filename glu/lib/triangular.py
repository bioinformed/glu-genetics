# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Compute triangular numbers, their inverses, and rectangular coordinate transformations'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


# Routines for computing triangular numbers and their inverses
from math import sqrt


def triangle(n):
  '''triangle(n) -> the n'th triangular number'''
  return n*(n+1)//2


def triangle_inv(x):
  '''triangle_inv(x) -> the smallest triangular number n such that triangle(n)<=x'''
  return int((sqrt(8*x+1)-1)//2)


def triangle_coords(n):
  '''triangle_coords(n) -> Rectangular coordinates x,y such that n=triangle(x)+y and 0<=y<=x'''
  x = triangle_inv(n)
  y = n-triangle(x)
  return x,y


def triangle_index(x,y):
  '''triangle_index(x,y) -> Returns n=triangle(x)+y, for 0<=y<=x, the inverse of triangle_coords(n)'''
  return triangle(x) + y


def strict_lower_triangle(n):
  '''strict_lower_triangle(n) -> the n'th strictly lower triangular number'''
  return n*(n-1)//2


def strict_lower_triangle_inv(x):
  '''strict_lower_triangle_inv(x) -> the smallest strictly lower triangular number n such that strict_lower_triangle(n)<=x'''
  return int((sqrt(8*x+1)+1)//2)


def strict_lower_triangle_coords(n):
  '''strict_lower_triangle_coords(n) -> Rectangular coordinates x,y such that n=strict_lower_triangle(x)+y and 0<=y<x'''
  x = strict_lower_triangle_inv(n)
  y = n-strict_lower_triangle(x)
  return x,y


def strict_lower_triangle_index(x,y):
  '''strict_lower_triangle_index(x,y) -> Returns n=strict_lower_triangle(x)+y, for 0<=y<x, the inverse of strict_lower_triangle_coords(n)'''
  return strict_lower_triangle(x) + y


def test():
  for i in range(1000):
    x,y = triangle_coords(i)
    assert 0<=y<=x and triangle(x)+y == i

    x,y = strict_lower_triangle_coords(i)
    print x,y
    assert 0<=y<x and strict_lower_triangle(x)+y == i

    assert i==strict_lower_triangle_index(*strict_lower_triangle_coords(i))
    assert i==triangle_index(*triangle_coords(i))


if __name__ == '__main__':
  test()
