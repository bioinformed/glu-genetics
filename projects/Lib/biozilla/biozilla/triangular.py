# Routines for computing triangular numbers and their inverses
from math import sqrt


def triangle(n):
  return n*(n+1)//2


def triangle_inv(x):
  return int((sqrt(8*x+1)-1)//2)
