# -*- coding: utf-8 -*-
'''
File:          ca.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Cochrane-Armitage Trend Test for 2xc tables

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from numpy import asanyarray,arange


def tabletrend(x):
  '''
  Calculate  Cochrane-Armitage Trend Test statistics for 2xc tables

  @param  x: data of 2xc tables
  @type   x: list
  @return  : Cochrane-Armitage Trend Test statistics
  @rtype   : float
  '''

  x = asanyarray(x, dtype=int)

  if len(x.shape) != 2 or x.shape[0] != 2:
    raise ValueError,'tabletrend requires a 2xc table'

  n_i   = x.sum(axis=0)
  n     = n_i.sum()
  r_i   = arange(len(n_i))
  r_bar = float((n_i*r_i).sum())/n
  s2    = (n_i*(r_i-r_bar)**2).sum()
  p1    = float(x[0].sum())/n

  t = (x[0]*(r_i-r_bar)/(p1*(1-p1)*s2)**0.5).sum()

  return t


def test():
  x = [[26,26,23,18,9],[6,7,9,14,23]]
  print tabletrend(x)


if __name__ == '__main__':
  test()
