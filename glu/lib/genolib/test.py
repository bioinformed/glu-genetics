# -*- coding: utf-8 -*-
'''
File:          test.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Run genolib doctest suite

Requires:      Python 2.5, glu

Revision:      $Id$
'''

def test():
  import os
  import doctest

  from   glu.lib.genolib import binary,bitarray,genoarray,io,merge,reprs,streams
  import streams as local_streams

  if not os.path.samefile(streams.__file__,local_streams.__file__):
    raise ImportError('Your PYTHONPATH is not set correctly to find this GLU tree')

  for module in binary,bitarray,genoarray,io,merge,reprs,streams:
    doctest.testmod(module,report=False)
  doctest.master.summarize()


if __name__ == '__main__':
  test()
