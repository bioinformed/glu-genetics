# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Helper functions for genolib'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id: $'


import  numpy as np
cimport numpy as np


def encode_ab_to_snps(np.ndarray abgenos):
  '''
  Recode a genotype array stored as a NumPy array of genotype strings with
  dtype='|S2' to a byte packed 2-bit representation using the following encoding:

  AB  Bits
  --  ----
      0b00
  AA  0b01
  AB  0b10
  BB  0b11

  where bitwise values are stored four to a byte from most to least
  significant bit.
  '''
  if abgenos.dtype.str!='|S2':
    raise ValueError('Invalid genotype input')

  cdef char *ab = abgenos.data
  cdef int    n = abgenos.shape[0] # number of genotypes
  cdef int    m = (n+3)//4         # bytes to allocate
  cdef int    k = n//4             # bytes filled completely
  cdef int    i                    # loop counter
  cdef char val = 0                # Current encoded byte value

  # Allocate empty byte array
  cdef np.ndarray snpgenos = np.empty(m, np.int8)
  cdef char *snps = snpgenos.data

  # Unroll assignment loop for each byte (4 genotypes)
  for i in range(k):
    val = 0

    # Genotype 4*i
    if   ab[0]=='A' and ab[1]=='A': val |= 0x40  # 0b01<<6
    elif ab[0]=='B' and ab[1]=='B': val |= 0xC0  # 0b11<<6
    elif ab[0]=='A' and ab[1]=='B': val |= 0x80  # 0b10<<6

    # Genotype 4*i+1
    if   ab[2]=='A' and ab[3]=='A': val |= 0x10  # 0b01<<4
    elif ab[2]=='B' and ab[3]=='B': val |= 0x30  # 0b11<<4
    elif ab[2]=='A' and ab[3]=='B': val |= 0x20  # 0b10<<4

    # Genotype 4*i+2
    if   ab[4]=='A' and ab[5]=='A': val |= 0x04  # 0b01<<2
    elif ab[4]=='B' and ab[5]=='B': val |= 0x0C  # 0b11<<2
    elif ab[4]=='A' and ab[5]=='B': val |= 0x08  # 0x10<<2

    # Genotype 4*i+3
    if   ab[6]=='A' and ab[7]=='A': val |= 0x01  # 0b01
    elif ab[6]=='B' and ab[7]=='B': val |= 0x03  # 0b11
    elif ab[6]=='A' and ab[7]=='B': val |= 0x02  # 0b10

    # Advance pointer and store result
    ab     += 8
    snps[i] = val

  # Handle last partially filled byte
  # (0-3 remaining genotypes not written by the unrolled loop above)

  i = n%4

  if i>0:
    val = 0

    # Genotype 4*k
    if     ab[0]=='A' and ab[1]=='A': val |= 0x40  # 0b01<<6
    elif   ab[0]=='B' and ab[1]=='B': val |= 0xC0  # 0b11<<6
    elif   ab[0]=='A' and ab[1]=='B': val |= 0x80  # 0b10<<6

    # Genotype 4*k+1
    if i>1:
      if   ab[2]=='A' and ab[3]=='A': val |= 0x10  # 0b01<<4
      elif ab[2]=='B' and ab[3]=='B': val |= 0x30  # 0b11<<4
      elif ab[2]=='A' and ab[3]=='B': val |= 0x20  # 0b10<<4

    # Genotype 4*k+2
    if i>2:
      if   ab[4]=='A' and ab[5]=='A': val |= 0x04  # 0b01<<2
      elif ab[4]=='B' and ab[5]=='B': val |= 0x0C  # 0b11<<2
      elif ab[4]=='A' and ab[5]=='B': val |= 0x08  # 0b10<<2

    # Store result
    snps[k] = val

  return snpgenos


def pair_generator(items):
  '''
  Generator for distinct pairs of items

  @param items: sequence of items
  @type  items: sequence
  @return     : distinct pairs of items
  @rtype      : generator of tuples

  >>> pairs = pair_generator(['A','B','C'])
  >>> for item1,item2 in pairs:
  ...   print item1,item2
  B A
  C A
  C B
  '''
  if not isinstance(items, (list,tuple)):
    items = list(items)
  n = len(items)
  for i in xrange(n):
    for j in xrange(0,i):
      yield items[i],items[j]
