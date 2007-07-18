# -*- coding: utf-8 -*-
'''
File:          sequence.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2007-02-06

Abstract:      Implements sequence alphabet, manipulation, and mapping of
               top/bottom sequence normalization and A/B SNP allele
               canonical mappings based on Illumina and dbSNP standards
               documented at ftp://ftp.ncbi.nih.gov/snp/database/Illumina_top_bot_strand.note.txt

Requires:

Revision:      $Id: sequence.py 494 2007-02-09 14:14:35Z jacobske $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv

from   itertools             import izip

comp  = {'A':'T','C':'G','T':'A','G':'C','-':'-','.':'.'}
iupac = {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'CG','W':'AT','K':'GT',
         'M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT','X':'ACGT', '.':'.', '-':'-'}

iupac_back = dict( (b,a) for a,b in iupac.items() if a!='X' )
iupac_back.update( (a,a) for a,b in iupac.items() )
iupac = dict( (a,set(b)) for a,b in iupac.items() )

def to_iupac(b):
  b = ''.join(''.join(iupac[c.upper()]) for c in b)
  b = set(b)
  b.discard('-')
  if not b:
    return '-'
  b = list(b)
  b.sort()
  return iupac_back[''.join(b)]


def complement_base(b):
  b = iupac[to_iupac(b.upper())]
  return to_iupac([ comp[c] for c in b ])


def complement(s):
  return ''.join([ complement_base(b) for b in s ])


def base_match(b1,b2):
  return iupac[b1.upper()].intersection(iupac[b2.upper()])


def sequence_match_raw(s1, s2):
  return [ base_match(b1,b2) for b1,b2 in zip(s1,s2) ]


def sequence_match(s1, s2):
  return [ ''.join(b) for b in sequence_match_raw(s1, s2) ]


def string_match(s1,s2):
  return sum(1 for c1,c2 in zip(s1,s2) if c1.upper() != c2.upper())


def top_or_bottom_comp(a,b):
  '''
  Return top or bottom strand without adjustment for allele order

  >>> top_or_bottom_comp('A','G')
  'Top'
  >>> top_or_bottom_comp('G','A')
  'Bot'
  >>> top_or_bottom_comp('C','T')
  'Bot'
  >>> top_or_bottom_comp('T','C')
  'Top'
  >>> top_or_bottom_comp('T','G')
  'Top'
  >>> top_or_bottom_comp('A','G')
  'Top'
  >>> top_or_bottom_comp('C','G')
  'Unknown'
  >>> top_or_bottom_comp('G','C')
  'Unknown'
  >>> top_or_bottom_comp('T','A')
  'Unknown'
  >>> top_or_bottom_comp('A','T')
  'Unknown'
  >>> top_or_bottom_comp('A','A')
  'Unknown'
  '''
  strand = 'Unknown'

  top = (a in 'AT') and (b in 'CG')
  bot = (b in 'AT') and (a in 'CG')

  if top and not bot:
    strand = 'Top'
  elif not top and bot:
    strand = 'Bot'

  return strand


def top_or_bottom(a,b):
  '''
  Return top or bottom strand with adjustment for allele order

  >>> top_or_bottom('A','G')
  'Top'
  >>> top_or_bottom('G','A')
  'Bot'
  >>> top_or_bottom('C','T')
  'Top'
  >>> top_or_bottom('T','C')
  'Bot'
  >>> top_or_bottom('T','G')
  'Bot'
  >>> top_or_bottom('A','G')
  'Top'
  >>> top_or_bottom('C','G')
  'Unknown'
  >>> top_or_bottom('G','C')
  'Unknown'
  >>> top_or_bottom('T','A')
  'Unknown'
  >>> top_or_bottom('A','T')
  'Unknown'
  >>> top_or_bottom('A','A')
  'Unknown'
  '''

  strand = top_or_bottom_comp(a,b)

  if a>b and strand == 'Top':
    strand = 'Bot'
  elif a<b and strand == 'Bot':
    strand = 'Top'

  return strand


def norm_snp_seq(seq):
  '''
  >>> norm_snp_seq('TACTC[G/A]GATCGT')
  ('Top', 'A', 'G')
  >>> norm_snp_seq('TACTC[A/G]GATCGT')
  ('Top', 'A', 'G')
  >>> norm_snp_seq('CGAGT[C/T]TTACCT')
  ('Bot', 'T', 'C')
  >>> norm_snp_seq('AGCCC[T/G]GATAGA')
  ('Bot', 'T', 'G')
  >>> norm_snp_seq('TACTT[A/G]GCCAAT')
  ('Top', 'A', 'G')
  >>> norm_snp_seq('TACTCGGCTCT[A/T]GCATTCGATCGT')
  ('Top', 'A', 'T')
  >>> norm_snp_seq('CGTTAGCATAT[A/T]TTTATCGATCGT')
  ('Bot', 'T', 'A')
  >>> norm_snp_seq('ACGATCGATAAA[A/T]ATATGCTAACG')
  ('Top', 'A', 'T')
  >>> norm_snp_seq('TACTGCTCTAT[C/G]GCCCATCAGTC')
  ('Top', 'C', 'G')
  >>> norm_snp_seq('TCCTTTTAATTATTTTTATAGATAATTTTCAACTGCAAGGTGATACATCA[C/T]GTTAATTACTGAATATGGAAAATGTTTGATAAGCCAAAAGGAAAAAACTG')
  ('Bot', 'T', 'C')
  '''
  seq = seq.upper()
  s = seq.index('[')
  e = seq.index(']')

  assert e-s == 4 and seq[s+2] in '/-'

  a = seq[s+1]
  b = seq[e-1]

  if b<a:
    a,b = b,a

  strand = top_or_bottom_comp(a,b)

  for l,r in izip( reversed(seq[:s]),seq[e+1:] ):
    if strand != 'Unknown':
      break
    strand = top_or_bottom_comp(l,r)

  if strand == 'Bot':
    a,b = b,a

  return strand,a,b


if __name__ == '__main__':
  import doctest
  doctest.testmod()
