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

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

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
  '''
  Convert the allele that was passed in from IUPC representation to traditional representation

  @param  b: allele
  @type   b: str
  @return  : allele in traditional representation
  @rtype   : str

  >>> b = 'CGT'
  >>> to_iupac(b)
  'B'
  >>> b = 'ACGT'
  >>> to_iupac(b)
  'N'

  '''
  b = ''.join(''.join(iupac[c.upper()]) for c in b)
  b = set(b)
  b.discard('-')
  if not b:
    return '-'
  b = list(b)
  b.sort()
  return iupac_back[''.join(b)]


def complement_base(b):
  '''
  Complement the base that was passed in

  @param  b: allele
  @type   b: str
  @return  : complemented allele
  @rtype   : str

  >>> b='h'
  >>> complement_base(b)
  'D'
  >>> b='a'
  >>> complement_base(b)
  'T'
  '''
  b = iupac[to_iupac(b.upper())]
  return to_iupac([ comp[c] for c in b ])


def complement(s):
  '''
  Complement the sequence that was passed in

  @param  s: nucleotide sequence
  @type   s: sequence
  @return  : complemented seqeuence
  @rtype   : sequence

  >>> s='ambhg'
  >>> complement(s)
  'TKVDC'
  '''
  return ''.join([ complement_base(b) for b in s ])


def base_match(b1,b2):
  '''
  Return the matched set between two iupac alleles

  @param  b1: first allele in iupac representation
  @type   b1: str
  @param  b2: second allele in iupac representation
  @type   b2: str
  @return   : matched base set
  @rtype    : str

  >>> b1='h'
  >>> b2='v'
  >>> base_match(b1,b2)
  set(['A', 'C'])
  '''
  return iupac[b1.upper()].intersection(iupac[b2.upper()])


def sequence_match_raw(s1, s2):
  '''
  Return the matched set between two iupac alleles

  @param  s1: first sequence in iupac representation
  @type   s1: str
  @param  s2: second sequence in iupac representation
  @type   s2: str
  @return   : list of matched base set
  @rtype    : list

  >>> s1='hm'
  >>> s2='vb'
  >>> sequence_match_raw(s1,s2)
  [set(['A', 'C']), set(['C'])]
  '''

  return [ base_match(b1,b2) for b1,b2 in zip(s1,s2) ]


def sequence_match(s1, s2):
  '''
  Return the matched set between two iupac alleles

  @param  s1: first sequence in iupac representation
  @type   s1: str
  @param  s2: second sequence in iupac representation
  @type   s2: str
  @return   : list of matched base str
  @rtype    : list

  >>> s1='hm'
  >>> s2='vb'
  >>> sequence_match(s1,s2)
  ['AC', 'C']
  '''
  return [ ''.join(b) for b in sequence_match_raw(s1, s2) ]


def string_match(s1,s2):
  '''
  Return the count of mismatched bases

  @param  s1: first sequence in iupac representation
  @type   s1: str
  @param  s2: second sequence in iupac representation
  @type   s2: str
  @return   : count of mismatched bases
  @rtype    : int

  >>> s1='athbc'
  >>> s2='agvba'
  >>> string_match(s1,s2)
  3
  '''

  return sum(1 for c1,c2 in zip(s1,s2) if c1.upper() != c2.upper())


def top_or_bottom_comp(a,b):
  '''
  Return top or bottom strand without adjustment for allele order

  @param    a: first allele
  @type     a: str
  @param    b: second allele
  @type     b: str
  @return    : strand name, 'Top' or 'Bot' or 'Unknown'
  @rtype     : str

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
  >>> top_or_bottom_comp('G','T')
  'Bot'
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

  @param    a: first allele
  @type     a: str
  @param    b: second allele
  @type     b: str
  @return    : strand name, 'Top' or 'Bot' or 'Unknown'
  @rtype     : str

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
  Determine the strand, either 'Top' or 'Bot', for the sequence that was passed in
  Return the alleles of the snp in the alphanumerical order if on 'Top' strand, otherwise the reversed order.

  @param  seq: nucleotide sequence with a snp
  @type   seq: str
  @return    : strand and snp alleles
  @rtype     : tuple of three strs

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
  >>> norm_snp_seq('CAGAAATCATA[C/G]AGAGAACCT')
  ('Top', 'C', 'G')
  >>> norm_snp_seq('GGACCCGCAA[G/A]GAGGGCGCGG')
  ('Top', 'A', 'G')
  >>> norm_snp_seq('CCGCGCCCTC[C/T]TTGCGGGTCC')
  ('Bot', 'T', 'C')
  >>> norm_snp_seq('GGTAGCCTGA[A/T]ACCCCCAAGA')
  ('Top', 'A', 'T')
  >>> norm_snp_seq('TCTTGGGGGT[A/T]TCAGGCTACC')
  ('Bot', 'T', 'A')
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
