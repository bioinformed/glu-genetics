# -*- coding: utf-8 -*-
'''
File:          seq_comp.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import time
import sys
import csv
import traceback

from   glu.lib.sequence import *
from   glu.lib.clust    import align


def norm_seq(seq,n):
  flank = 'N'*n
  seq = flank + seq.upper() + flank
  s = seq.index('[')
  e = seq.index(']')
  b = set(seq[s+1:e])
  b.discard('/')
  l = seq[s-n:s]
  r = seq[e+1:e+n+1]
  return l+to_iupac(b)+r


def seq_align(seq1,seq2):
  seq3 = complement(reversed(seq2))
  s1,s2 = align([seq1,seq2])
  m = seq_match(s1,s2)

  if m < len(seq1)/20:
    return s1,s2,m

  s3,s4 = align([seq1,seq3])
  k = seq_match(s3,s4)

  if k < m:
    return s3,s4,k
  else:
    return s1,s2,m


def seq_match(seq1,seq2):
  '''
  Return the count of mismatches between two IUPAC sequences

  @param  seq1: first sequence in iupac representation
  @type   seq1: str
  @param  seq2: second sequence in iupac representation
  @type   seq2: str
  @return     : count of mismatched bases
  @rtype      : int

  >>> s1='hm'
  >>> s2='vb'
  >>> seq_match(s1,s2)
  0
  >>> s1='gc'
  >>> s2='at'
  >>> seq_match(s1,s2)
  2
  '''
  return sequence_match(seq1,seq2).count('')


def re_snp(seq,n):
  l = seq[:n]
  s = seq[n]
  r = seq[n+1:]
  s = '/'.join(iupac[seq[n]])
  return '%s[%s]%s' % (l,s,r)


def merge_bases(b1,b2):
  return to_iupac(iupac[b1].union(b2))


def match_bases(b1,b2):
  b1 = iupac[b1]
  b2 = iupac[b2]
  if b1 == b2:
    return '.'
  else:
    return to_iupac(b1.union(b2))


def conform_left_base(b1,b2):
  b1 = iupac[b1]
  b2 = iupac[b2]
  return len(b2) == 1 and len(b1) != 1


def conform_left(s1,s2):
  return sum(conform_left_base(b1,b2) for b1,b2 in zip(s1,s2))


def find_snp(seq,x):
  '''
  Retrieve a snp at the specified location in the sequence

  @param  seq: nucleotide sequence
  @type   seq: sequence
  @param  seq: location in the sequence
  @type   seq: int
  @return    : snp at the specified location
  @rtype     : str
  '''
  c = 0
  for i,b in enumerate(seq):
    if b != '-':
      if c==x:
        return i
      c+=1

  raise "Cannot find SNP at location %d in sequence %s" % (seq,x)


def _test():
  import doctest
  return doctest.testmod()


def main():
  seq = csv.reader(file('CGF_sequence.csv'))
  #rs  = csv.reader(file('design_2006-02-24.csv'))
  rs  = csv.reader(file('Hap300K_design.csv'))

  seq_header = seq.next()
  rs_header = rs.next()

  rs_name_i  = rs_header.index('SNP_Name')
  rs_seq_i   = rs_header.index('Sequence')
  #rs_score_i = rs_header.index('SNP_Score')

  seq_name_i  = rs_header.index('SNP_Name')
  seq_seq_i   = rs_header.index('Sequence')

  rs = dict( (row[rs_name_i],row) for row in rs )

  x = 40
  results = []
  for i,row1 in enumerate(seq):
    if (i+1)%1000 == 0:
      print >> sys.stderr, '%s: Processing SNP #%d' % (time.asctime(),i+1)

    name = row1[seq_name_i]

    if name not in rs:
      print 'Missing locus:',name
      continue

    row2 = rs[name]

    seq1   = row1[seq_seq_i].replace('-','N')
    seq2   = row2[rs_seq_i].replace('-','N')
    score = ''
    #score  = float(row2[rs_score_i])
    #score = float(row2[score_i])

    try:
      seq1 = norm_seq(seq1,x)
      seq2 = norm_seq(seq2,x)
      s1 = seq1[x]
      s2 = seq2[x]
      seq1,seq2,m = seq_align(seq1, seq2)
      y = find_snp(seq1,x)
      z = find_snp(seq2,x)
      assert seq1[y] in (s1,complement_base(s1)), '%s != %s' % (seq1[y],s1)
      assert seq2[z] in (s2,complement_base(s2)), '%s != %s' % (seq2[z],s2)
      m1 = seq_match(seq1[:y+1], seq2[:y+1])
      m2 = seq_match(seq1[z:], seq2[z:])
      m = m1*m2+m1+m2
      n = conform_left(seq1,seq2)
      k = conform_left(seq2,seq1)
      cseq = [ match_bases(a,b) for a,b in zip(seq1,seq2) ]
      cseq[y] = merge_bases(seq1[y],seq2[y])
      if y != z:
        cseq[z] = merge_bases(seq1[z],seq2[z])
      cseq = ''.join(cseq)
      seq1 = re_snp(seq1,y)
      seq2 = re_snp(seq2,z)
      cseq = re_snp(cseq,y)
      if m or n or k:
        results.append( (name,m,n,k,score,seq1,seq2,cseq) )
    except KeyboardInterrupt:
      raise
    except:
      print >> sys.stderr, 'Bad sequence for snp %s' % name
      print >> sys.stderr, seq1
      print >> sys.stderr, seq2
      print >> sys.stderr
      traceback.print_exc(file=sys.stderr)

  results.sort(key=lambda r: r[1:5],reverse=True)

  for i,(name,m,n,k,score,seq1,seq2,cseq) in enumerate(results):
    print '\t'.join(str(i) for i in [i+1,name,m,n,score,'CGF',seq1])
    print '\t\t\t%d\t\tIllumina\t%s' % (k,seq2)
    print '\t\t\t\t\tConsensus\t%s' % cseq
    print

if __name__ == '__main__':
  _test()
  if 0:
    import profile, pstats
    profile.run('main()', 'snpselect.prof')
    stats = pstats.Stats('snpselect.prof')
    stats.strip_dirs()
    stats.sort_stats('time', 'calls')
    stats.print_stats(25)
  else:
    main()
