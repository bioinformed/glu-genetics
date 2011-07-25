# -*- coding: utf-8 -*-
'''
File:          clust.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Multiple sequence alignment using ClustalW

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import os
import time
from Bio.Clustalw import MultipleAlignCL, do_alignment


def align(seqs):
  '''
  Multiple seqeuence alignment using ClustalW

  @param seqs: list of biological sequences
  @type  seqs: list of strs
  @return    : aligned sequences
  @rtype     : list of strs
  '''
  name = 'tmpFoo%d' % time.time()
  tfile = open(name+'.fasta','w')
  try:
    for i,seq in enumerate(seqs):
      tfile.write('>%d\n' % (i+1))
      tfile.write('%s\n' % seq)
    tfile.flush()

    cline = MultipleAlignCL(tfile.name)
    cline.set_output(name+'.aln')

    alignment = do_alignment(cline)

    return [ rec.seq.tostring() for rec in alignment.get_all_seqs() ]
  finally:
    os.unlink(name + '.fasta')
    os.unlink(name + '.dnd')
    os.unlink(name + '.aln')


def main():
  seq1 = 'GGATGTACTACCAGTTAATTGTTCATGCCTTTAAAAAAAC[C/T]TGTATRCTTCTAAATGTTAATACYATTATATCNTTTCATT'
  seq2 = 'GGATGTACTACCAGTTAATTGTTCATGCCTTTAAAAAAAC[C/T]TGTATGCTTCTAAATGTTAATACCATTATATCTTTCATTG'

  for seq in align([seq1,seq2]):
    print seq


if __name__ == '__main__':
  main()
