import os
import time
from Bio.Clustalw import MultipleAlignCL, do_alignment


def align(seqs):
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