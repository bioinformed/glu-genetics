# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert between sequence file formats'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   itertools         import chain

from   Bio.Seq           import Seq
from   Bio               import SeqIO
from   Bio.SeqRecord     import SeqRecord
from   Bio.Alphabet      import single_letter_alphabet

from   glu.lib.fileutils import autofile, hyphen, guess_format, parse_augmented_filename


INPUT_EXTS = ['ace', 'clustal', 'fa', 'fasta', 'fastq', 'fq',
              'fastq-sanger', 'fastq-solexa', 'fastq-illumina', 'genbank',
              'gb', 'ig', 'nexus', 'phd', 'phylip', 'pir', 'stockholm',
              'sff', 'sff-trim', 'swiss', 'qual', 'qseq' ]


OUTPUT_EXTS = ['clustal', 'fa', 'fasta', 'fastq', 'fq',
               'fastq-sanger', 'fastq-solexa', 'fastq-illumina', 'genbank',
               'gb', 'nexus', 'phd', 'phylip', 'stockholm',
               'sff', 'qual' ]


FORMAT_REMAP = {'fa' : 'fasta',
                'fq' : 'fastq',
                'gb' : 'genbank' }


def guess_informat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  f = guess_format(filename, INPUT_EXTS)
  return FORMAT_REMAP.get(f,f)


def guess_informat_list(filenames):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  formats = set( guess_informat(f) for f in filenames )
  formats.discard(None)
  return formats.pop() if len(formats) == 1 else None


def guess_outformat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  f = guess_format(filename, OUTPUT_EXTS)
  return FORMAT_REMAP.get(f,f)


def write_sequence(sequences, filename, outformat=None):
  outformat = outformat or guess_outformat(filename) or 'fasta'
  outformat = FORMAT_REMAP.get(outformat,outformat)
  filename  = parse_augmented_filename(filename, {})
  outfile   = autofile(hyphen(filename,sys.stdout),'wb')

  return SeqIO.write(sequences, outfile, outformat)


def sequence_writer(filename, outformat=None):
  outformat = outformat or guess_outformat(filename) or 'fasta'

  if outformat.endswith('.gz'):
    outformat = outformat[:-3]

  outformat = FORMAT_REMAP.get(outformat,outformat)
  filename  = parse_augmented_filename(filename, {})
  outfile   = autofile(hyphen(filename,sys.stdout),'wb')

  return SeqIO._FormatToWriter[outformat](outfile)


def read_sequence(filename, informat=None):
  informat = informat or guess_informat(filename)
  informat = FORMAT_REMAP.get(informat,informat)
  filename = parse_augmented_filename(filename, {})

  if not informat:
    raise ValueError('Input format must be specified for filename %s' % filename)

  return SeqIO.parse(autofile(filename,'rb'), informat)


def guess_sequence_format(filename, informat=None):
  informat = informat or guess_informat(filename)
  informat = FORMAT_REMAP.get(informat,informat)
  filename = parse_augmented_filename(filename, {})
  return format,filename


def trim_quality(seq,q):
  import numpy as np

  quals  = np.array(seq.letter_annotations['phred_quality'],dtype=int)

  if len(seq):
    scores = np.clip(quals-q,0,1e10).cumsum()
    start  = np.where(scores>0)[0]

    if not len(start):
      return seq[:0]

    start  = start[0]

    stop   = np.where(scores>2)[0]

    if not len(stop):
      return seq[:0]

    stop   = min(stop[-1],scores.argmax())+1
    seq    = seq[start:stop]

  return seq


def trim_sequences(sequences, options):
  trimleft    = options.trimleft    or 0
  trimquality = options.trimquality or 0
  minlen      = options.minlen      or 0

  for seq in sequences:
    if trimleft:
      seq = seq[trimleft:]

    if trimquality:
      seq = trim_quality(seq, trimquality)

    if len(seq)>=minlen:
      yield seq


def drop_short(sequences, minlen):
  for seq in sequences:
    if len(seq)>=minlen:
      yield seq


SOLEXA_SCORE_OFFSET = 64

def QseqIterator(handle, alphabet = single_letter_alphabet):
  """Parse Illumina QSEQ files
  """
  import csv

  q_mapping = {}
  for letter in range(256):
      q_mapping[chr(letter)] = letter-SOLEXA_SCORE_OFFSET

  filtermap = {'0':'Y','1':'N'}

  try:
      filename = handle.name
      parts    = filename.split('/')[-1].split('_')
      read_num = parts[2]
  except (IndexError,AttributeError):
      read_num = '?'

  reader = csv.reader(handle,dialect='excel-tab')

  for row in reader:
      # Row:
      #   0 - Instrument name
      #   1 - Run number
      #   2 - Lane number
      #   3 - Tile number
      #   4 - X coordinate
      #   5 - Y coordinate
      #   6 - index
      #   7 - read number
      #   8 - sequence
      #   9 - qc filtered?
      id = name = '%s:%s:%s:%s:%s:%s:%s %s:%s:%s:%s' % (row[0],row[1],'?',row[2],row[3],row[4],row[5],
                                                        read_num,filtermap[row[10]],0,row[6])

      seq = row[8].replace('.','N')

      record = SeqRecord(Seq(seq, alphabet),
                         id=id, name=name, description='')

      qualities = [q_mapping[letter] for letter in row[9]]
      if qualities and (min(qualities) < 0 or max(qualities) > 62):
          raise ValueError("Invalid character in quality string")

      #Dirty trick to speed up this line:
      #record.letter_annotations["phred_quality"] = qualities
      dict.__setitem__(record._per_letter_annotations,
                       "phred_quality", qualities)
      yield record


import Bio.SeqIO
SeqIO._FormatToIterator['qseq'] = QseqIterator


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('seqfile', nargs='+', help='Input sequence file(s)')

  parser.add_argument('-f', '--informat', metavar='FORMAT',
                    help='Input sequence format.  Formats include: '
                         'ace, clustal, embl, fasta, fastq/fastq-sanger, fastq-solexa, fastq-illumina, '
                         'genbank/gb, ig (IntelliGenetics), nexus, phd, phylip, pir, stockholm, '
                         'sff, sff-trim, swiss (SwissProt), tab (Agilent eArray), qual')
  parser.add_argument('-F', '--outformat', metavar='FORMAT',
                    help='Output sequence format (default=fasta).  As above, except ace, ig, '
                         'pir, sff-trim, swiss.')
  parser.add_argument('--trimleft', type=int, metavar='N', default=0,
                    help='Trim N leading bases')
  parser.add_argument('--trimquality', type=float, metavar='Q', default=0,
                    help='Trim sequences based on quality valuesN leading bases')
  parser.add_argument('--minlen', type=int, metavar='N',
                    help='Drop sequences of length less than N')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if options.informat is None:
    options.informat = guess_informat_list(options.seqfile)

  sequences = chain.from_iterable(read_sequence(filename, options.informat) for filename in options.seqfile)

  if options.trimleft or options.trimquality:
    # Also drops short sequences
    sequences = trim_sequences(sequences, options)
  elif options.minlen:
    sequences = drop_short(sequences, options.minlength)

  count     = write_sequence(sequences, options.output, options.outformat)

  print >> sys.stderr, 'Processed %d records' % count


if __name__=='__main__':
  main()
