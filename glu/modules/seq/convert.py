# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert between sequence file formats'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   itertools         import chain

from   Bio               import SeqIO

from   glu.lib.fileutils import autofile, hyphen, guess_format, parse_augmented_filename


INPUT_EXTS = ['ace', 'clustal', 'fa', 'fasta', 'fastq', 'fq',
              'fastq-sanger', 'fastq-solexa', 'fastq-illumina', 'genbank',
              'gb', 'ig', 'nexus', 'phd', 'phylip', 'pir', 'stockholm',
              'sff', 'sff-trim', 'swiss', 'qual' ]


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


def write_sequence(sequences, filename, outformat):
  outformat = guess_outformat(filename) or outformat or 'fasta'
  outformat = FORMAT_REMAP.get(outformat,outformat)
  filename  = parse_augmented_filename(filename, {})
  outfile   = autofile(hyphen(filename,sys.stdout),'wb')

  return SeqIO.write(sequences, outfile, outformat)


def read_sequence(filename, informat):
  informat = guess_informat(filename) or informat
  informat = FORMAT_REMAP.get(informat,informat)
  filename = parse_augmented_filename(filename, {})

  if not informat:
    raise ValueError('Input format must be specified for filename %s' % filename)

  return SeqIO.parse(autofile(filename,'rb'), informat)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] [input files..]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--informat', dest='informat', metavar='FORMAT', default='sff',
                    help='Input sequence format (default=sff).  Formats include: '
                         'ace, clustal, embl, fasta, fastq/fastq-sanger, fastq-solexa, fastq-illumina, '
                         'genbank/gb, ig (IntelliGenetics), nexus, phd, phylip, pir, stockholm, '
                         'sff, sff-trim, swiss (SwissProt), tab (Agilent eArray), qual')
  parser.add_option('-F', '--outformat', dest='outformat', metavar='FORMAT', default='fasta',
                    help='Output sequence format (default=fasta).  As above, except ace, ig, '
                         'pir, sff-trim, swiss.')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output file')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    sys.exit(2)

  if options.informat is None:
    options.informat = guess_informat_list(args)

  sequences = chain.from_iterable(read_sequence(filename, options.informat) for filename in args)
  count     = write_sequence(sequences, options.output, options.outformat)

  print >> sys.stderr, 'Processed %d records' % count


if __name__=='__main__':
  main()
