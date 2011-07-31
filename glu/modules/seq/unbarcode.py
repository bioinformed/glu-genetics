# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Partition reads by embedded barcode sequences'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   collections             import namedtuple, defaultdict

import numpy as np

from   glu.lib.utils           import pair_generator
from   glu.lib.fileutils       import autofile, hyphen, table_reader
from   glu.lib.seqlib.edits    import levenshtein_distance, hamming_distance
from   glu.lib.seqlib.edits    import levenshtein_sequence, hamming_sequence

from   glu.modules.seq.convert import read_sequence, write_sequence, guess_informat_list


def read_barcodes(filename, trimb=0):
  barcodes      = table_reader(filename)
  header        = next(barcodes)
  barcode_tuple = namedtuple('barcode_tuple', 'id seq')
  barcodes      = ( (id,seq.upper()) for id,seq in barcodes if id and seq )

  if trimb:
    barcodes    = ( (id,seq[:trimb]) for id,seq in barcodes )

  barcodes      = map(barcode_tuple._make, barcodes)

  return barcodes


def read_sequence_and_barcode(filename, format, location, trim5, trim3, trunc, barcode_len):
  seqs = read_sequence(filename, format)

  # FIXME: Add auto-detect

  if location=='454':
    start = 4 if trim5 is None else trim5
    end   = start+barcode_len

    for seq in seqs:
      sseq = str(seq.seq)

      if trunc is not None:
        sseq = sseq[:trunc]

      if trim3:
        sseq = sseq[:-trim3]

      if len(sseq)<end:
        continue

      #assert sseq[:4]=='tcag'
      barcode = sseq[start:end].upper()
      yield barcode,seq

  elif location=='illumina':
    for seq in seqs:
      barcode = seq.id.split('#')[1].split('/')[0]
      yield barcode,seq


def barcode_len(barcodes):
  barcode_lens = set(len(m.seq) for m in barcodes)

  if len(barcode_lens)!=1:
    lens = ','.join(map(str,sorted(barcode_lens)))
    raise ValueError('Barcodes must all be of equal length.  Found: %s' % lens)

  return max(barcode_lens)


def compute_barcode_distances(out,barcodes,distance_metric):
  blen         = barcode_len(barcodes)
  barcodestats = defaultdict(lambda: [0]*(blen+1))

  for b1,b2 in pair_generator(barcodes):
    d = distance_metric(b1.seq,b2.seq)
    barcodestats[b1.id][d] += 1
    barcodestats[b2.id][d] += 1

  for barcode in barcodes:
    stat = barcodestats[barcode.id]
    nz,  = np.nonzero(stat)
    nz   = nz[0] if len(nz) else ' '
    out.write('%s %s %s %s\n' % (barcode.id,barcode.seq,nz,stat))


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('barcodes',             help='Tabular or delimited file of barcode ID and sequences')
  parser.add_argument('sequences', nargs='*', help='Input sequence file(s)')

  parser.add_argument('-f', '--informat', metavar='FORMAT',
                    help='Input sequence format.  Formats include: '
                         'ace, clustal, embl, fasta, fastq/fastq-sanger, fastq-solexa, fastq-illumina, '
                         'genbank/gb, ig (IntelliGenetics), nexus, phd, phylip, pir, stockholm, '
                         'sff, sff-trim, swiss (SwissProt), tab (Agilent eArray), qual')
  parser.add_argument('-F', '--outformat', metavar='FORMAT',
                    help='Output sequence format.  As above, except ace, ig, '
                         'pir, sff-trim, swiss.')
  parser.add_argument('--destdir', default='.',
                    help='Destination directory for output files.  Write to input file directory by default.')
  parser.add_argument('--distance', metavar='TYPE', default='levenshtein',
                    help='Distance metric for barcodes: levenshtein (default) or hamming')
  parser.add_argument('--location', metavar='LOC', default='454',
                    help='Barcode location: 454 (default) or illumina')
  parser.add_argument('--trim5', metavar='N', type=int,
                    help="Trim N 5' bases in reads prior to matching barcode")
  parser.add_argument('--trim3', metavar='N', type=int,
                    help="Trim N 3' bases in reads prior to matching barcode")
  parser.add_argument('--trimb', metavar='N', type=int,
                    help="Trim barcode to N bases")
  parser.add_argument('--trunc', metavar='N', type=int,
                    help="Truncate reads to N bases prior to matching barcode")
  parser.add_argument('--maxerror', metavar='N', type=int,
                    help='Maximum allowed errors')
  parser.add_argument('--mindist', metavar='N', type=int, default=1,
                    help='Minimum edit distance to next best match')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if options.informat is None:
    options.informat = guess_informat_list(options.sequences)

  distance = options.distance.lower()
  if distance and 'levenshtein'.startswith(distance):
    distance_metric = levenshtein_distance
  elif distance and 'hamming'.startswith(distance):
    distance_metric = hamming_distance
  else:
    raise ValueError('Unknown distance metric: %s' % options.distance)

  location = options.location.lower()
  if location not in ('454','illumina'):
    raise ValueError('Unknown barcode location: %s' % options.location)

  barcodes     = read_barcodes(options.barcodes, options.trimb)
  blen         = barcode_len(barcodes)

  out = autofile(hyphen(options.output,sys.stdout),'wb')

  compute_barcode_distances(out,barcodes,distance_metric)

  bad      = 0
  match    = [0]*(blen+1)
  mismatch = [0]*(blen+1)
  hits     = 0
  misses   = 0
  cache    = {}

  barcodestats = defaultdict(lambda: [0]*(blen+1))
  substats = defaultdict(int)
  posstats = defaultdict(int)
  decoded  = defaultdict(list)

  for filename in options.sequences:
    seqs = read_sequence_and_barcode(filename, options.informat, location, options.trim5,
                                               options.trim3, options.trunc, blen)

    for barcode,seq in seqs:
      if len(barcode)!=blen:
        bad += 1
        continue

      if barcode not in cache:
        #distances = ( (m,distance_metric(barcode,m.seq,compress=False)) for m in barcodes )
        #distances = [ (len(ops),m,ops) for m,ops in distances ]
        distances = [ (distance_metric(barcode,m.seq),m) for m in barcodes ]
        distances.sort()

        min_dist,min_barcode   = distances[0]
        min_count              = sum(1 for d,m in distances if d==min_dist)
        next_dist,next_barcode = distances[1]

        if (options.maxerror is not None and min_dist>options.maxerror) or min_count>1 \
          or min_dist+options.mindist>next_dist:
          min_barcode = ''

        #if min_barcode:
        #  print min_barcode,barcode,'distance %d<=%d maxerror' % (min_dist,options.maxerror)

        cache[barcode] = min_dist,min_count,min_barcode
        misses += 1
      else:
        min_dist,min_count,min_barcode = cache[barcode]
        hits += 1

      #print min_barcode,'distance %d<=%d maxerror' % (min_dist,options.maxerror)

      if min_barcode:
        if options.outformat:
          decoded[min_barcode.id].append(seq)

        match[min_dist] += 1
        barcodestats[min_barcode.id][min_dist] += 1
        #if min_dist==1:
        #  op,pos,src,dst = min_ops[0]
        #  if op=='S':
        #    substats[pos+1,src,dst]+= 1
        #    #posstats[pos,min_barcode.seq[pos]]+= 1
      else:
        mismatch[min_dist] += 1

      #print seq.id,barcode,min_dist,min_count,min_barcode

  if options.outformat:
    for mid,sequences in decoded.iteritems():
      filename = '%s/%s.%s' % (options.destdir,mid,options.outformat)
      write_sequence(sequences, filename, options.outformat)
      sequences[:] = []

    decoded.clear()

  if not hits+misses:
    return

  hit_rate = hits/(hits+misses)*100 if hits+misses else 0
  out.write('      BAD BARCODES: %d\n' % bad)
  out.write('   MATCH DISTANCES: %d : %s\n' % (sum(match),match))
  out.write('MISMATCH DISTANCES: %d : %s\n' % (sum(mismatch),mismatch))
  out.write('    CACHE HIT RATE: %.3f%%\n\n' % hit_rate)

  #for pos in sorted(posstats):
  #  print pos,posstats[pos]
  #for i in xrange(1,blen+1):
  #  out.write('%d %s\n' % (i,' '.join('%d' % sum(substats[i,b1,b2] for b2 in 'ACGTN') for b1 in 'ACGTN')))
  #for sub in sorted(substats):
  #  print sub,substats[sub]

  total = 0
  for barcode in barcodes:
    bstats = barcodestats[barcode.id]
    if options.maxerror:
      bstats = bstats[:options.maxerror+1]
    total += sum(bstats)
    out.write('%s %d %s\n' % (barcode.id,sum(bstats),bstats))

  out.write('TOTAL    MATCHES: %d\n' % total)
  out.write('TOTAL MISMATCHES: %d\n' % sum(mismatch))


if __name__=='__main__':
  import doctest
  doctest.testmod()
  main()
