# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Partition reads by embedded barcode sequences'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   collections          import namedtuple, defaultdict

import numpy as np

from   Bio                  import SeqIO

from   glu.lib.utils        import pair_generator
from   glu.lib.fileutils    import autofile, hyphen, table_reader
from   glu.lib.seqlib.edits import levenshtein_sequence, hamming_sequence

from   glu.modules.seq.convert import read_sequence


def read_barcodes(filename):
  barcodes      = table_reader(filename)
  header    = next(barcodes)
  barcode_tuple = namedtuple('barcode_tuple', 'id seq')
  barcodes      = map(barcode_tuple._make, ( (id,seq.upper()) for id,seq in barcodes if id and seq) )
  return barcodes


def read_sequence_and_barcode(filename, location, barcode_len):
  seqs = read_sequence(filename, None)

  # FIXME: Add auto-detect

  if location=='454':
    end = 4+barcode_len
    for seq in seqs:
      sseq = str(seq.seq)
      assert sseq[:4]=='tcag'
      barcode = sseq[4:end].upper()
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


def compute_barcode_distances(barcodes,distance_metric):
  blen         = barcode_len(barcodes)
  barcodestats = defaultdict(lambda: [0]*(blen+1))

  for b1,b2 in pair_generator(barcodes):
      d = len(distance_metric(b1.seq,b2.seq,compress=False))
      barcodestats[b1.id][d] += 1
      barcodestats[b2.id][d] += 1

  for barcode in barcodes:
    stat = barcodestats[barcode.id]
    nz,  = np.nonzero(stat)
    nz   = nz[0] if len(nz) else ' '
    print barcode.id,barcode.seq,nz,stat



def option_parser():
  import optparse

  usage = 'usage: %prog [options] [barcode file] [input files..]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--distance', dest='distance', metavar='TYPE', default='levenshtein',
                    help='Distance metric for barcodes: levenshtein (default) or hamming')
  parser.add_option('--location', dest='location', metavar='LOC', default='454',
                    help='Barcode location: 454 (default) or illumina')
  parser.add_option('--maxerror', dest='maxerror', metavar='N', type='int',
                    help='Maximum allowed errors')
  parser.add_option('--mindist', dest='mindist', metavar='N', type='int', default=1,
                    help='Minimum edit distance to next best match')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output file')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)<1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  distance = options.distance.lower()
  if distance and 'levenshtein'.startswith(distance):
    distance_metric = levenshtein_sequence
  elif distance and 'hamming'.startswith(distance):
    distance_metric = hamming_sequence
  else:
    raise ValueError('Unknown distance metric: %s' % options.distance)

  location = options.location.lower()
  if location not in ('454','illumina'):
    raise ValueError('Unknown barcode location: %s' % options.location)

  barcodes     = read_barcodes(args[0])
  blen         = barcode_len(barcodes)

  compute_barcode_distances(barcodes,distance_metric)

  bad      = 0
  match    = [0]*(blen+1)
  mismatch = [0]*(blen+1)
  hits     = 0
  misses   = 0
  cache    = {}

  barcodestats = defaultdict(lambda: [0]*(blen+1))
  substats = defaultdict(int)
  posstats = defaultdict(int)

  out = autofile(hyphen(options.output,sys.stdout),'wb')
  for filename in args[1:]:
    for barcode,seq in read_sequence_and_barcode(filename, location, blen):

      if len(barcode)!=blen:
        bad += 1
        continue

      if barcode not in cache:
        distances = ( (m,distance_metric(barcode,m.seq,compress=False)) for m in barcodes )
        distances = [ (len(ops),m,ops) for m,ops in distances ]
        distances.sort()

        min_dist,min_barcode,min_ops    = distances[0]
        min_count                   = sum(1 for d,m,ops in distances if d==min_dist)
        next_dist,next_barcode,next_ops = distances[1]

        if (options.maxerror and min_dist>=options.maxerror) or min_count>1 \
          or min_dist+options.mindist>next_dist:
          min_barcode,min_ops='',[]

        #print min_barcode,'%d<%d' % (min_dist,options.maxerror)

        cache[barcode] = min_dist,min_count,min_barcode,min_ops
        misses += 1

      else:
        min_dist,min_count,min_barcode,min_ops = cache[barcode]
        hits += 1

      if min_barcode:
        match[min_dist] += 1
        barcodestats[min_barcode.id][min_dist] += 1
        if min_dist==1:
          op,pos,src,dst = min_ops[0]
          if op=='S':
            substats[pos+1,src,dst]+= 1
            #posstats[pos,min_barcode.seq[pos]]+= 1
      else:
        mismatch[min_dist] += 1

      #print seq.id,barcode,min_dist,min_count,min_barcode

  if not hits+misses:
    return

  hit_rate = hits/(hits+misses)*100 if hits+misses else 0
  print '      BAD BARCODES:',bad
  print '   MATCH DISTANCES:',match
  print 'MISMATCH DISTANCES:',mismatch
  print '    CACHE HIT RATE: %.3f%%' % hit_rate

  #for pos in sorted(posstats):
  #  print pos,posstats[pos]

  for i in xrange(1,blen+1):
    print i,[ sum(substats[i,b1,b2] for b2 in 'ACGTN') for b1 in 'ACGTN' ]

  #for sub in sorted(substats):
  #  print sub,substats[sub]

  for barcode in barcodes:
    print barcode.id,barcode.seq,barcodestats[barcode.id]


if __name__=='__main__':
  import doctest
  doctest.testmod()
  main()
