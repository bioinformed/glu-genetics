# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import time
import sys

from   collections                  import defaultdict

from   pysam                        import Fastafile

from   glu.lib.fileutils            import table_writer, table_reader, autofile, hyphen, map_reader
from   glu.lib.progressbar          import progress_loop

from   glu.lib.seqlib.edits         import reduce_match
from   glu.lib.seqlib.intervaltree  import IntervalTree


def build_segmodel(model,samples,phenofile=None):
  phenomap = None
  if phenofile:
    phenos = table_reader(options.phenos)
    phenomap = {}
    for row in phenos:
      if len(row)>1 and row[0] and row[1]:
        sample           = row[0]
        constraint       = row[1]
        cost             = int(row[2] or 1) if len(row)>2 else 1
        phenomap[sample] = constraint,cost

  segmodel   = []

  if model:
    if phenomap is None:
      phenomap = {}
      for sample in samples:
        code = sample[-1]
        if code in 'ACU':
          phenomap[sample] = code

    if model=='rec':
      cost = 1
      for sample,code in phenomap.items():
        if code=='A':
          phenomap[sample] = 'HOM',cost
        elif code=='C':
          phenomap[sample] = 'VAR',cost
        elif code=='U':
          phenomap[sample] = 'REF',cost

    else:
      cost = 1 if model!='domm1' else 0.999999

      for sample,code in phenomap.items():
        if code in 'AC':
          phenomap[sample] = 'VAR',cost
        elif code=='U':
          phenomap[sample] = 'REF',1

    segmodel = [ (i,)+phenomap[s] for i,s in enumerate(samples) if s in phenomap ]

  return segmodel


def check_segmodel(segmodel,genos,maxfail=1):
  fail = 0
  for i,constraint,cost in segmodel:
    g = genos[i]
    if constraint=='VAR' and g not in ('0/1','1/0','1/1'):
      fail += cost
    elif constraint=='HOM' and g!='1/1':
      fail += cost
    elif constraint=='HET' and g not in ('0/1','1/0'):
      fail += cost
    elif constraint=='REF' and g!='0/0':
      fail += cost
    elif constraint=='NOTREF' and g=='0/0':
      fail += cost
    elif constraint=='MISS' and g!='./.':
      fail += cost

    if fail>=maxfail:
      return False

  return True


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('--includeid', action='store_true',
                   help='Include variants with non-missing ID.')
  parser.add_argument('--excludeid', action='store_true',
                   help='Exclude variants with non-missing ID.')

  parser.add_argument('--includefilter', action='append', metavar='TERMS',
                   help='Include variants that contain all of the specified filter terms.')
  parser.add_argument('--excludefilter', action='append', metavar='TERMS',
                   help='Exclude variants that contain any of the the specified filter terms.')

  parser.add_argument('--phenos', metavar='NAME',
                      help='Mapping from individual ID to phenotype state')
  parser.add_argument('--model',   metavar='NAME', default='',
                      help='Segregation model (dom,rec)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  filename = options.variants
  variants = table_reader(filename,hyphen=sys.stdin)
  out      = autofile(hyphen(options.output,sys.stdout),'w')
  model    = options.model.lower() if not options.phenos else 'custom'

  includefilter = set()
  for opts in options.includefilter or []:
    includefilter |= set(o.strip().upper() for o in opts.split(','))

  excludefilter = set()
  for opts in options.excludefilter or []:
    excludefilter |= set(o.strip().upper() for o in opts.split(','))

  for row in variants:
    if not row or row[0].startswith('##'):
      out.write('\t'.join(row))
      out.write('\n')
      continue

    elif row[0].startswith('#'):
      out.write('\t'.join(row))
      out.write('\n')
      header = list(row)
      header[0] = header[0][1:]

      samples    = [ s.split('.')[0] for s in header[9:] ]
      segmodel   = build_segmodel(model,samples,options.phenos)

      if segmodel:
        print
        print '-'*80
        print 'GENETIC MODEL CONSTRAINTS:'
        for i,code,cost in segmodel:
          print '  %s: %s (cost=%f)' % (samples[i],code,cost)
        print '-'*80
        print

      continue

    chrom = row[0]
    end   = int(row[1])
    start = end-1
    names = row[2]

    if names=='.':
      names = []
    else:
      names = names.split(',')

    if options.includeid and not names:
      continue

    if options.excludeid and names:
      continue

    ref   = row[3]
    var   = row[4].split(',')[0]

    filter = row[6]

    if filter=='.':
      filter = []
    else:
      filter = filter.split(';')

    if includefilter or excludefilter:
      normfilter = set(f.upper() for f in filter)
      if includefilter and len(normfilter&includefilter)!=len(includefilter):
        continue

      if excludefilter and (normfilter&excludefilter):
        continue

    info = row[7].split(';')
    genos = [ g.split(':')[0] for g in row[9:] ]

    if not check_segmodel(segmodel,genos):
      continue

    out.write('\t'.join(row))
    out.write('\n')


if __name__=='__main__':
  main()
