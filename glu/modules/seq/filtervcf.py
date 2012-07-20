# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import time
import sys

from   glu.lib.fileutils            import table_reader
from   glu.lib.progressbar          import progress_loop

from   glu.lib.seqlib.vcf           import VCFReader, VCFWriter


def build_segmodel(model,samples,phenofile=None):
  constraintmap = {}

  if phenofile:
    phenos = table_reader(options.phenos)
    for row in phenos:
      if len(row)>1 and row[0] and row[1]:
        sample           = row[0]
        constraint       = row[1]
        cost             = int(row[2] or 1) if len(row)>2 else 1
        score            = int(row[3] or 0) if len(row)>3 else 0
        constraintmap[sample] = constraint,cost,score

  segmodel = []
  maxcost  = 1
  minscore = 0

  if model:
    phenomap = {}
    for sample in samples:
      code = sample[-1]
      if code in 'ACU':
        phenomap[sample] = code

    if model=='aff2':
      cost     = 0
      score    = 0.99
      minscore = 1
      for sample,code in phenomap.items():
        if code in 'AC':
          constraintmap[sample] = 'VAR',cost,score

    elif model=='rec':
      cost  = 1
      score = 0
      for sample,code in phenomap.items():
        if code=='A':
          constraintmap[sample] = 'HOM',cost,score
        elif code=='C':
          constraintmap[sample] = 'VAR',cost,score
        elif code=='U':
          constraintmap[sample] = 'REF',cost,score

    else:
      cost  = 1 if model!='domm1' else 0.999999
      score = 0

      for sample,code in phenomap.items():
        if code in 'AC':
          constraintmap[sample] = 'VAR',cost,score
        elif code=='U':
          constraintmap[sample] = 'REF',1,score

    segmodel = [ (i,)+constraintmap[s] for i,s in enumerate(samples) if s in constraintmap ]

  return maxcost,minscore,segmodel


def check_segmodel(segmodel,genos):
  maxcost,minscore,constraints = segmodel
  cost = score = 0

  for i,constraint,c_cost,c_score in constraints:
    g = genos[i]

    if g=='.':
      a=b='.'
    else:
      a,b = g.split('/')

    if constraint=='VAR':
      match = a not in '0.' or b not in '0.'
    elif constraint=='HOM':
      match = a==b and a not in '0.'
    elif constraint=='HET':
      match = a!=b and a not in '0.' and b not in '0.'
    elif constraint=='REF':
      match = a==b=='0'
    elif constraint=='NOTREF':
      match = a not in '.0' and b not in '.0'
    elif constraint=='MISS':
      match = a==b=='.'
    elif constraint=='NOTMISS':
      match = a!='.' and b!='.'
    else:
      raise ValueError('Unknown constraint: %s' % constraint)

    if match:
      score += c_score
    else:
      cost  += c_cost

    if cost>=maxcost:
      return False

  return score>=minscore


def build_infomap(v):
  infomap = {}
  for inf in v.info:
    if '=' in inf:
      key,value = inf.split('=',1)
    else:
      key,value = inf,''

    infomap[key] = value

  return infomap



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
  parser.add_argument('--maxrefvar', metavar='RATE', type=int, default=0,
                      help='Maximum number of variant outgroup samples (unlimited=0, default)')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  filename = options.variants
  vcf      = VCFReader(options.variants,hyphen=sys.stdin,normalize_indels=False)
  out      = VCFWriter(options.output, vcf.metadata, vcf.samples)
  model    = options.model.lower() if not options.phenos else 'custom'
  segmodel = build_segmodel(model,vcf.samples,options.phenos)

  if segmodel[2]:
    print
    print '-'*80
    print 'MAX COST  =',segmodel[0]
    print 'MIN SCORE =',segmodel[1]
    print 'GENETIC MODEL CONSTRAINTS:'
    for i,code,cost,score in segmodel[2]:
      print '  %s: %s (cost=%f, score=%f)' % (vcf.samples[i],code,cost,score)
    print '-'*80
    print

  includefilter = set()
  for opts in options.includefilter or []:
    includefilter |= set(o.strip().upper() for o in opts.split(','))

  excludefilter = set()
  for opts in options.excludefilter or []:
    excludefilter |= set(o.strip().upper() for o in opts.split(','))

  for v in vcf:
    if options.includeid and not v.names:
      continue

    if options.excludeid and v.names:
      continue

    if includefilter or excludefilter:
      normfilter = set(f.upper() for f in v.filter)
      if includefilter and len(normfilter&includefilter)!=len(includefilter):
        continue

      if excludefilter and (normfilter&excludefilter):
        continue

    if options.maxrefvar>0:
      info = build_infomap(v)
      if int(info.get('REFVAR_OUTGROUP_COUNT',0))>options.maxrefvar:
        continue

    genos = [ g[0] for g in v.genos ]
    if not check_segmodel(segmodel,genos):
      continue

    out.write_locus(v)


if __name__=='__main__':
  main()
