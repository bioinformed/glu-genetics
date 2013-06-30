# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   glu.lib.utils                import unique
from   glu.lib.fileutils            import table_writer, table_reader, list_reader, autofile, hyphen
from   glu.lib.progressbar          import progress_loop


from   glu.lib.seqlib.vcf           import VCFReader
from   glu.lib.seqlib.cga           import cga_reader
from   glu.lib.seqlib.vannotator    import VariantAnnotator
from   glu.lib.seqlib.cgfvariants   import CGFVariants
from   glu.lib.seqlib.kaviar        import kaviar_reader
from   glu.lib.seqlib.refvariants   import ReferenceVariants


def flatten_vcf(vcf,records=None):
  filters  = []
  info     = []
  samples  = vcf.samples or []

  for item in vcf.metadata.get('FILTER',[]):
    assert item[:9]=='##FILTER='
    item   = item[9:]
    assert item[0]=='<' and item[-1]=='>'
    item   = item[1:-1]
    filter = item.split(',',1)[0]
    if filter.startswith('ID='):
      filter = filter[3:]
    filters.append(filter)

  for item in vcf.metadata.get('INFO',[]):
    assert item[:7]=='##INFO='
    item   = item[7:]
    assert item[0]=='<' and item[-1]=='>'
    item   = item[1:-1]
    inf    = item.split(',',1)[0]
    if inf.startswith('ID='):
      inf = inf[3:]
    info.append(inf)

  header = ( ['CHROM','REF_START','REF_STOP','IDS','REF','VAR','QUAL']
           + [ f.upper() for f in filters ]
           + [ i.upper() for i in info   ]
           + samples
           + [ '%s_depth' % s for s in samples ] )

  yield header

  if records is None:
    records = iter(vcf)

  for v in records:
    # FORMAT: chrom start end names ref var filter info format genos

    infomap = {}
    for inf in v.info:
      if '=' in inf:
        key,value = inf.split('=',1)
      else:
        key,value = inf,''

      infomap[key] = value

    row = ( [ v.chrom, str(v.start), str(v.end), ','.join(v.names), v.ref, ','.join(v.var), v.qual ]
          + [ 'Y' if f in v.filter else '' for f in filters ]
          + [ infomap.get(i,'') for i in info ]
          + [ g[0] for g in v.genos ]
          + [ '/'.join(g[1].split(',')) if len(g)>1 else '' for g in v.genos ] )

    yield row


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input VCF variant file')

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  vcf      = VCFReader(options.variants,sys.stdin)
  results  = flatten_vcf(vcf)

  out      = table_writer(options.output,hyphen=sys.stdout)
  out.writerows(results)
