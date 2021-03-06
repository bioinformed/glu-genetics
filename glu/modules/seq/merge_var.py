from __future__ import division

__gluindex__  = True
__abstract__  = 'Merge multiple CGI mastervar files together to create a unified list of variant individuals'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os
import sys

from   operator                     import attrgetter
from   itertools                    import groupby, izip

from   glu.lib.utils                import unique
from   glu.lib.imerge               import imerge
from   glu.lib.fileutils            import list_reader, table_writer
from   glu.lib.seqlib.cga           import cga_reader


COMPONENTS = ['CDS','INTRON','DONOR','ACCEPTOR','TSS-UPSTREAM',
              'SPAN5','SPAN3','SPAN','UTR5','UTR3','UTR']


IMPACTS = ['NO-CHANGE','SYNONYMOUS',
           'MISSENSE','NONSENSE','NONSTOP',
           'DELETE','INSERT','DELETE+','INSERT+','FRAMESHIFT',
           'MISSTART','DISRUPT',
           'UNKNOWN-VNC','UNKNOWN-INC','UNKNOWN-TR']


UNKNOWN = set(['UNKNOWN-VNC','UNKNOWN-INC','UNKNOWN-TR','',None])

IMPACT_MAP = {'UNKNOWN-VNC' : '',
              'UNKNOWN-INC' : '',
              'UNKNOWN-TR'  : '',
              None          : ''}


CHROMOSOMES = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
               'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
               'chr20','chr21','chr22','chrX','chrY','chrM']


def mastervar_variants(filename):
  attrs,header,records = cga_reader(filename,extra='individual',skip_ref=True)
  name = os.path.basename(filename)
  name = name.split('.')[0]

  def _records():
    for v in records:
      v.individual = name
      yield v

  return name,_records()


def function_records(gene):
  if not gene:
    return ''
  elif len(gene)==1:
    g = gene[0]
    return '%s:%s:%s' % (g.symbol,g.component,IMPACT_MAP.get(g.impact,g.impact))
  else:
    return ';'.join(unique('%s:%s:%s' % (g.symbol,g.component,IMPACT_MAP.get(g.impact,g.impact)) for g in gene))


def merge_sparse(loci):
  yield ['CHROMOSOME','BEGIN','END','REFERENCE','ALLELES','XREFS','FUNCTION','INDIVIDUALS','GENOTYPES']

  for (chromosome,begin,end),locus in loci:
    locus        = list(locus)
    exemplar     = locus[0]
    allelemap    = dict( [(exemplar.reference or '',0), ('?','?'), ('=','=')] )
    alleles      = []
    xrefs        = []
    function     = []
    individuals  = []
    genotypes    = []

    for var in locus:
      a1,a2 = var.allele1Seq or '',var.allele2Seq or ''

      if a1 not in allelemap:
        alleles.append(a1)
        allelemap[a1] = len(alleles)
        function.append(function_records(var.allele1Gene))
        xrefs.append(var.allele1XRef or '')

      if a2 not in allelemap:
        alleles.append(a2)
        allelemap[a2] = len(alleles)
        function.append(function_records(var.allele2Gene))
        xrefs.append(var.allele2XRef or '')

      individuals.append(var.individual)
      genotypes.append('%s/%s' % (allelemap[a1],allelemap[a2]))

  yield [chromosome,begin,end,
                    exemplar.reference,
                    ','.join(alleles),
                    ','.join(xrefs),
                    ','.join(function),
                    ','.join(individuals),
                    ','.join(genotypes)]


def merge_dense(loci,names):
  yield ['CHROMOSOME','BEGIN','END','REFERENCE','ALLELES','XREFS','FUNCTION'] + names

  n      = len(names)
  indmap = dict( (s,i) for i,s in enumerate(names) )

  for (chromosome,begin,end),locus in loci:
    locus        = list(locus)
    exemplar     = locus[0]
    allelemap    = dict( [(exemplar.reference or '',0), ('?','?'), ('=','=')] )
    alleles      = []
    xrefs        = []
    function     = []
    genotypes    = ['0/0']*n

    for var in locus:
      a1,a2 = var.allele1Seq or '',var.allele2Seq or ''

      if a1 not in allelemap:
        alleles.append(a1)
        allelemap[a1] = len(alleles)
        function.append(function_records(var.allele1Gene))
        xrefs.append(var.allele1XRef or '')

      if a2 not in allelemap:
        alleles.append(a2)
        allelemap[a2] = len(alleles)
        function.append(function_records(var.allele2Gene))
        xrefs.append(var.allele2XRef or '')

      genotypes[ indmap[var.individual] ] = '%s/%s' % (allelemap[a1],allelemap[a2])

    yield [chromosome,begin,end,
                      exemplar.reference,
                      ','.join(alleles),
                      ','.join(xrefs),
                      ','.join(function)] + genotypes


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', nargs='+', help='Input variant files')

  parser.add_argument('-F', '--outformat', metavar='FORMAT', default='dense',
                      help='Output format: sparse or dense (default)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser     = option_parser()
  options    = parser.parse_args()

  refmap     = dict( (r,i) for i,r in enumerate(CHROMOSOMES) )
  variants   = []
  names      = []

  for filename in options.variants:
    name,vars = mastervar_variants(filename)
    names.append(name)
    variants.append(vars)

  names.sort()

  def mergekey(v):
    return refmap[v.chromosome],v.begin,v.end

  variants   = imerge(variants, key=mergekey)

  out        = table_writer(options.output,hyphen=sys.stdout)

  keyfunc    = attrgetter('chromosome','begin','end')
  loci       = groupby(variants,keyfunc)

  outformat  = options.outformat.lower()
  if outformat=='dense':
    rows = merge_dense(loci,names)
  elif outformat=='sparse':
    rows = merge_sparse(loci)
  else:
    raise ValueError('Unknown output format: %s' % outformat)

  for row in rows:
    out.writerow(row)


if __name__=='__main__':
  if 1:
    main()
  else:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof,stream=sys.stderr)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
