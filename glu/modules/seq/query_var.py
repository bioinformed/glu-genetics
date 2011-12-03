# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Query merged CGI variants'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                     import itemgetter

from   glu.lib.utils                import unique
from   glu.lib.fileutils            import table_writer, table_reader, autofile, hyphen
from   glu.lib.progressbar          import progress_loop

from   glu.lib.genedb               import open_genedb
from   glu.lib.genedb.queries       import query_gene_by_name
from   glu.lib.seqlib.vannotator    import VariantAnnotator
from   glu.lib.seqlib.cgfvariants   import CGFVariants

import pysam


COMPONENTS = ['CDS','INTRON','DONOR','ACCEPTOR','TSS-UPSTREAM',
              'SPAN5','SPAN3','SPAN','UTR5','UTR3','UTR']


IMPACTS = ['NO-CHANGE','SYNONYMOUS',
           'MISSENSE','NONSENSE','NONSTOP',
           'DELETE','INSERT','DELETE+','INSERT+','FRAMESHIFT',
           'MISSTART','DISRUPT',
           'UNKNOWN-VNC','UNKNOWN-INC','UNKNOWN-TR']


def as_set(s):
  if s is None:
    return s
  return set(s)


def split_list(s):
  if s is None:
    return None
  return s.split(',')


def query_variants(filename,genedb,genes=None,regions=None,
                                   includeind=None,excludeind=None,
                                   includeknown=None,excludeknown=None,
                                   includeimpact=None,excludeimpact=None,
                                   includecomponent=None,excludecomponent=None):

  vars     = table_reader(filename)
  header   = next(vars)
  vars.close()

  vars     = pysam.Tabixfile(filename)

  queries    = []

  for gene in genes or []:
    g = query_gene_by_name(genedb,gene)
    queries.append( (gene, '%s:%d-%d' % (g[2],g[3]+1,g[4])) )

  for r in regions or []:
    queries.append( (r,r) )

  if not queries:
    queries.append( ('','') )

  includeind   = as_set(includeind)
  excludeind   = as_set(excludeind)
  filterind    = includeind is not None or excludeind is not None

  includeknown = as_set(includeknown)
  excludeknown = as_set(excludeknown)
  filterknown  = includeknown is not None or excludeknown is not None

  includecomp  = as_set(includecomponent)
  excludecomp  = as_set(excludecomponent)
  filtercomp   = includecomp is not None or excludecomp is not None

  includeimpact  = as_set(includeimpact)
  excludeimpact  = as_set(excludeimpact)
  filterimpact   = includeimpact is not None or excludeimpact is not None

  # CHROMOSOME      BEGIN   END     REFERENCE       ALLELES XREFS   GENES   COMPONENTS      IMPACTS INDIVIDUALS     GENOTYPES
  xref_idx = header.index('XREFS')
  comp_idx = header.index('COMPONENTS')
  imp_idx  = header.index('IMPACTS')
  ind_idx  = header.index('INDIVIDUALS')

  yield ['REGION']+header

  for rname,region in queries:
    for var in vars.fetch(region=region):
      var = var.rstrip().split('\t')

      if filterind:
        inds = var[ind_idx].split(',')
        if includeind is not None and     includeind.isdisjoint(inds):
          continue
        if excludeind is not None and not excludeind.isdisjoint(inds):
          continue

      if filterknown:
        xrefs = [ r.split(':')[0].split('.')[0] for x in var[xref_idx].split(',') for r in x.split(';') ]
        if includeknown is not None and includeknown.isdisjoint(xrefs):
          continue
        if excludeknown is not None and excludeknown.issuperset(xrefs):
          continue

      if filtercomp:
        comps = var[comp_idx].split(',')
        if includecomp is not None and includecomp.isdisjoint(comps):
          continue
        if excludecomp is not None and excludecomp.issuperset(comps):
          continue

      if filterimpact:
        impacts = var[imp_idx].split(',')
        if includeimpact is not None and includeimpact.isdisjoint(impacts):
          continue
        if excludeimpact is not None and excludeimpact.issuperset(impacts):
          continue

      yield rname,var


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  #parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
  #                    help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--cgfvariants',   metavar='NAME',
                      help='CGFvariant database annotation')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')

  filters = parser.add_argument_group('Variant filtering (comma delimited values)')

  filters.add_argument('--gene', metavar='NAME',
                      help='Gene symbols or IDs')
  filters.add_argument('--region', metavar='REGION',
                      help='Region chr:begin-end')
  filters.add_argument('--includeind', metavar='NAME',
                      help='Include variants variant in individual')
  filters.add_argument('--excludeind', metavar='NAME',
                      help='Exclude variants variant in individual')
  filters.add_argument('--includeknown', metavar='NAME',
                      help='Include variants known in database')
  filters.add_argument('--excludeknown', metavar='NAME',
                      help='Exclude variants known in database')
  filters.add_argument('--includecomponent', metavar='NAME',
                      help='Gene component to include')
  filters.add_argument('--excludecomponent', metavar='NAME',
                      help='Gene component to exclude')
  filters.add_argument('--includeimpact', metavar='NAME',
                      help='Functional impact to include')
  filters.add_argument('--excludeimpact', metavar='NAME',
                      help='Functional impact to exclude')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  #vs     = VariantAnnotator(options.genedb, options.reference)
  #cv     = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  genedb  = open_genedb(options.genedb)
  out     = table_writer(options.output,hyphen=sys.stdout)

  vars = query_variants(options.variants, genedb,
                        genes            = split_list(options.gene),
                        regions          = split_list(options.region),
                        includeind       = split_list(options.includeind),
                        excludeind       = split_list(options.excludeind),
                        includeknown     = split_list(options.includeknown),
                        excludeknown     = split_list(options.excludeknown),
                        includeimpact    = split_list(options.includeimpact),
                        excludeimpact    = split_list(options.excludeimpact),
                        includecomponent = split_list(options.includecomponent),
                        excludecomponent = split_list(options.excludecomponent))

  header = next(vars)
  out.writerow(header)

  for rname,var in vars:
    out.writerow([rname]+var)


if __name__=='__main__':
  main()
