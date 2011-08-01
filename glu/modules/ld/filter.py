# -*- coding: utf-8 -*-

__gluindex__  = True
__program__   = 'TagZilla LD filter'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Sequentially filter a list of SNPs based on an LD threshold'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools               import islice

from   glu.lib.fileutils       import table_reader, table_writer, resolve_column_header
from   glu.lib.genolib         import load_genostream, geno_options
from   glu.lib.genolib.ld      import count_haplotypes, estimate_ld


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('snplist',   help='Tabular or delimited file of SNP names')
  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-r', '--r2threshold', metavar='N', type=float, default=0.80,
                          help='Minimum r-squared threshold (default=0.80)')
  parser.add_argument('-m', dest='maxdist', metavar='BASES', default=200000, type=int,
                    help='Maximum distance in bases between loci to apply LD check.  default=200000')
  parser.add_argument('--lheader', default='Locus',
                    help='Locus header column name or number (default=Locus)')
  parser.add_argument('-L', '--limit', metavar='N', type=int, default=0,
                          help='Filter the top N loci (default=0 for unlimited)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output LD filter results to FILE')

  return parser


def close(loc1,loc2,maxdist):
  chr1,loc1 = loc1.chromosome,loc1.location
  chr2,loc2 = loc2.chromosome,loc2.location
  if None in (chr1,chr2,loc1,loc2):
    return True
  return chr1==chr2 and abs(loc1-loc2) <= maxdist


def main():
  parser = option_parser()
  options= parser.parse_args()

  rows   = table_reader(options.snplist,hyphen=sys.stdin)
  header = rows.next()
  index  = resolve_column_header(header,options.lheader)

  if options.limit:
    rows = islice(rows,options.limit)

  rows   = list(rows)
  snps   = set(row[index] for row in rows)
  genos  = load_genostream(options.genotypes,format=options.informat,genorepr=options.ingenorepr,
                           genome=options.loci,phenome=options.pedigree,
                           transform=options, hyphen=sys.stdin)
  genos  = genos.transformed(include_loci=snps).as_ldat().materialize()
  genome = genos.genome

  out    = table_writer(options.output,hyphen=sys.stdout)

  out.writerow(header+['LDFILTER_RANK','LDFILTER_RANK_TAKEN','LDFILTER_TAKEN','LDFILTER_REASON','LDFILTER_DETAILS'])


  maxdist     = options.maxdist
  r2threshold = options.r2threshold
  taken       = []
  takenset    = set()
  skippedld   = 0
  skippedt    = 0
  missingdata = 0
  genos       = dict(genos)

  for i,row in enumerate(rows):
    locus = row[index]

    if locus in takenset:
      skippedt += 1
      out.writerow(row+[i+1,'','SKIP','ALREADY_TAKEN',''])
      continue

    if locus not in genos:
      takenset.add(locus)
      missingdata += 1
      out.writerow(row+[i+1,'','TAKE','SNP NOT IN REFERENCE DATA',''])
      continue

    geno  = genos[locus]
    loc   = genome.get_locus(locus)
    near  = []
    snpld = []

    for tlocus in taken:
      tgeno = genos[tlocus]
      tloc  = genome.get_locus(tlocus)

      if not close(loc,tloc,maxdist):
        continue

      near.append(tlocus)

      counts    = count_haplotypes(geno,tgeno)
      r2,dprime = estimate_ld(*counts)

      if r2>=r2threshold:
        snpld.append( (tlocus,r2) )

    if not near:
      taken.append(locus)
      takenset.add(locus)
      out.writerow(row+[i+1,len(taken),'TAKE','No SNPs nearby when selected',''])
    elif not snpld:
      taken.append(locus)
      takenset.add(locus)
      out.writerow(row+[i+1,len(taken),'TAKE','Below LD threshold', '%d SNPs nearby' % len(near) ])
    else:
      skippedld += 1
      firstloc,firstr2 = snpld[0]
      out.writerow(row+[i+1,'','SKIP','Above LD threshold',
                               'First: %s r2=%.2f, Others: %d' % (firstloc,firstr2,(len(snpld)-1))])

  print >> sys.stderr, 'Skip Taken    :',skippedt
  print >> sys.stderr, 'Skip LD       :',skippedld
  print >> sys.stderr, 'Taken         :',len(takenset)
  print >> sys.stderr, '  Has ref data:',len(taken)
  print >> sys.stderr, '  No ref data :',missingdata


if __name__ == '__main__':
  main()
