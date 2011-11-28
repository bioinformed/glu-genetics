# -*- coding: utf-8 -*-

__gluindex__  = True
__program__   = 'TagZilla LD filter'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__abstract__  = 'Sequentially filter a list of SNPs based on an LD threshold'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools                 import islice

from   glu.lib.fileutils         import table_reader, table_writer, resolve_column_header

from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.transform import _intersect_options, _union_options
from   glu.lib.genolib.ld        import count_haplotypes, estimate_ld


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('snplist',   help='Tabular or delimited file of SNP names')
  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-m', dest='maxdist', metavar='BASES', default=200000, type=int,
                    help='Maximum distance in bases between loci to apply LD check.  default=200000')
  parser.add_argument('--lheader', default='Locus',
                    help='Locus header column name or number (default=Locus)')
  parser.add_argument('-L', '--limit', metavar='N', type=int, default=0,
                          help='Filter the top N loci (default=0 for unlimited)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output LD filter results to FILE')
  parser.add_argument('-O', '--detailout', metavar='FILE',
                    help='Output LD filter details to FILE')
  parser.add_argument('--pheader',    default='score p-value',
                          help='p value header column name or number (default=score p-value)')
  parser.add_argument('--p1hit', metavar='P', type=float, default=1.0,
                          help='Filter loci with p value lower than one-hit threshold (default=1.0)')
  parser.add_argument('--p2hit', metavar='P', type=float, default=0.0,
                          help='Filter loci with p value lower than two-hit threshold (default=0.0)')
  parser.add_argument('--r1hit', metavar='N', type=float, default=0.80,
                          help='Minimum one hit r-squared threshold (default=0.80)')
  parser.add_argument('--r2hit', metavar='N', type=float, default=0.80,
                          help='Minimum two hit r-squared threshold (default=0.80)')
  parser.add_argument('--includedesign', metavar='FILE', action='append',
                    help='List of loci that may be taken')
  parser.add_argument('--excludedesign', metavar='FILE', action='append',
                    help='List of loci that may not be taken')

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

  maxdist     = options.maxdist
  take1r2     = options.r1hit
  take2r2     = options.r2hit
  take1p      = options.p1hit
  take2p      = options.p2hit

  if take1p<take2p:
    raise ValueError('Two-hit p-value threshold must be less than the one-hit threshold (%f<%f)' % (take1p,take2p))

  rows   = table_reader(options.snplist,hyphen=sys.stdin)
  header = rows.next()
  lindex = resolve_column_header(header,options.lheader)
  pindex = resolve_column_header(header,options.pheader) if options.pheader else -1

  if options.limit:
    rows = islice(rows,options.limit)

  if pindex>=0 and take1p<1.0:
    snps = set()
    new_rows = []
    for row in rows:
      pval = float(row[pindex] or 1.0)
      if pval>take1p:
        break
      new_rows.append(row)
      snps.add(row[lindex])
    rows = new_rows
  else:
    rows   = list(rows)
    snps   = set(row[lindex] for row in rows)

  genos  = load_genostream(options.genotypes,format=options.informat,genorepr=options.ingenorepr,
                           genome=options.loci,phenome=options.pedigree,
                           transform=options, hyphen=sys.stdin)
  genos  = genos.transformed(include_loci=snps).as_ldat().materialize()
  genome = genos.genome

  out    = table_writer(options.output,hyphen=sys.stdout)

  out.writerow(header+['LDFILTER_RANK','LDFILTER_RANK_TAKEN','LDFILTER_TAKEN','LDFILTER_REASON','LDFILTER_DETAILS'])

  taken       = []
  takenset    = set()
  twohitset   = set()
  twohit      = 0
  skippedld   = 0
  skippedt    = 0
  excluded    = 0
  missingdata = 0
  genos       = dict(genos)

  include = _intersect_options(options.includedesign or [])
  exclude =     _union_options(options.excludedesign or []) or set()

  if include is None:
    include = snps
  else:
    exclude |= snps-include

  r2threshold = take2r2
  pval = 1.0

  for i,row in enumerate(rows):
    locus = row[lindex]

    if pindex>=0:
      pval = float(row[pindex] or 1.0)

      if pval>take1p:
        break
      elif pval>take2p:
        r2threshold = take1r2

    if locus in takenset:
      skippedt += 1
      out.writerow(row+[i+1,'','SKIP','ALREADY_TAKEN',''])
      continue

    if locus in exclude:
      excluded += 1
      out.writerow(row+[i+1,'','SKIP','SNP EXCLUDED FROM CONSIDERATION',''])
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

    hits2  = set(tlocus for tlocus,r2 in snpld if r2>=take2r2)
    hits2 &= twohitset

    if not near:
      taken.append(locus)
      takenset.add(locus)
      out.writerow(row+[i+1,len(taken), 'TAKE', 'No SNPs nearby when selected', ''])
    elif not snpld:
      taken.append(locus)
      takenset.add(locus)
      if pval<=take2p:
        twohitset.add(locus)
      out.writerow(row+[i+1,len(taken), 'TAKE', 'Below LD threshold', '%d SNPs nearby' % len(near)])
    elif hits2:
      twohit += 1
      taken.append(locus)
      takenset.add(locus)
      twohitset -= hits2
      hits2 = ','.join(sorted(hits2))
      out.writerow(row+[i+1,len(taken), 'TAKE', 'Second hit', 'First hit(s): %s' % hits2])
    else:
      skippedld += 1
      firstloc,firstr2 = snpld[0]
      out.writerow(row+[i+1,'','SKIP','Above LD threshold',
                               'First: %s r2=%.2f, Others: %d' % (firstloc,firstr2,(len(snpld)-1))])

  print >> sys.stderr, 'Skip Excluded     :',excluded
  print >> sys.stderr, 'Skip already taken:',skippedt
  print >> sys.stderr, 'Skip LD           :',skippedld
  print >> sys.stderr, 'Take              :',len(takenset)
  print >> sys.stderr, '  Has ref data    :',len(taken)
  print >> sys.stderr, '  No ref data     :',missingdata
  print >> sys.stderr, '  2-hit extra     :',twohit
  print >> sys.stderr, 'Missing 2 hit     :',len(twohitset)

  if options.detailout:
    out = table_writer(options.detailout)

    out.writerow(['LOCUS','STATUS'])
    for locus in sorted(twohitset):
      out.writerow([locus,'Missing second surrogate'])


if __name__ == '__main__':
  main()
