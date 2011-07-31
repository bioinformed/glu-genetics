# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Detect non-Mendelian inheritance among parents and their offspring'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   operator  import itemgetter
from   itertools import repeat,izip

from   glu.lib.utils             import percent
from   glu.lib.fileutils         import table_writer
from   glu.lib.genolib           import load_genostream,geno_options


def parent_offspring_concordance(parent1, parent2, child, locusstats):
  '''
  >>> from glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,build_model
  >>> model = build_model('AB')
  >>> NN,AA,AB,BB = model.genotypes
  >>> def lstats(n): return [ [0,0] for i in xrange(n) ]

  All null:

  >>> p1,p2,c,l = None,None,None,lstats(1)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (0, 0)
  >>> l
  [[0, 0]]

  Child null:

  >>> p1,p2,c,l = [AA],[AB],None,lstats(1)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (0, 0)
  >>> l
  [[0, 0]]

  Parents null:

  >>> p1,p2,c,l = None,None,[AB],lstats(1)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (0, 0)
  >>> l
  [[0, 0]]

  Missings:

  >>> p1,p2,c,l = [NN,AA,NN],[NN,AB,NN],[NN,NN,AB],lstats(3)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (0, 0)
  >>> l
  [[0, 0], [0, 0], [0, 0]]

  Parent concordant:

  >>> p1,p2,c,l = [AA,AA,NN,NN],[NN,NN,AA,AA],[AA,AB,AA,AB],lstats(4)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (4, 4)
  >>> l
  [[1, 1], [1, 1], [1, 1], [1, 1]]

  Homozygous concordant:

  >>> p1,p2,c,l = [AA,BB],[BB,AA],[AB,AB],lstats(2)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (2, 2)
  >>> l
  [[1, 1], [1, 1]]

  Homozygous discordant:

  >>> p1,p2,c,l = [BB,BB,AA,AA],[AA,AA,BB,BB],[AA,BB,AA,BB],lstats(4)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (0, 4)
  >>> l
  [[0, 1], [0, 1], [0, 1], [0, 1]]

  Heterozygous concordant:

  >>> p1,p2,c,l = [AB,AB,AB],[AB,AB,AB],[AA,AB,BB],lstats(3)
  >>> parent_offspring_concordance(p1,p2,c,l)
  (3, 3)
  >>> l
  [[1, 1], [1, 1], [1, 1]]
  '''
  # Must have informative child and at least one parent
  if not child or not (parent1 and parent2):
    return 0,0

  if not parent1:
    parent1 = repeat(None)
  if not parent2:
    parent2 = repeat(None)

  concordant = comparisons = 0
  for p1,p2,c,locusstat in izip(parent1,parent2,child,locusstats):
    # Must have informative child and at least one parent genotype
    if not c or not (p1 or p2):
      continue

    locusstat[1] += 1
    comparisons  += 1

    # If ensure that p1 is not missing
    if not p1 and p2:
      p1,p2 = p2,p1

    a,b = c

    # Check Parent1 -> Offspring case
    if p1 and not p2:
      if a in p1 or b in p1:
        concordant   += 1
        locusstat[0] += 1

    # Check Parent1,Parent2 -> Offspring case
    elif (a in p1 and b in p2) or (b in p1 and a in p2):
      concordant   += 1
      locusstat[0] += 1

  return concordant,comparisons


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genotypes', help='Input genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output Mendelian concordance by sample')
  parser.add_argument('-O', '--locout', metavar='FILE',
                    help='Output Mendelian concordance by locus')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  genos   = load_genostream(options.genotypes,format=options.informat,
                            genorepr=options.ingenorepr,
                            genome=options.loci,
                            phenome=options.pedigree,
                            transform=options)

  if genos.samples:
    pset = set()
    for sample in genos.samples:
      phenos = genos.phenome.get_phenos(sample)
      if phenos.parent1 or phenos.parent2:
        pset.add(sample)
        pset.add(phenos.parent1)
        pset.add(phenos.parent2)
    genos = genos.transformed(includesamples=pset)

  genos = genos.as_sdat().materialize()

  # Index samples
  samples = dict(genos)

  # Initialize statistics
  samplestats = []
  locusstats  = [ [0,0] for i in xrange(len(genos.loci)) ]

  # Check all parent child relationships
  for child in samples:
    phenos  = genos.phenome.get_phenos(child)
    parent1 = phenos.parent1
    parent2 = phenos.parent2

    i,n = parent_offspring_concordance(samples.get(parent1), samples.get(parent2),
                                       samples.get(child), locusstats)
    if i!=n:
      samplestats.append( (child,parent1,parent2,i,n,percent(i,n)) )

  # Build and sort resulting statistics
  locusstats = sorted( (locus,i,n,percent(i,n)) for locus,(i,n) in zip(genos.loci,locusstats) )
  samplestats.sort(key=itemgetter(5))

  # Produce output
  sampleout = table_writer(options.output,hyphen=sys.stdout)
  sampleout.writerow( ['CHILD','PARENT1','PARENT2','CONCORDANT','TOTAL','RATE'] )
  sampleout.writerows(samplestats)

  if options.locout:
    locout = table_writer(options.locout)
    locout.writerow( ['LOCUS','CONCORDANT','TOTAL','RATE'] )
    locout.writerows(locusstats)


if __name__ == '__main__':
  main()
