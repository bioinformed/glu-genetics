# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute the genotype correlation between a reference and comparison genotype set'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools                 import izip

from   glu.lib.utils             import Counter
from   glu.lib.fileutils         import table_writer
from   glu.lib.progressbar       import progress_loop
from   glu.lib.association       import estimate_maf

from   glu.lib.genolib           import load_genostream, geno_options
from   glu.lib.genolib.transform import GenoTransform
from   glu.lib.genolib.ld        import count_haplotypes, estimate_ld


def load(options,filename1,filename2):
  sys.stderr.write('Loading %s...\n' % filename1)

  trans    = GenoTransform.from_options(options)

  genos1   = load_genostream(filename1,format=options.informat,genorepr=options.ingenorepr,
                                       genome=options.loci,phenome=options.pedigree,
                                       transform=options).as_ldat().materialize()

  loci1    = set(genos1.loci)
  samples1 = set(genos1.samples)

  if trans.loci.include is None:
    trans.loci.include  = loci1
  else:
    trans.loci.include &= loci1

  if trans.samples.include is None:
    trans.samples.include  = samples1
  else:
    trans.samples.include &= samples1

  sys.stderr.write('Loading %s...\n' % filename2)

  genos2   = load_genostream(filename2,format=options.informat,genorepr=options.ingenorepr,
                                       genome=options.loci,phenome=options.pedigree,
                                       transform=options).as_ldat()

  loci2    = set(genos2.loci)
  samples2 = set(genos2.samples)

  loci     = sorted(loci1&loci2)
  samples  = sorted(samples1&samples2)

  genos1   = genos1.transformed(includeloci=loci,      orderloci=loci,
                                includesamples=samples,ordersamples=samples)

  genos2   = genos2.transformed(includeloci=loci,      orderloci=loci,
                                includesamples=samples,ordersamples=samples)

  return genos1,genos2


def option_parser():
  import optparse

  usage = 'usage: %prog [options] group1 group2'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True,filter=True)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='output table file name')
  parser.add_option('-P', '--progress', dest='progress', action='store_true',
                    help='Show analysis progress bar, if possible')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    sys.exit(2)

  # Build population labels
  genos1,genos2 = load(options,args[0],args[1])
  genos         = izip(genos1,genos2)

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['LOCUS','SAMPLES','MISSING1','MAF1','MISSING2','MAF2','CONCORDANCE','COMPS','R2','HAPLOS'])

  if options.progress:
    genos = progress_loop(genos, length=len(genos1.loci), units='loci')

  for (lname1,locus1),(lname2,locus2) in genos:
    assert lname1==lname2

    counts1  = Counter(locus1)
    counts2  = Counter(locus2)

    maf1     = estimate_maf(counts1)
    missing1 = counts1.get( (None,None), 0 ) / len(locus1)
    maf2     = estimate_maf(counts2)
    missing2 = counts2.get( (None,None), 0 ) / len(locus2)

    haps     = count_haplotypes(locus1, locus2)
    r2,dp    = estimate_ld(*haps)

    concord  = comps = 0
    for a,b in izip(locus1,locus2):
      if a and b:
        if a==b:
          concord += 1
        comps += 1

    concord = concord/comps if comps else 1

    out.writerow([lname1,len(locus1),missing1,maf1,missing2,maf2,concord,comps,r2,str(haps)])


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  main()
