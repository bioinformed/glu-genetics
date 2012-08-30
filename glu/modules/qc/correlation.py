# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Compute the genotype correlation between a reference and comparison genotype set'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools                     import izip

import numpy as np

from   glu.lib.utils                  import Counter
from   glu.lib.fileutils              import table_writer, guess_format
from   glu.lib.progressbar            import progress_loop
from   glu.lib.association            import estimate_maf

from   glu.lib.genolib                import load_genostream, geno_options
from   glu.lib.genolib.transform      import GenoTransform
from   glu.lib.genolib.ld             import count_diplotypes, estimate_ld

from   glu.lib.genolib.formats.wtccc  import load_wtccc_dosage
from   glu.lib.genolib.formats.beagle import load_beagle_dosage


DOSAGE_FORMATS=('beagle','wtccc')


def load_genos(options,filename1,filename2):
  sys.stderr.write('Loading %s...\n' % filename1)

  trans    = GenoTransform.from_options(options)

  genos1   = load_genostream(filename1,format=options.informat,genorepr=options.ingenorepr,
                                       genome=options.loci,phenome=options.pedigree,
                                       transform=options).as_ldat()

  if not genos1.loci or not genos1.samples:
    genos1 = genos1.materialize()

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

  if not genos2.loci or not genos2.samples:
    genos2 = genos2.materialize()

  sys.stderr.write('Merging loci...\n')

  loci2    = set(genos2.loci)
  samples2 = set(genos2.samples)

  loci     = [ l for l in loci1    if l in loci2    ]
  samples  = [ s for s in samples1 if s in samples2 ]

  genos1   = genos1.transformed(includeloci=loci,      orderloci=loci,
                                includesamples=samples,ordersamples=samples,
                                repack=True)

  genos2   = genos2.transformed(includeloci=loci,      orderloci=loci,
                                includesamples=samples,ordersamples=samples,
                                repack=True)

  return genos1,genos2


def genos_to_dosage(locus):
  indices      = locus.indices()
  mask         = indices==0
  dosage       = (indices-1.0)/2.0
  dosage[mask] = np.nan

  return dosage


def load_dosage(options,filename1,filename2):
  sys.stderr.write('Loading %s...\n' % filename1)

  trans    = GenoTransform.from_options(options)

  format1  = guess_format(filename1, ['wtccc'])
  format2  = guess_format(filename2, ['wtccc'])

  if format1 not in DOSAGE_FORMATS:
    genos1   = load_genostream(filename1,format=options.informat,genorepr=options.ingenorepr,
                                         genome=options.loci,phenome=options.pedigree,
                                         transform=options).as_ldat()

    if not genos1.loci or not genos1.samples:
      genos1 = genos1.materialize()

    loci1    = set(genos1.loci)
    samples1 = genos1.samples

    if trans.loci.include is None:
      trans.loci.include  = loci1
    else:
      trans.loci.include &= loci1

  elif format1=='wtccc':
    samples1,dosage1 = load_wtccc_dosage(filename1,'wtccc')
    loci1 = None
  elif format1=='beagle':
    samples1,dosage1 = load_beagle_dosage(filename1,'beagle')
    loci1 = None

  if trans.samples.include is None:
    trans.samples.include  = set(samples1)
  else:
    trans.samples.include &= set(samples1)

  sys.stderr.write('Loading %s...\n' % filename2)

  if format2 not in DOSAGE_FORMATS:
    genos2   = load_genostream(filename2,format=options.informat,genorepr=options.ingenorepr,
                                         genome=options.loci,phenome=options.pedigree,
                                         transform=options).as_ldat()

    if not genos2.loci or not genos2.samples:
      genos2 = genos2.materialize()

    loci2    = set(genos2.loci)
    samples2 = genos2.samples

  elif format2=='wtccc':
    samples2,dosage2 = load_wtccc_dosage(filename2,'wtccc')
    loci2 = None
  elif format2=='beagle':
    samples2,dosage2 = load_beagle_dosage(filename2,'beagle')
    loci2 = None

  sys.stderr.write('Merging loci...\n')

  sset = set(samples1)&set(samples2)

  if loci1 is not None and loci2 is not None:
    samples  = [ s for s in samples1 if s in sset  ]
    loci     = [ l for l in loci1    if l in loci2 ]

    genos1   = genos1.transformed(includeloci=loci,      orderloci=loci,
                                  includesamples=samples,ordersamples=samples,
                                  repack=True)

    genos2   = genos2.transformed(includeloci=loci,      orderloci=loci,
                                  includesamples=samples,ordersamples=samples,
                                  repack=True)

    dosage1 = ( (lname,genos_to_dosage(genos)) for lname,genos in genos1)
    dosage2 = ( (lname,genos_to_dosage(genos)) for lname,genos in genos2)

  elif loci1 is not None:
    smap     = [ (s,i) for i,s in enumerate(samples2) if s in sset ]
    samples,indices2 = zip(*smap)
    indices2 = np.array(indices2,dtype=int)

    dosage2  = [ (lname,dosage[indices2]) for (lname,chrom,loc,a,b,dosage) in dosage2 if lname in loci1 ]
    loci     = [ lname for lname,dosage in dosage2 ]

    genos1   = genos1.transformed(includeloci=loci,      orderloci=loci,
                                  includesamples=samples,ordersamples=samples,
                                  repack=True)

    dosage1 = ( (lname,genos_to_dosage(genos)) for lname,genos in genos1)

  elif loci2 is not None:
    smap     = [ (s,i) for i,s in enumerate(samples1) if s in sset ]
    samples,indices1 = zip(*smap)
    indices1 = np.array(indices1,dtype=int)
    dosage1  = [ (lname,dosage[indices1]) for (lname,chrom,loc,a,b,dosage) in dosage1 if lname in loci2 ]
    loci     = [ lname for lname,dosage in dosage1 ]

    genos2   = genos2.transformed(includeloci=loci,      orderloci=loci,
                                  includesamples=samples,ordersamples=samples,
                                  repack=True)

    dosage2 = ( (lname,genos_to_dosage(genos)) for lname,genos in genos2)

  else:
    raise RuntimeError('Comparison of two dosage files is not yet supported')

  return dosage1,dosage2


def trend_r2(locus1,locus2):
  model1   = locus1.descriptor[0]
  genos1   = model1.genotypes
  model2   = locus2.descriptor[0]
  genos2   = model2.genotypes

  assert genos1[2].heterozygote()
  assert genos2[2].heterozygote()

  indices1 = locus1.indices()
  indices2 = locus2.indices()

  mask = (indices1!=0)&(indices2!=0)

  if mask.sum()<2:
    return 0.

  corr = np.corrcoef(indices1[mask],indices2[mask])

  if len(corr) and np.isfinite(corr[0,1]):
    r2 = corr[0,1]**2
  else:
    r2 = 0.

  return r2


def trend_r2_dips(dips):
   '''
      0: AA CC : 0 0
      1: AA CD : 0 1
      2: AA DD : 0 2
      3: AB CC : 1 0
      4: AB CD : 1 1
      5: AB DD : 1 2
      6: BB CC : 2 0
      7: BB CD : 2 1
      8: BB DD : 2 2
   '''
   xvals = np.array([0,0,0,1,1,1,2,2,2],dtype=int)
   yvals = np.array([0,1,2,0,1,2,0,1,2],dtype=int)
   dips  = np.array(dips, dtype=float)

   n     = dips.sum()

   if not n:
     return 0

   x     = (dips * xvals      ).sum()/n
   xx    = (dips * xvals**2   ).sum()/n
   y     = (dips * yvals      ).sum()/n
   yy    = (dips * yvals**2   ).sum()/n
   xy    = (dips * xvals*yvals).sum()/n

   covxy = xy - x*y
   varx  = xx - x*x
   vary  = yy - y*y

   if varx<=0 or vary<=0:
     return 0

   r2    = covxy**2 / varx / vary

   return r2


def correlation_genos(options,filename1,filename2):
  genos1,genos2 = load_genos(options,filename1,filename2)
  genos         = izip(genos1,genos2)

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['LOCUS','SAMPLES','MISSING1','MAF1','MISSING2','MAF2','CONCORDANCE','COMPS','R2_EM','R2_TREND'])

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

    dips     = count_diplotypes(locus1, locus2)

    #r2_em,dp = estimate_ld(dips)
    r2_em,dp = estimate_ld(locus1,locus2)
    #r2_trend= trend_r2(locus1,locus2)
    r2_trend = trend_r2_dips(dips)

    concord  = comps = 0
    for a,b in izip(locus1,locus2):
      if a and b:
        if a==b:
          concord += 1
        comps += 1

    concord = concord/comps if comps else 1

    out.writerow([lname1,len(locus1),missing1,maf1,missing2,maf2,concord,comps,r2_em,r2_trend])


def correlation_dosage(options,filename1,filename2):

  dosage1,dosage2 = load_dosage(options,filename1,filename2)
  dosage          = izip(dosage1,dosage2)

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(['LOCUS','SAMPLES','MISSING1','MAF1','MISSING2','MAF2','R2'])

  for (lname1,dosage1),(lname2,dosage2) in dosage:
    assert lname1==lname2

    mask1    = np.isfinite(dosage1)
    mask2    = np.isfinite(dosage2)
    mask     = mask1&mask2
    n1       = mask1.sum()
    n2       = mask2.sum()
    maf1     = dosage1[mask1].sum()/n1
    maf2     = dosage2[mask2].sum()/n2
    maf1     = min(maf1,1-maf1)
    maf2     = min(maf2,1-maf2)
    missing1 = 1 - n1/len(dosage1)
    missing2 = 1 - n2/len(dosage2)
    r2       = 0.

    if mask.sum()>=2:
      corr = np.corrcoef(dosage1[mask],dosage2[mask])

      if len(corr) and np.isfinite(corr[0,1]):
        r2 = corr[0,1]**2

    out.writerow([lname1,len(dosage1),missing1,maf1,missing2,maf2,r2])


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('group1', help='Group 1 genotype file')
  parser.add_argument('group2', help='Group 2 genotype file')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='output table file name')
  parser.add_argument('-P', '--progress', action='store_true',
                    help='Show analysis progress bar, if possible')

  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  format1 = guess_format(options.group1, ['wtccc','beagle'])
  format2 = guess_format(options.group2, ['wtccc','beagle'])

  np.seterr(all='ignore')

  if format1 in DOSAGE_FORMATS or format2 in DOSAGE_FORMATS:
    correlation_dosage(options,options.group1,options.group2)
  else:
    correlation_genos(options,options.group1,options.group2)


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  main()
