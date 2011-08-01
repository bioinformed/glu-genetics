# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Check concordance between a reference and comparison genotype set'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                  import itemgetter
from   collections               import defaultdict

from   glu.lib.fileutils         import map_reader, table_writer
from   glu.lib.remap             import remap_alleles, remap_category
from   glu.lib.hwp               import hwp_exact_biallelic
from   glu.lib.genolib           import GenotripleStream, load_genostream
from   glu.lib.genolib.transform import load_rename_alleles_file


SAMPLE_HEADER = ['REFKEY','COMPKEY','CONCORD','DISCORD_HET_HET','DISCORD_HET_HOM','DISCORD_HOM_HET',
                 'DISCORD_HOM_HOM','CONCORDANCE_RATE']

LOCUS_HEADER  = SAMPLE_HEADER + ['REF_HWP','COMP_HWP','CONCORD_GENO_PAIR','DISCORD_GENO_PAIR',
                                 'ALLELE_MAP_CATEGORY','ALLELE_MAPS']

def geno_pair_mode(g1,g2):
  return g1.homozygote()*2 + g2.homozygote()


class SampleConcordStat(object):
  def __init__(self):
    self.stats = defaultdict(lambda: [0,0,0,0,0])

  def update(self, refgeno, refsample, compgeno, compsample):
    values = self.stats[refsample,compsample]
    if refgeno==compgeno:
      values[4]    += 1
    else:
      mode = geno_pair_mode(refgeno,compgeno)
      values[mode] += 1


class LocusConcordStat(object):
  def __init__(self):
    self.stats = defaultdict(lambda: defaultdict(int))

  def update(self, refgeno, reflocus, compgeno, complocus):
    self.stats[reflocus,complocus][refgeno,compgeno] += 1


def generate_sample_output(sampleconcord):
  samplestats = []
  for sample,stats in sampleconcord.stats.iteritems():
    concord = stats[4]
    discord = stats[:4]

    samplestat = [sample[0], sample[1], concord]
    samplestat.extend(discord)
    samplestat.append( '%.6f' % (float(concord)/(concord+sum(discord))) )
    samplestats.append(samplestat)

  samplestats.sort(key=itemgetter(7))

  return samplestats


def output_sample_concordstat(filename, sampleconcord):
  f = table_writer(filename)
  f.writerow(SAMPLE_HEADER)
  samplestats = generate_sample_output(sampleconcord)
  totals = [ sum(samplestat[i] for samplestat in samplestats) for i in xrange(2,7) ]
  grand_total = sum(totals)

  if grand_total:
    totals.append( '%.6f' % (float(totals[0])/grand_total) )
  else:
    totals.append('')

  samplestat = ['*','*'] + totals
  samplestats.append(samplestat)
  f.writerows(samplestats)


# FIXME: consolidate with the one in glu.lib.hwp
def count_genos(genos):
  hom1 = hom2 = het = 0
  for g,n in genos.iteritems():
    if g.heterozygote():
      het  = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)
  return hom1,het,hom2


def generate_locus_output(locusconcord,allelemaps):
  locusstats = []
  for locus,stats in locusconcord.stats.iteritems():
    concord   = 0
    discord   = [0,0,0,0]
    refgenos  = {}
    compgenos = {}

    for (g1,g2),n in stats.iteritems():
      refgenos[g1]  = refgenos.get(g1,0)  + n
      compgenos[g2] = compgenos.get(g2,0) + n
      if g1==g2:
        concord       += n
      else:
        mode = geno_pair_mode(g1,g2)
        discord[mode] += n

    locusstat = [locus[0], locus[1], concord]
    locusstat.extend(discord)
    locusstat.append( '%.6f' % (float(concord)/(concord+sum(discord))) )

    refhwp  = hwp_exact_biallelic(*count_genos(refgenos))
    locusstat.append( '%.6f' % refhwp )
    comphwp = hwp_exact_biallelic(*count_genos(compgenos))
    locusstat.append( '%.6f' % comphwp )

    concordgenos = []
    discordgenos = []
    for (g1,g2),n in stats.iteritems():
      cdat = (''.join(g1),''.join(g2),n)
      if g1==g2:
        concordgenos.append(cdat)
      else:
        discordgenos.append(cdat)

    locusstat.append( '|'.join( '%s->%s:%d' % (g1,g2,n) for g1,g2,n in concordgenos ))
    locusstat.append( '|'.join( '%s->%s:%d' % (g1,g2,n) for g1,g2,n in discordgenos ))

    if allelemaps is not None:
      amap = allelemaps.get(locus[1])
      if amap is not None:
        locusstat.append(remap_category(amap))
        locusstat.append('|'.join( '%s->%s' % (a or '',b or '') for a,b in amap.iteritems() if a or b))
      else:
        locusstat.extend(['',''])

    locusstats.append(locusstat)

  locusstats.sort(key=itemgetter(7))

  return locusstats


def output_locus_concordstat(filename, locusconcord, allelemaps):
  f = table_writer(filename)
  f.writerow(LOCUS_HEADER)
  locusstats = generate_locus_output(locusconcord,allelemaps)

  totals = [ sum(locusstat[i] for locusstat in locusstats) for i in xrange(2,7) ]
  grand_total = sum(totals)

  if grand_total:
    totals.append( '%.6f' % (float(totals[0])/sum(totals)) )
  else:
    totals.append('')

  locusstat = ['*','*'] + totals + ['']
  locusstats.append(locusstat)

  f.writerows(locusstats)


def load_reference_genotypes(filename, format, locusset, sampleset):
  # FIXME: needs genome and phenome pass-through
  genos = load_genostream(filename,format=format)
  genos = genos.transformed(include_samples=sampleset, include_loci=locusset)
  return genos


def load_comparison_genotypes(filename, format, locusset, sampleset, lmapfile, smapfile):
  # FIXME: needs genome and phenome pass-through
  genos = load_genostream(filename,format=format)
  genos = genos.transformed(rename_samples=smapfile,   include_samples=smapfile,
                            rename_loci=lmapfile,      include_loci=lmapfile)
  genos = genos.transformed(include_samples=sampleset, include_loci=locusset)
  return genos


def invert_dict(d):
  r = defaultdict(list)
  for key,value in d.iteritems():
    r[value].append(key)
  return dict(r)


def fix_eqsets(refgenos,sampleeq,locuseq):
  # Construct identity sample map if one is not specified
  if sampleeq is None:
    sampleeq = dict( (s,[s]) for s in refgenos.samples )
  else:
    # Otherwise filter equivalence set to include only samples that exist
    samples = set(refgenos.samples)
    for refsample,eqsamples in sampleeq.items():
      eqsamples = [ s for s in sampleeq[refsample] if s in samples ]
      if eqsamples:
        sampleeq[refsample] = eqsamples
      else:
        del sampleeq[refsample]

  # Construct identity locus map if one is not specified
  if locuseq is None:
    locuseq = dict( (l,[l]) for l in refgenos.loci )
  else:
    # Otherwise filter equivalence set to include only loci that exist
    loci = set(refgenos.loci)
    for reflocus,eqloci in locuseq.items():
      eqloci = [ l for l in eqloci if l in loci ]
      if eqloci:
        locuseq[reflocus] = eqloci
      else:
        del locuseq[reflocus]

  return sampleeq,locuseq


def concordance(refgenos,compgenos,sampleeq,locuseq,sampleconcord,locusconcord):
  if refgenos.format not in ('sdat','ldat'):
    refgenos = refgenos.as_ldat()

  refgenos = refgenos.materialize()

  sampleeq,locuseq = fix_eqsets(refgenos,sampleeq,locuseq)

  if refgenos.format == 'ldat':
    samples = dict( (s,i) for i,s in enumerate(refgenos.samples) )
    loci    = dict(refgenos)

    # Assumes missing genotypes are prefiltered
    for sample,locus,compgeno in compgenos:
      if sample not in sampleeq or locus not in locuseq:
        continue

      for refsample in sampleeq[sample]:
        for reflocus in locuseq[locus]:
          refgeno = loci[reflocus][samples[refsample]]

          if refgeno:
            sampleconcord.update(refgeno, refsample, compgeno, sample)
            locusconcord.update(refgeno, reflocus,  compgeno, locus)

  elif refgenos.format == 'sdat':
    loci    = dict( (l,i) for i,l in enumerate(refgenos.loci) )
    samples = dict(refgenos)

    # Assumes missing genotypes are prefiltered
    for sample,locus,compgeno in compgenos:
      if sample not in sampleeq or locus not in locuseq:
        continue

      for refsample in sampleeq[sample]:
        for reflocus in locuseq[locus]:
          refgeno = samples[refsample][loci[reflocus]]

          if refgeno:
            sampleconcord.update(refgeno, refsample, compgeno, sample)
            locusconcord.update(refgeno, reflocus,  compgeno, locus)


def compute_allele_maps(locusconcord):
  for (reflocus,complocus),stats in locusconcord.stats.iteritems():
    concord,bestmap = remap_alleles(stats)
    yield complocus,bestmap


def output_allele_maps(amap,mapfile):
  w = table_writer(mapfile)
  for locus,map in amap:
    if any(1 for a,b in map.iteritems() if a!=b):
      old,new = zip(*map.iteritems())
      w.writerow([locus, ','.join(old), ','.join(new)])


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  usage = 'Usage: %prog [options] reference comparison...'
  parser.add_argument('reference',             help='Reference genotype file')
  parser.add_argument('comparison', nargs='+', help='Comparison genotype file(s)')

  parser.add_argument('-f', '--refformat',  metavar='FILE', default=None,
                     help='The file format for reference genotype data')
  parser.add_argument('-F', '--compformat', metavar='FILE', default=None,
                     help='The file format for other(comparison) genotype data')
  parser.add_argument('-r', '--remap',     metavar='FILE',
                     help='Determine and output the optimal allele mapping based on greatest concordance')
  parser.add_argument('-a', '--allelemap', metavar='FILE',
                     help='A list of loci to remap the comparison data alleles to the reference data alleles')
  parser.add_argument('-o',           dest='sampleout', metavar='FILE',
                     help='Output the concordance statistics by sample to FILE')
  parser.add_argument('-O',           dest='locusout',  metavar='FILE',
                     help='Output the concordance statistics by locus to FILE')
  parser.add_argument('--samplemap',  metavar='FILE',
                     help='Map the sample ids for the comparison data to the set of ids in the sample equivalence map')
  parser.add_argument('--locusmap',   metavar='FILE',
                     help='Map the locus ids for the comparison data to the set of ids in the locus equivalence map')
  parser.add_argument('--sampleeq',   metavar='FILE',
                     help='Equivalence mapping between the sample ids from the comparison data and the reference data')
  parser.add_argument('--locuseq',    metavar='FILE',
                     help='Equivalence mapping between the locus ids from the comparison data and the reference data')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  # Load equivalence maps as many-to-many mappings between final reference
  # names and final comparison names.  Since we match comparison to reference,
  # the mappings must be inverted before use.
  #
  # FIXME: The current code only implements a 1-to-many mapping that does
  #        not represent the true transitive equivalence state, and is
  #        insufficient for reference sets that contain duplicates.  The
  #        solution is to create an equivalence mapping between two
  #        equivalence sets.
  locuseq = eqlocus = None
  if options.locuseq:
    locuseq = map_reader(options.locuseq)
    eqlocus = invert_dict(locuseq)

  sampleeq = eqsample = None
  if options.sampleeq:
    sampleeq = map_reader(options.sampleeq)
    eqsample = invert_dict(sampleeq)

  refgenos  = load_reference_genotypes(options.reference,options.refformat,locuseq,sampleeq)

  compgenos = [ load_comparison_genotypes(filename, options.compformat, eqlocus, eqsample,
                options.locusmap, options.samplemap) for filename in options.comparison ]

  compgenos = GenotripleStream.from_streams(compgenos).transformed(filter_missing=True)

  allelemaps = None
  if options.allelemap:
    allelemaps = load_rename_alleles_file(options.allelemap)
    compgenos = compgenos.transformed(rename_alleles=allelemaps)

  sampleconcord = SampleConcordStat()
  locusconcord  = LocusConcordStat()

  concordance(refgenos,compgenos,eqsample,eqlocus,sampleconcord,locusconcord)

  if options.remap:
    sys.stderr.write('Computing best allele mappings...')
    amap = compute_allele_maps(locusconcord)
    output_allele_maps(amap,options.remap)
    sys.stderr.write('Done.\n')
    # FIXME: Until multipass analysis is implemented, we must stop here.
    return

  if options.sampleout:
    output_sample_concordstat(options.sampleout,sampleconcord)

  if options.locusout:
    output_locus_concordstat(options.locusout,locusconcord,allelemaps)


if __name__=='__main__':
  main()
