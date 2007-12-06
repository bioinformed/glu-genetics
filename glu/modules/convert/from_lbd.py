# -*- coding: utf-8 -*-
'''
File:          from_lbd.py

Authors:       Brian Staats (staatsb@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-15

Abstract:      Parses an Illumina Locus by DNA report file and outputs a
               genotype matrix or triple stream.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   operator                  import itemgetter
from   itertools                 import islice,chain,izip,groupby

from   numpy                     import array,zeros

from   glu.lib.utils             import izip_exact
from   glu.lib.fileutils         import autofile,hyphen
from   glu.lib.sections          import read_sections
from   glu.lib.sequence          import norm_snp_seq,complement_base

from   glu.lib.genolib.locus     import Genome, Nothing
from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.io        import save_genostream
from   glu.lib.genolib.genoarray import model_from_alleles


def load_lbd_file(filename):
  def parse_gentrain():
    row = data.next()
    assert row[1] == 'Gentrain Scores'
    assert row[2] == ''
    return map(float,row[3:])

  def parse_loci():
    row = data.next()
    assert row[1] == 'locusNames'
    assert row[2] == ''
    return row[3:]

  def skiprows(n):
    for i in xrange(n):
      data.next()

  def sample_generator():
    for key,rows in groupby(data,itemgetter(0,1,2,3,4)):
      rows = list(rows)
      assert len(rows) == 2
      genos,scores = rows
      assert  genos[6] == 'calls'
      assert scores[6] == 'Score_Call'

      sampleid = genos[0]

      # FIXME: Why do we have blanks?  Is this a local artifact or an
      # Illumina one?
      if sampleid == 'blank':
        continue

      genos  = genos[8:]
      scores = map(float,islice(scores,8,None))

      yield sampleid,genos,scores

  data = csv.reader(autofile(filename), dialect='excel')

  skiprows(11)
  gentrain = parse_gentrain()
  skiprows(3)
  loci    = parse_loci()
  skiprows(5)
  samples = sample_generator()

  assert len(gentrain) == len(loci)

  return loci,gentrain,samples


def load_abmap(file_or_name,skip=0):
  '''
  Creates a dictionary representing a mapping from A and B probes for each
  locus to allele names.
  '''
  lfile = autofile(file_or_name)
  f = csv.reader(lfile, dialect='excel-tab')

  if skip:
    f = islice(f,skip,None)

  for row in f:
    if len(row) < 3 or '' in row[:3]:
      continue

    locus,a,b = row[:3]

    yield locus,(a,b)


def load_illumina_manifest(filename):
  '''
  Load an Illumina assay manifest file, parsing out the various sections
  '''
  ifile = csv.reader(autofile(filename),dialect='excel')
  sections = read_sections(ifile)

  heading,contents = sections.next()

  # OPA manifest
  if heading == 'data':
    heading,contents = sections.next()
    assert heading == 'Heading'

    headings = dict(c[:2] for c in islice(contents,10) if len(c) > 1)
    assert headings['Assay Format'] == 'Golden Gate'

  # Infinium manifest
  elif heading == 'Heading':
    contents = dict(c[:2] for c in contents if len(c) > 1)
    assert contents['Assay Format'] in ('Infinium','Infinium II')

    heading,contents = sections.next()
    assert heading == 'Assay'

  return contents


def find_index(header,headings):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass
  raise ValueError,'Cannot find heading index'


def parse_manifest(manifest,genome,abmap,targetstrand='customer'):
  '''
  Parse an Illumina manifest to obtain mappings from A/B probes to alleles
  relative to a specific genomic strand
  '''
  targetstrand = targetstrand.lower()
  assert targetstrand in ('top','bottom','forward','reverse','customer','anticustomer','design','antidesign')

  manifest    = iter(manifest)
  header      = manifest.next()
  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  snp_idx     = find_index(header,['SNP'])
  chr_idx     = find_index(header,['Chr'])
  loc_idx     = find_index(header,['MapInfo'])
  cstrand_idx = find_index(header,['CustomerStrand'])
  dstrand_idx = find_index(header,['IlmnStrand','Ilmn Strand'])
  topseq_idx  = find_index(header,['TopGenomicSeq'])
  max_idx     = max(topseq_idx,dstrand_idx,cstrand_idx,snp_idx,name_idx,assayid_idx)

  for assay in manifest:
    assay  += ['']*(max_idx-len(assay)+1)
    locus   = assay[name_idx]
    snp     = assay[snp_idx]
    chr     = assay[chr_idx] or None
    loc     = assay[loc_idx]
    cstrand = assay[cstrand_idx].lower()
    dstrand = assay[dstrand_idx].lower()
    topseq  = assay[topseq_idx]

    assert cstrand in ('top','bot')
    assert dstrand in ('top','bot')
    assert (snp[0],snp[2],snp[4]) == ('[','/',']')

    if loc:
      loc = int(loc)

    # Alleles on the design strand
    aa,bb = snp[1],snp[3]

    try:
      tstrand,a,b = norm_snp_seq(topseq)
      tstrand     = tstrand.lower()
      assert tstrand == 'top'
    except ValueError:
      tstrand,a,b = dstrand,aa,bb

    if dstrand!=tstrand:
      aa = complement_base(aa)
      bb = complement_base(bb)

    if (a,b) != (aa,bb):
      raise ValueError('Sequence alleles do not match assay alleles')

    gstrand = assay[assayid_idx].split('_')[2]
    assert gstrand in 'FRU'

    if gstrand != 'U':
      # Get the strand orientation of the design sequence
      # Alleles are forward strand if the tstrand matches the design strand
      # and the design is on the forward strand or the converse of both
      # conditions is true.
      forward = (tstrand != dstrand) ^ (gstrand == 'F')
      strand  = '+' if forward else '-'
    else:
      if targetstrand in ('forward','reverse') and gstrand == 'U':
        raise ValueError("Unknown strand for assay '%s'" % locus)
      strand = Nothing

    flip =    ((targetstrand == 'customer'     and tstrand != cstrand)
           or  (targetstrand == 'anticustomer' and tstrand == cstrand)
           or  (targetstrand == 'design'       and tstrand != dstrand)
           or  (targetstrand == 'antidesign'   and tstrand == dstrand)
           or  (targetstrand == 'top'          and tstrand != 'top'  )
           or  (targetstrand == 'bottom'       and tstrand != 'bot'  )
           or  (targetstrand == 'forward'      and not forward       )
           or  (targetstrand == 'reverse'      and     forward       ))

    if flip:
      a,b = complement_base(a),complement_base(b)
      strand = {Nothing:Nothing,'+':'-','-':'+'}[strand]

    genome.merge_locus(locus, chromosome=chr, location=loc, strand=strand)
    abmap[locus] = (a,b)


def filter_gc(samples, gcthreshold):
  '''
  Filter genotypes by a minimum gc score threshold.  Any that do not meet
  that threshold are set to missing
  '''
  for sampleid,genos,gcscores in samples:
    genos = [ (g if gc>gcthreshold else 'U') for g,gc in izip(genos,gcscores) ]
    yield sampleid,genos,gcscores


def build_models(loci,abmap,genome):
  '''
  Build genotype locus models given a mapping derived from either the
  Illumina manifest or a user-specified A/B map file.
  '''
  genomap={}
  modelcache = {}
  models = []

  for locus in loci:
    a,b = abmap.get(locus, ('A','B') )

    key = tuple(sorted([a,b]))
    model = modelcache.get(key)

    if not model:
      model = modelcache[key] = model_from_alleles(key,max_alleles=2)

    genome.merge_locus(locus, model, fixed=True)
    models.append(model)

    genomap[locus,'A'] = model[a,a]
    genomap[locus,'B'] = model[b,b]
    genomap[locus,'H'] = model[a,b]
    genomap[locus,'U'] = model[None,None]

  return genomap,models


def encode_genotypes(loci, samples, genomap):
  '''
  Encode Illumina's A/B/H/U genotype calls into the appropriate genotype
  objects
  '''
  for sampleid,genos,gcscores in samples:
    genos = [genomap[locus,geno] for locus,geno in izip(loci,genos)]
    yield sampleid,genos


class GCSummary(object):
  def __init__(self, loci, samples):
    self.loci     = loci
    self.samples  = samples

  def __iter__(self):
    self.samplestats = samplestats = []
    locusstats1 = zeros(len(self.loci))
    locusstats2 = zeros(len(self.loci))

    n = len(self.loci)

    for sampleid,genos,gcscores in self.samples:
      gc = array(gcscores,dtype=float)
      samplestats.append( (sampleid, gc.mean(), gc.std()) )

      # Track first two uncentered moments
      locusstats1 += gc
      locusstats2 += gc**2

      yield sampleid,genos,gcscores

    self.locusstats = []
    if samplestats:
      m   = len(samplestats)
      mu  = locusstats1/m

      # Compute std.dev from sqrt(E(X**2) - E(X)**2), with compensation for
      # the inherant numerical problems with the approach
      var = locusstats2/m - mu**2
      var[var < 0] = 0
      std = var**0.5

      self.locusstats = [ (locus,mui,stdi) for locus,mui,stdi,gcmax in izip(self.loci,mu,std) ]


def option_parser():
  import optparse

  usage = 'usage: %prog [options] lbdfile...'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output genotype file name')
  parser.add_option('-F','--outformat',  dest='outformat', metavar='string',
                    help='Output genotype format')
  parser.add_option('-G', '--outgenorepr', dest='outgenorepr', metavar='REP',
                    help='Output genotype representation')
  parser.add_option('-m', '--abmap', dest='abmap', metavar='FILE',
                    help='Mappings from A and B probes to other allele codings')
  parser.add_option('-M', '--manifest', dest='manifest', metavar='FILE',
                    help='Illumina manifest file')
  parser.add_option('-s', '--targetstrand', dest='targetstrand', metavar='T', default='customer',
                    help='Target strand based on Illumina manifest file: top, bottom, forward, '
                         'reverse, customer (default), anticustomer, design, antidesign')
  parser.add_option('-t', '--gcthreshold', dest='gcthreshold', type='float', metavar='N', default=0,
                    help='Genotypes with GC score less than N set to missing')
  parser.add_option('--samplestats', dest='samplestats', metavar='FILE',
                    help='Output per sample average GC statistics to FILE')
  parser.add_option('--locusstats',  dest='locusstats',  metavar='FILE',
                    help='Output per locus average GC statistics to FILE')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  genome = Genome()
  abmap = {}

  if options.manifest:
    print >> sys.stderr, 'Processing Illumina manifest file...',
    manifest = load_illumina_manifest(options.manifest)
    parse_manifest(manifest,genome,abmap,targetstrand=options.targetstrand)
    print >> sys.stderr, 'done.'

  if options.abmap:
    abmap.update(load_abmap(options.abmap))

  args = iter(args)
  loci,gentrain,samples = load_lbd_file(args.next())

  loci = list(loci)
  for arg in args:
    more_loci,more_gentrain,more_samples = load_lbd_file(arg,gcthreshold=options.gcthreshold)

    if list(more_loci) != loci:
      raise RuntimeError,'Genotype headers do not match'

    samples = chain(samples,more_samples)

  if options.gcthreshold > 0:
    samples = filter_gc(samples, options.gcthreshold)

  if options.samplestats or options.locusstats:
    # Summary is both object and new stream
    summary = samples = GCSummary(loci,samples)

  genomap,models = build_models(loci,abmap,genome)
  samples = encode_genotypes(loci, samples, genomap)
  genos = GenomatrixStream(samples, 'sdat', loci=loci, models=models, genome=genome)

  save_genostream(options.output,genos,format=options.outformat,genorepr=options.outgenorepr,hyphen=sys.stdout)

  if options.samplestats:
    out = csv.writer(autofile(options.samplestats,'w'),dialect='tsv')
    out.writerow(['SAMPLE','GC_MEAN','GC_STDDEV'])
    out.writerows( [ (s,'%.2f' % gc, '%.3f' % dev)
                      for s,gc,dev in summary.samplestats ] )

  if options.locusstats:
    out = csv.writer(autofile(options.locusstats,'w'),dialect='tsv')
    out.writerow(['LOCUS','GENTRAIN','GC_MEAN','GC_STDDEV'])
    out.writerows( [ (l,'%.2f' % gt, '%.2f' % gc, '%.3f' % dev)
                      for (l,gc,dev),gt in izip_exact(summary.locusstats,gentrain) ] )


if __name__ == '__main__':
  main()
