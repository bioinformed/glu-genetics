# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Convert an Illumina Locus by DNA report file into a GLU genotype file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                  import getitem
from   itertools                 import chain, imap, izip, repeat

from   glu.lib.utils             import izip_exact, is_str
from   glu.lib.fileutils         import table_reader, table_writer
from   glu.lib.illumina          import read_Illumina_LBD, create_Illumina_abmap

from   glu.lib.genolib.locus     import Genome, load_genome
from   glu.lib.genolib.phenos    import Phenome, load_phenome
from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.io        import save_genostream, geno_options
from   glu.lib.genolib.transform import GenoTransform
from   glu.lib.genolib.genoarray import build_model, build_descr, GenotypeArray, GenotypeArrayDescriptor


def load_abmap(file_or_name,skip=0):
  '''
  Creates a dictionary representing a mapping from A and B probes for each
  locus to allele names.
  '''
  f = table_reader(file_or_name,skip=skip)

  for row in f:
    if len(row) < 3 or '' in row[:3]:
      continue

    locus,a,b = row[:3]

    yield locus,(a,b)


def filter_gc(samples, gcthreshold):
  '''
  Filter genotypes by a minimum gc score threshold.  Any that do not meet
  that threshold are set to missing
  '''
  for sampleid,genos,gcscores in samples:
    genos[gcscores<=gcthreshold] = 'U'
    yield sampleid,genos,gcscores


def encode_genotypes(loci, samples, genome, phenome):
  '''
  Encode Illumina's A/B/H/U genotype calls to a homogeneous A/B model
  '''
  model = build_model(alleles='AB',genotypes=[('A','A'),('A','B'),('B','B')],max_alleles=2)

  genomap = dict( [('A',model['A','A']),
                   ('B',model['B','B']),
                   ('H',model['A','B']),
                   ('U',model[None,None]) ])

  descr = build_descr(model,len(loci))

  for lname in loci:
    genome.merge_locus(lname, model)

  def _encode(descr):
    for sampleid,genos,gcscores in samples:
      genos = GenotypeArray(descr, imap(getitem, repeat(genomap), genos))
      yield sampleid,genos

  return GenomatrixStream(_encode(descr), 'sdat', loci=loci, models=list(descr),
                                          genome=genome, phenome=phenome, packed=True)


def recode_genotypes(genos, abmap):
  '''
  Recode packed genotypes according to supplied ab-map by writing a new
  descriptor and copying the encoded genotype bitstrings.
  '''

  assert genos.packed

  modelcache = {}
  genome = Genome()

  models = list(genos.models)

  for i,lname in enumerate(genos.loci):
    model     = models[i]
    old_locus = genos.genome.loci[lname]

    if lname in abmap:
      alleles = abmap[lname]
      new_model = modelcache.get(alleles)
      if new_model is None:
        a,b = alleles
        model = build_model(genotypes=[(a,a),(a,b),(b,b)],max_alleles=2)

      models[i] = model

    genome.merge_locus(lname, model, old_locus.chromosome, old_locus.location, old_locus.strand)

  def _recode():
    descr = GenotypeArrayDescriptor(models)
    for sample,row in genos:
      new_row = GenotypeArray(descr)
      new_row.data = row.data
      yield sample,new_row

  return genos.clone(_recode(), models=models, genome=genome, materialized=False)


class GCSummary(object):
  def __init__(self, loci, samples):
    self.loci     = loci
    self.samples  = samples

  def __iter__(self):
    import numpy as np

    self.samplestats = samplestats = []
    locusstats1 = np.zeros(len(self.loci))
    locusstats2 = np.zeros(len(self.loci))
    locuscounts = np.zeros(len(self.loci),dtype=int)

    for sampleid,genos,gcscores in self.samples:
      gc = np.array(gcscores,dtype=float)

      mask = np.isfinite(gc)
      gc[~mask] = 0

      # Track first two uncentered moments
      locusstats1 += gc
      locusstats2 += gc**2

      # Track the number of finite values (in Python, False=0, True=1)
      locuscounts += mask

      gc = gc[mask]
      samplestats.append( (sampleid, gc.mean(), gc.std()) )

      yield sampleid,genos,gcscores

    self.locusstats = []
    if samplestats:
      mu  = locusstats1/locuscounts

      # Compute std.dev from sqrt(E(X**2) - E(X)**2), with compensation for
      # the inherent numerical problems with the approach
      var = locusstats2/locuscounts - mu**2
      var[var < 0] = 0
      std = var**0.5

      self.locusstats = [ (locus,mui,stdi) for locus,mui,stdi in izip(self.loci,mu,std) ]


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('lbdfile', nargs='+', help='LBD file(s)')

  geno_options(parser,filter=True,transform=True,output=True)

  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                      help='Output genotype file name')
  parser.add_argument('-m', '--abmap', metavar='FILE',
                      help='Mappings from A and B probes to other allele codings')
  parser.add_argument('-M', '--manifest', metavar='FILE',
                      help='Illumina manifest file (BPM or CSV)')
  parser.add_argument('-s', '--targetstrand', metavar='T', default='customer',
                      help='Target strand based on Illumina manifest file: ab, top, bottom, forward, '
                           'reverse, customer (default), anticustomer, design, antidesign')
  parser.add_argument('-t', '--gcthreshold', type=float, metavar='N', default=0,
                      help='Genotypes with GC score less than N set to missing')
  parser.add_argument('-w', '--warnings', action='store_true',
                      help='Emit warnings and A/B calls for SNPs with invalid manifest data')
  parser.add_argument('--samplestats', metavar='FILE',
                      help='Output per sample average GC statistics to FILE')
  parser.add_argument('--locusstats',  metavar='FILE',
                      help='Output per locus average GC statistics to FILE')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  genome = Genome()

  if options.loci is None:
    genome = Genome()
  elif is_str(options.loci):
    genome = load_genome(options.loci)

  if options.pedigree is None:
    phenome = Phenome()
  elif is_str(phenome):
    phenome = load_phenome(options.pedigree)

  abmap = {}

  errorhandler = None
  if options.warnings:
    def errorhandler(msg):
      sys.stderr.write('WARNING: %s\n' % msg)

  if options.manifest:
    sys.stderr.write('Processing Illumina manifest file...')
    abmap = create_Illumina_abmap(options.manifest,genome,targetstrand=options.targetstrand,errorhandler=errorhandler)
    sys.stderr.write('done.\n')

  if options.abmap:
    abmap.update(load_abmap(options.abmap))

  transform = GenoTransform.from_options(options)

  lbdfiles = iter(options.lbdfile)
  loci,gentrain,samples = read_Illumina_LBD(lbdfiles.next(),transform)

  loci = list(loci)
  for filename in lbdfiles:
    more_loci,more_gentrain,more_samples = read_Illumina_LBD(filename,transform)

    if list(more_loci) != loci:
      raise RuntimeError('Genotype headers do not match')

    samples = chain(samples,more_samples)

  if options.gcthreshold > 0:
    samples = filter_gc(samples, options.gcthreshold)

  if options.samplestats or options.locusstats:
    # Summary is both object and new stream
    summary = samples = GCSummary(loci,samples)

  genos = encode_genotypes(loci, samples, genome, phenome)

  if options.targetstrand != 'ab' and abmap:
    genos = recode_genotypes(genos, abmap)

  # Late removal of excluded loci and insurance on removal of samples.
  # renaming is supported, but not really a good idea in most cases
  genos = genos.transformed(transform)

  save_genostream(options.output,genos,format=options.outformat,genorepr=options.outgenorepr,hyphen=sys.stdout)

  if options.samplestats:
    out = table_writer(options.samplestats)
    out.writerow(['SAMPLE','GC_MEAN','GC_STDDEV'])
    out.writerows( [ (s,'%.4f' % gc, '%.4f' % dev)
                      for s,gc,dev in summary.samplestats ] )

  if options.locusstats:
    out = table_writer(options.locusstats)
    out.writerow(['LOCUS','GENTRAIN','GC_MEAN','GC_STDDEV'])
    out.writerows( [ (l,'%.4f' % gt, '%.4f' % gc, '%.4f' % dev)
                      for (l,gc,dev),gt in izip_exact(summary.locusstats,gentrain) ] )


if __name__ == '__main__':
  main()
