# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Convert an Illumina Locus by DNA report file into a GLU genotype file'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys
import csv

from   operator                  import itemgetter,getitem
from   itertools                 import islice,chain,imap,izip,groupby,repeat

from   numpy                     import array,zeros

from   glu.lib.utils             import izip_exact, is_str
from   glu.lib.fileutils         import autofile,table_reader,table_writer
from   glu.lib.sections          import read_sections
from   glu.lib.sequence          import norm_snp_seq,complement_base

from   glu.lib.genolib.locus     import Genome, Nothing, load_genome
from   glu.lib.genolib.phenos    import Phenome, load_phenome
from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.io        import save_genostream, geno_options, GenoTransform
from   glu.lib.genolib.genoarray import build_model, build_descr, GenotypeArray, GenotypeArrayDescriptor


def load_lbd_file(filename,options):
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

  data = csv.reader(autofile(filename), dialect='csv')

  skiprows(11)
  gentrain = parse_gentrain()
  skiprows(3)
  loci    = parse_loci()
  skiprows(5)
  samples = sample_generator()

  if options.samples.include is not None:
    samples = (s for s in samples if s[0]     in options.samples.include)

  if options.samples.exclude is not None:
    samples = (s for s in samples if s[0] not in options.samples.exclude)

  assert len(gentrain) == len(loci)

  return loci,gentrain,samples


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


def load_illumina_manifest(filename):
  '''
  Load an Illumina assay manifest file, parsing out the various sections
  '''
  ifile = csv.reader(autofile(filename),dialect='csv')
  sections = read_sections(ifile)

  while 1:
    heading,contents = sections.next()
    attrs = dict(c[:2] for c in islice(contents,10) if len(c) > 1)
    if 'Assay Format' in attrs:
      break

  format = attrs['Assay Format']
  # OPA manifest
  if format == 'Golden Gate':
    pass

  # Infinium
  # Known formats: Infinium,Infinium II,Infinium 2,Infinium HD Super
  elif format.startswith('Infinium'):
    heading,contents = sections.next()
    assert heading == 'Assay'
  else:
    raise ValueError('Unknown manifest format: %s' % format)

  return contents


def find_index(header,headings,optional=False):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass

  if not optional:
    raise ValueError('Cannot find heading index')


def parse_manifest(manifest,genome,abmap,targetstrand='customer',errorhandler=None):
  '''
  Parse an Illumina manifest to obtain mappings from A/B probes to alleles
  relative to a specific genomic strand
  '''
  if errorhandler is None:
    def errorhandler(msg):
      raise ValueError(msg)

  indelmap = {'D':'-','I':'+'}
  NA = ('N','A')
  targetstrand = targetstrand.lower()
  assert targetstrand in ('top','bottom','forward','reverse','customer','anticustomer','design','antidesign')

  manifest    = iter(manifest)
  header      = manifest.next()
  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  snp_idx     = find_index(header,['SNP'])
  chr_idx     = find_index(header,['Chr'])
  loc_idx     = find_index(header,['MapInfo'])
  cstrand_idx = find_index(header,['CustomerStrand','SourceStrand'])
  dstrand_idx = find_index(header,['IlmnStrand','Ilmn Strand'])
  topseq_idx  = find_index(header,['TopGenomicSeq'])
  probea_idx  = find_index(header,['AlleleA_ProbeSeq'],optional=True)
  probeb_idx  = find_index(header,['AlleleB_ProbeSeq'],optional=True)

  max_idx     = max(topseq_idx,dstrand_idx,cstrand_idx,snp_idx,name_idx,
                    assayid_idx,probea_idx,probeb_idx)

  for assay in manifest:
    assay  += ['']*(max_idx-len(assay)+1)
    locus   = assay[name_idx]
    snp     = assay[snp_idx]
    chr     = assay[chr_idx] or None
    loc     = assay[loc_idx]
    cstrand = assay[cstrand_idx].lower()
    dstrand = assay[dstrand_idx].lower()
    topseq  = assay[topseq_idx]
    gstrand = assay[assayid_idx].split('_')[-2]

    if chr=='Mt':
      chr='M'

    if cstrand not in cstrand in ('top','bot','p','m'):
      errorhandler('Invalid customer strand %s for %s' % (cstrand,locus))
      continue

    if dstrand not in ('top','bot','p','m'):
      errorhandler('Invalid design strand %s for %s' % (dstrand,locus))
      continue

    if gstrand not in 'FRU':
      errorhandler('Unknown gstrand %s for %s' % (gstrand,locus))
      continue

    if len(snp) != 5 or (snp[0],snp[2],snp[4]) != ('[','/',']'):
      errorhandler('Invalid SNP alleles %s for %s' % (snp,locus))
      continue

    a,b = snp[1],snp[3]

    if loc:
      loc = int(loc)

    # Handle CNV probes, which can simply be skipped
    if (a,b) == NA:
      genome.merge_locus(locus, chromosome=chr, location=loc)
      continue

    # Handle indels, which are strand neutral
    elif cstrand in 'pm' and dstrand in 'pm':
      a = indelmap[a]
      b = indelmap[b]
      genome.merge_locus(locus, chromosome=chr, location=loc)
      abmap[locus] = (a,b)
      continue

    # Otherwise, we have a SNP with A and B alleles on the design strand
    elif None not in (probea_idx,probeb_idx):
      probea = assay[probea_idx]
      probeb = assay[probeb_idx]
      if probea and probeb and (a,b) != (probea[-1],probeb[-1]):
        errorhandler('Design alleles do not match probes (%s,%s) != (%s,%s) for %s'
                         % (a,b,probea[-1],probeb[-1],locus))
        genome.merge_locus(locus, chromosome=chr, location=loc)
        continue

    try:
      tstrand,aa,bb = norm_snp_seq(topseq)
      tstrand       = tstrand.lower()
      if tstrand != 'top':
        errorhandler('Supplied sequence is not correctly normalize to top strand for %s' % locus)
        genome.merge_locus(locus, chromosome=chr, location=loc)
        continue

    except ValueError:
      tstrand,aa,bb = dstrand,a,b

    if dstrand!=tstrand:
      a,b = complement_base(a),complement_base(b)

    # a,b and aa,bb should both be normalized to tstrand and must be equal
    if (a,b) != (aa,bb):
      errorhandler('Assay alleles do not match sequence alleles (%s/%s != %s/%s) for %s' % (a,b,aa,bb,locus))
      genome.merge_locus(locus, chromosome=chr, location=loc)
      continue

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

      self.locusstats = [ (locus,mui,stdi) for locus,mui,stdi in izip(self.loci,mu,std) ]


def option_parser():
  import optparse

  usage = 'usage: %prog [options] lbdfile...'
  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,filter=True,transform=True,output=True)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output genotype file name')
  parser.add_option('-m', '--abmap', dest='abmap', metavar='FILE',
                    help='Mappings from A and B probes to other allele codings')
  parser.add_option('-M', '--manifest', dest='manifest', metavar='FILE',
                    help='Illumina manifest file')
  parser.add_option('-s', '--targetstrand', dest='targetstrand', metavar='T', default='customer',
                    help='Target strand based on Illumina manifest file: ab, top, bottom, forward, '
                         'reverse, customer (default), anticustomer, design, antidesign')
  parser.add_option('-t', '--gcthreshold', dest='gcthreshold', type='float', metavar='N', default=0,
                    help='Genotypes with GC score less than N set to missing')
  parser.add_option('-w', '--warnings', action='store_true', dest='warnings',
                    help='Emit warnings and A/B calls for SNPs with invalid manifest data')
  parser.add_option('--samplestats', dest='samplestats', metavar='FILE',
                    help='Output per sample average GC statistics to FILE')
  parser.add_option('--locusstats',  dest='locusstats',  metavar='FILE',
                    help='Output per locus average GC statistics to FILE')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help(sys.stderr)
    return

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
    manifest = load_illumina_manifest(options.manifest)
    parse_manifest(manifest,genome,abmap,targetstrand=options.targetstrand,errorhandler=errorhandler)
    sys.stderr.write('done.\n')

  if options.abmap:
    abmap.update(load_abmap(options.abmap))

  filter = GenoTransform.from_options(options)

  args = iter(args)
  loci,gentrain,samples = load_lbd_file(args.next(),filter)

  loci = list(loci)
  for arg in args:
    more_loci,more_gentrain,more_samples = load_lbd_file(arg,filter)

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
  genos = genos.transformed(transform=filter)

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
