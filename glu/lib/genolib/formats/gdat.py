# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GDAT file format read support'
__copyright__ = 'Copyright (c) 2010. BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

__all__ = ['load_gdat']

__genoformats__ = [
  #      LOADER                      SAVER                    WRITER            PFORMAT   ALIAS   EXTS
  ('load_gdat',                      None,                    None,             'sdat',   None,  'gdat') ]


from   itertools                 import izip
from   contextlib                import closing

import numpy as np
import h5py

from   glu.lib.utils             import is_str
from   glu.lib.fileutils         import parse_augmented_filename,get_arg,trybool,tryfloat, \
                                        compressed_filename,namefile,table_writer
from   glu.lib.genolib.locus     import Genome,Locus
from   glu.lib.genolib.phenos    import Phenome
from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,build_model
from   glu.lib.genolib.helpers   import encode_ab_to_snps


class GCSummary(object):
  def __init__(self, data, loci):
    self.loci = loci
    self.data = data

  def __iter__(self):
    self.locusstats  = []
    self.samplestats = samplestats = []

    locusstats1 = np.zeros(len(self.loci))
    locusstats2 = np.zeros(len(self.loci))
    locuscounts = np.zeros(len(self.loci),dtype=int)

    for sampleid,genos,gcscores in self.data:
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

    if not samplestats:
      return

    mu  = locusstats1/locuscounts

    # Compute std.dev from sqrt(E(X**2) - E(X)**2), with compensation for
    # the inherent numerical problems with the approach
    var = locusstats2/locuscounts - mu**2
    var[var < 0] = 0
    std = var**0.5

    self.locusstats = [ (locus,mui,stdi) for locus,mui,stdi in izip(self.loci,mu,std) ]


def load_models(gdat,ignoreloci=False):
  '''
  Load models from an HDF5 binary gdat file

  @param   gfile: gdat file
  @type    gfile: h5py HDF5 file instance
  '''
  print '!!! In load_models'

  model_cache = {}
  genome   = Genome()
  loci     = []
  models   = []

  for row in gdat['SNPs'][:].tolist():
    name,chromosome,location,alleles_forward = row[:4]

    model = model_cache.get(alleles_forward)

    if model is None:
      a,b     = alleles_forward
      genos   = [(a,a),(a,b),(b,b)] if a and b else [(a or b,a or b)]
      model_cache[alleles_forward] = model = build_model(genotypes=genos, max_alleles=2)

    loci.append(name)
    models.append(model)

    if ignoreloci:
      genome.loci[name] = Locus(name, model)
    else:
      if location == -1:
        strand   = None
        location = None
      else:
        strand   = '+'

      genome.loci[name] = Locus(name, model, chromosome, location, strand)

  return loci,genome,models


#######################################################################################


def load_gdat(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load the genotype matrix data from file.
  Note that the first row is header and the rest rows are genotypes,
  and the file is tab delimited.

  @param     filename: a file name or file object
  @type      filename: str or file object
  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param    chunksize: size of chunks to write/compress in bytes
  @type     chunksize: int
  @param      scratch: the buffer space available to use while reading or writing a binary file.
  @type       scratch: int
  @return:             format and sequence of column names followed by
                       tuples of row label and row data
  @rtype:              tuple of string and generator
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename     = parse_augmented_filename(filename,args)

  ignoreloci   =  trybool(get_arg(args,['ignoreloci']))
  gcthreshold  = tryfloat(get_arg(args,['gc','gcthreshold'])) or 0
  samplestats  =          get_arg(args,['samplestats'])
  locusstats   =          get_arg(args,['locusstats'])

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if not is_str(filename):
    raise ValueError('Invalid filename')

  if compressed_filename(filename):
    raise ValueError('Binary genotype files must not have a compressed extension')

  gdat         = h5py.File(filename,'r')
  attrs        = gdat.attrs
  format_found = attrs.get('GLU_FORMAT')
  gdat_version = attrs.get('GLU_VERSION')
  snp_count    = attrs.get('SNPCount')
  sample_count = attrs.get('SampleCount')

  if format_found!='gdat':
    raise ValueError('Input file "%s" does not appear to be in %s format.  Found %s.' \
                        % (namefile(filename),format,format_found))

  if gdat_version!=1:
    raise ValueError('Unknown gdat file version: %s' % gdat_version)

  if snp_count!=len(gdat['SNPs']):
    raise ValueError('Inconsistent gdat SNP metadata. gdat file may be corrupted.')

  if sample_count!=len(gdat['Genotype']):
    raise ValueError('Inconsistent gdat sample metadata. gdat file may be corrupted.')

  loci,file_genome,models = load_models(gdat,ignoreloci)

  phenome = Phenome()
  samples = gdat['Samples'][:].tolist()

  if gcthreshold<=0 and not samplestats and not locusstats:
    def _load():
      with closing(gdat):
        # Yield an initial dummy value to ensure that the generator starts,
        # so that gdat is closed properly when it shuts down
        yield

        descr = GenotypeArrayDescriptor(models)

        for name,genos in izip(samples,gdat['Genotype']):
          # Recode directly to binary as __:00,AA:01,AB:10,BB:11
          g = GenotypeArray(descr)
          g.data = encode_ab_to_snps(genos)

          yield name,g
  else:
    def _load():

      def genogc_iter():
        gc_scale = gdat['GC'].attrs['SCALE']
        gc_nan   = gdat['GC'].attrs['NAN']

        for name,genos,gc in izip(samples,gdat['Genotype'],gdat['GC']):
          mask     = gc==gc_nan
          gc       = gc.astype(float)
          gc      /= gc_scale
          gc[mask] = np.nan

          yield name,genos,gc

      with closing(gdat):
        # Yield an initial dummy value to ensure that the generator starts,
        # so that gdat is closed properly when it shuts down
        yield

        descr  = GenotypeArrayDescriptor(models)

        genogc = genogc_iter()

        if samplestats or locusstats:
          summary = genogc = GCSummary(genogc, loci)

        for name,genos,gc in genogc:
          if gcthreshold>0:
            mask        = gc<=gcthreshold
            mask       |= ~np.isfinite(gc)
            genos[mask] = '  '

          # Recode directly to binary as __:00,AA:01,AB:10,BB:11
          g = GenotypeArray(descr)
          g.data = encode_ab_to_snps(genos)

          yield name,g

      if samplestats:
        out = table_writer(samplestats)
        out.writerow(['SAMPLE','GC_MEAN','GC_STDDEV'])
        out.writerows( [ (s,'%.4f' % gc, '%.4f' % dev)
                          for s,gc,dev in summary.samplestats ] )

      if locusstats:
        out = table_writer(locusstats)
        out.writerow(['LOCUS','GC_MEAN','GC_STDDEV'])
        out.writerows( [ (l,'%.4f' % gc, '%.4f' % dev)
                          for l,gc,dev in summary.locusstats ] )

  # Create the loader and fire it up by requesting the first dummy element
  _loader = _load()
  _loader.next()

  genos = GenomatrixStream(_loader,'sdat',samples=samples,loci=loci,models=models,genome=file_genome,
                                          phenome=phenome,unique=True,packed=True)

  if genome:
    genos = genos.transformed(recode_models=genome)

  return genos


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
