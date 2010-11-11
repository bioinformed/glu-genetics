# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert an Illumina Genotype Final Report (GFR) to a GLU raw genotype (GDAT) file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import csv

from   itertools                 import imap, islice, groupby, chain, izip
from   operator                  import itemgetter

import h5py
import numpy as np

from   glu.lib.fileutils         import autofile
from   glu.lib.illumina          import create_Illumina_abmap

from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.transform import GenoTransform


def progress_bar(samples, sample_count):
  try:
    from glu.lib.progressbar import progress_loop
  except ImportError:
    return samples

  return progress_loop(samples, length=sample_count, units='samples', update_interval=1)


def gfr_header(indata):
  num_snps = num_samples = manifest = None

  row = next(indata,None)

  if not row:
    return None,None,None,None

  elif row[0]=='[Header]':
    while 1:
      row = next(indata)

      if not row:
        continue
      elif row[0]=='Num SNPs':
        num_snps = int(row[1])
      elif row[0]=='Num Samples':
        num_samples = int(row[1])
      elif row[0]=='Content':
        manifest = row[1]
      elif row[0] == '[Data]':
        break

    header = next(indata)

    return header,indata,num_snps,num_samples,manifest

  elif 'SNP Name' in row and 'Sample ID' in row:
    header = row
    sample_idx = header.index('Sample ID')
    rows = [next(indata)]
    sample_id = rows[0][sample_idx]
    num_snps  = 1
    for row in indata:
      rows.append(row)
      if row[sample_idx] != sample_id:
        num_snps = len(rows)-1
        break

    return header,chain(rows,indata),num_snps,None,None


def read_gfr(filename):
  infile = autofile(filename)
  indata = csv.reader(infile,dialect='excel-tab')
  header,indata,num_snps,num_samples,manifest = gfr_header(indata)

  fields     = ['Sample ID','SNP Name','GC Score','Allele1 - Forward','Allele2 - Forward',
                'X','Y','X Raw','Y Raw','B Allele Freq','Log R Ratio']
  indices    = [ header.index(f) for f in fields ]

  def data():
    sample_idx,snp_idx,gc_idx,a1_idx,a2_idx,x_idx,y_idx,x_raw_idx,y_raw_idx,baf_idx,lrr_idx = indices

    alleles  = 'ACGTBID-'
    geno_map = dict( ((a1,a2), (a1+a2).replace('-',' ').replace('I','+').replace('D','-')) for a1 in alleles for a2 in alleles )

    samples_seen = set()
    for sample_id,data in groupby(indata,itemgetter(sample_idx)):
      data = list(data)

      assert len(data)==num_snps

      snp_names = map(itemgetter(snp_idx),data)
      if not samples_seen:
        snps = snp_names
      else:
        assert snps==snp_names

      assert sample_id not in samples_seen

      samples_seen.add(sample_id)

      geno  = imap(itemgetter(a1_idx,a2_idx),data)
      gc    = np.fromiter(imap(itemgetter(gc_idx),   data),    count=num_snps, dtype=np.float32)*10000
      x     = np.fromiter(imap(itemgetter(x_idx),    data),    count=num_snps, dtype=np.float32)*1000
      y     = np.fromiter(imap(itemgetter(y_idx),    data),    count=num_snps, dtype=np.float32)*1000
      x_raw = np.fromiter(imap(itemgetter(x_raw_idx),data),    count=num_snps, dtype=np.int16  )
      y_raw = np.fromiter(imap(itemgetter(y_raw_idx),data),    count=num_snps, dtype=np.int16  )
      baf   = np.fromiter( (d[baf_idx] or '0' for d in data),  count=num_snps, dtype=np.float32)*10000
      lrr   = np.fromiter( (d[lrr_idx] or '0' for d in data),  count=num_snps, dtype=np.float32)*10000

      yield sample_id,snps,geno,gc,x,y,x_raw,y_raw,baf,lrr

    infile.close()

  return num_snps,num_samples,manifest,data()


def create_gdat(filename, num_snps, num_samples=None):
  gdat         = h5py.File(filename, 'w')
  comp         = dict(compression='gzip',compression_opts=5)
  sample_chunk = 16
  if num_samples:
    sample_chunk = min(sample_chunk,num_samples)
  snp_chunk    = min(1024*16,num_snps)
  defchunks    = (sample_chunk,snp_chunk)
  defshape     = (None,num_snps)
  n            = num_samples or 1
  mdims        = (n, num_snps)
  shuffle      = True

  vstr         = h5py.special_dtype(vlen=str)
  snptype      = [ ('name',vstr),('chromosome',vstr),('location',np.uint32), ('alleles_forward','S2') ]

  gdat.create_dataset('SNPs',      (num_snps,), snptype,                 **comp)
  gdat.create_dataset('Samples',   (n,),        vstr, maxshape=(None,),  **comp)
  gdat.create_dataset('Genotype',  mdims, 'S2',
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('GC',        mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('X',         mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('Y',         mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('X_raw',     mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('Y_raw',     mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('BAF',       mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  gdat.create_dataset('LRR',       mdims, np.int16,
                                   maxshape=defshape,chunks=defchunks,shuffle=shuffle,**comp)
  return gdat


def gdat_writer(gdat, gfr_data, genome, abmap, transform):
  sample_chunk = 16

  gdat_SNPs     = gdat['SNPs']
  gdat_Samples  = gdat['Samples']
  gdat_Genotype = gdat['Genotype']
  gdat_GC       = gdat['GC']
  gdat_X        = gdat['X']
  gdat_Y        = gdat['Y']
  gdat_X_raw    = gdat['X_raw']
  gdat_Y_raw    = gdat['Y_raw']
  gdat_BAF      = gdat['BAF']
  gdat_LRR      = gdat['LRR']

  num_snps      = len(gdat['SNPs'])

  id_chunk      = []
  geno_chunk    = np.empty( (sample_chunk,num_snps), dtype='S2')
  gc_chunk      = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  x_chunk       = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  y_chunk       = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  x_raw_chunk   = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  y_raw_chunk   = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  baf_chunk     = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  lrr_chunk     = np.empty( (sample_chunk,num_snps), dtype=np.int16)

  include       = transform.samples.include
  exclude       = transform.samples.exclude
  rename        = transform.samples.rename

  # Resize storage to accommodate num_samples
  def resize(num_samples):
    d = (num_samples,num_snps)
    gdat_Samples.resize( (num_samples,) )
    gdat_Genotype.resize(d)
    gdat_GC.resize(d)
    gdat_X.resize(d)
    gdat_Y.resize(d)
    gdat_X_raw.resize(d)
    gdat_Y_raw.resize(d)
    gdat_BAF.resize(d)
    gdat_LRR.resize(d)

  # Save currently buffered samples
  def save_chunk(start):
    n    = len(id_chunk)
    stop = start + n

    if len(gdat_Samples)<stop:
      resize(stop)

    gdat_Samples[start:stop]  = id_chunk[:n]
    gdat_Genotype[start:stop] = geno_chunk[:n]
    gdat_GC[start:stop]       = gc_chunk[:n]
    gdat_X[start:stop]        = x_chunk[:n]
    gdat_Y[start:stop]        = y_chunk[:n]
    gdat_X_raw[start:stop]    = x_raw_chunk[:n]
    gdat_Y_raw[start:stop]    = y_raw_chunk[:n]
    gdat_BAF[start:stop]      = baf_chunk[:n]
    gdat_LRR[start:stop]      = lrr_chunk[:n]

  locusmap = genome.loci
  mapcache = {}

  saved = 0
  for sample_id,snp_names,geno,gc,x,y,x_raw,y_raw,baf,lrr in gfr_data:
    j = len(id_chunk) % sample_chunk

    if id_chunk and j==0:
      save_chunk(saved)
      saved   += sample_chunk
      id_chunk = []

    print 'Processing sample %d: %s' % (saved+j+1,sample_id),

    if include is not None and sample_id not in include:
      print ', skipped...'
      continue

    if exclude is not None and sample_id in exclude:
      print ', skipped...'
      continue

    if rename is not None:
      sample_id = rename.get(sample_id,sample_id)

    if not saved and not len(id_chunk):
      genomap  = []
      snp_recs = []

      for i,name in enumerate(snp_names):
        ab  = abmap[name]
        loc = locusmap[name]
        snp_recs.append( (name,loc.chromosome,loc.location,''.join(ab)) )

        gmap = mapcache.get(ab)
        if gmap is None:
          a,b = ab
          gmap = mapcache[ab] = { (a,a):'AA', (a,b):'AB', (b,a):'AB', (b,b):'BB', ('-','-'):'  ' }
        genomap.append(gmap)

      gdat_SNPs[:] = snp_recs

    id_chunk.append(sample_id)

    geno_chunk[j]  = [ gmap[g] for gmap,g in izip(genomap,geno) ]
    gc_chunk[j]    = gc
    x_chunk[j]     = x
    y_chunk[j]     = y
    x_raw_chunk[j] = x_raw
    y_raw_chunk[j] = y_raw
    baf_chunk[j]   = baf
    lrr_chunk[j]   = lrr

    print

  if id_chunk:
    save_chunk(saved)
    saved += len(id_chunk)

  # Ensure no extra space has been allocated
  resize(saved)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] manifest.bpm gfr.txt[.gz/.bz2] -o output.gdat'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--includesamples', dest='includesamples', metavar='FILE', action='append',
                    help='List of samples to include, all others will be skipped')
  parser.add_option('--excludesamples', dest='excludesamples', metavar='FILE', action='append',
                    help='List of samples to exclude, only samples not present will be kept')
  parser.add_option('--renamesamples', dest='renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output genotype file name')
  parser.add_option('-w', '--warnings', action='store_true', dest='warnings',
                    help='Emit warnings and A/B calls for SNPs with invalid manifest data')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=2 or not options.output:
    parser.print_help(sys.stderr)
    sys.exit(2)

  errorhandler = None
  if options.warnings:
    def errorhandler(msg):
      sys.stderr.write('WARNING: %s\n' % msg)

  sys.stderr.write('Loading Illumina manifest file...')
  genome = Genome()
  abmap  = create_Illumina_abmap(args[0],genome,targetstrand='forward',
                                         errorhandler=errorhandler)
  sys.stderr.write('done.\n')

  transform  = GenoTransform.from_object(options)

  num_snps,num_samples,manifest,gfr_data = read_gfr(args[1])
  gdat       = create_gdat(options.output, num_snps, num_samples)

  #gfr_data  = progress_bar(gfr_data,num_samples)
  gdat_writer(gdat, gfr_data, genome, abmap, transform)

  gdat.close()


if __name__=='__main__':
  main()
