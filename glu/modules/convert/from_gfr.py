import time

import sys
import csv

from   itertools         import imap, islice, groupby, chain, izip
from   operator          import itemgetter

import h5py
import numpy as np

from   glu.lib.fileutils import autofile
from   glu.lib.illumina  import IlluminaManifest

from   glu.lib.genolib.transform import GenoTransform


def progress_bar(samples, sample_count):
  try:
    from glu.lib.progressbar import progress_loop
  except ImportError:
    return samples

  update_interval = 1

  return progress_loop(samples, length=sample_count, units='samples', update_interval=update_interval)


def gfr_header(indata):
  num_snps = num_samples = None

  row = next(indata,None)

  if not row:
    return None,None,None

  elif row[0]=='[Header]':
    while 1:
      row = next(indata)

      if not row:
        continue
      elif row[0]=='Num SNPs':
        num_snps = int(row[1])
      elif row[0]=='Num Samples':
        num_samples = int(row[1])
      elif row[0] == '[Data]':
        break

    header = next(indata)

    return header,indata,num_snps,num_samples

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

    return header,chain(rows,indata),num_snps,None


def option_parser():
  import optparse

  usage = 'usage: %prog [options] gfr.txt[.gz/.bz2] -o output.gdat'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--includesamples', dest='includesamples', metavar='FILE', action='append',
                    help='List of samples to include, all others will be skipped')
  parser.add_option('--excludesamples', dest='excludesamples', metavar='FILE', action='append',
                    help='List of samples to exclude, only samples not present will be kept')
  parser.add_option('--renamesamples', dest='renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output genotype file name')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1 or not options.output:
    parser.print_help(sys.stderr)
    sys.exit(2)

  transform = GenoTransform.from_object(options)

  infile = autofile(args[0])
  indata= csv.reader(infile,dialect='excel-tab')
  header,indata,num_snps,num_samples = gfr_header(indata)

  out          = h5py.File(options.output)
  comp         = dict(compression='gzip',compression_opts=5)
  sample_chunk = 16
  snp_chunk    = 1024*16
  shuffle      = True
  n            = num_samples or 1

  vstr         = h5py.special_dtype(vlen=str)
  snps         = out.require_dataset('SNPs',      (num_snps,), vstr,                   **comp)
  samples      = out.require_dataset('Samples',   (n,),        vstr, maxshape=(None,), **comp)
  geno         = out.require_dataset('Genotype',  (n, num_snps), 'S2',
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  gc           = out.require_dataset('GC',        (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  x            = out.require_dataset('X',         (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  y            = out.require_dataset('Y',         (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  x_raw        = out.require_dataset('X_raw',     (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  y_raw        = out.require_dataset('Y_raw',     (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  baf          = out.require_dataset('BAF',       (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)
  lrr          = out.require_dataset('LRR',       (n, num_snps), np.int16,
                                     maxshape=(None,num_snps),chunks=(sample_chunk,snp_chunk),
                                     shuffle=shuffle,**comp)

  fields     = ['Sample ID','SNP Name','GC Score','Allele1 - Forward','Allele2 - Forward',
                'X','Y','X Raw','Y Raw','B Allele Freq','Log R Ratio']
  indices    = [ header.index(f) for f in fields ]
  sample_idx,snp_idx,gc_idx,a1_idx,a2_idx,x_idx,y_idx,x_raw_idx,y_raw_idx,baf_idx,lrr_idx = indices

  indata     = groupby(indata,itemgetter(sample_idx))
  #indata     = progress_bar(indata,num_samples)

  id_chunk     = []
  geno_chunk   = np.empty( (sample_chunk,num_snps), dtype='S2')
  gc_chunk     = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  x_chunk      = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  y_chunk      = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  x_raw_chunk  = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  y_raw_chunk  = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  baf_chunk    = np.empty( (sample_chunk,num_snps), dtype=np.int16)
  lrr_chunk    = np.empty( (sample_chunk,num_snps), dtype=np.int16)

  alleles = 'ACGTBID-'
  geno_map = dict( ((a1,a2), (a1+a2).replace('-',' ').replace('I','+').replace('D','-')) for a1 in alleles for a2 in alleles )

  samples_seen = set()
  include      = transform.samples.include
  exclude      = transform.samples.exclude
  rename       = transform.samples.rename

  # Resize storage to accommodate num_samples
  def resize(num_samples):
    d = (num_samples,num_snps)
    samples.resize( (num_samples,) )
    geno.resize(d)
    x.resize(d)
    gc.resize(d)
    y.resize(d)
    x_raw.resize(d)
    y_raw.resize(d)
    baf.resize(d)
    lrr.resize(d)

  # Save currently buffered samples
  def save_chunk(start):
    n    = len(id_chunk)
    stop = start + n

    if num_samples is None:
      resize(stop)

    samples[start:stop] = id_chunk[:n]
    gc[start:stop]      = gc_chunk[:n]
    geno[start:stop]    = geno_chunk[:n]
    x[start:stop]       = x_chunk[:n]
    y[start:stop]       = y_chunk[:n]
    x_raw[start:stop]   = x_raw_chunk[:n]
    y_raw[start:stop]   = y_raw_chunk[:n]
    baf[start:stop]     = baf_chunk[:n]
    lrr[start:stop]     = lrr_chunk[:n]

  saved = 0
  for sample_id,data in indata:
    j = len(id_chunk) % sample_chunk

    if id_chunk and j==0
      print 'Saving chunk',
      t = time.time()

      save_chunk(saved)
      saved   += sample_chunk
      id_chunk = []

      print ': %.2fs' % (time.time()-t)

    print 'Processing sample %d: %s' % (saved+j+1,sample_id),
    t = time.time()

    if include is not None and sample_id not in include:
      print ', skipped...'
      continue

    if exclude is not None and sample_id in exclude:
      print ', skipped...'
      continue

    if rename is not None:
      sample_id = rename.get(sample_id,sample_id)

    data = list(data)
    print ', materialize=%.2fs' % (time.time()-t),
    t1 = time.time()
    assert len(data)==num_snps

    snp_names = map(itemgetter(snp_idx),data)
    if not saved and not id_chunk:
      snps[:] = snps1 = snp_names
    else:
      assert snps1==snp_names

    assert sample_id not in samples_seen

    id_chunk.append(sample_id)
    samples_seen.add(sample_id)

    geno_chunk[j]  = [ geno_map[g] for g in imap(itemgetter(a1_idx,a2_idx),data) ]
    gc_chunk[j]    = np.fromiter(imap(itemgetter(gc_idx),   data), count=num_snps, dtype=np.float32)*10000
    x_chunk[j]     = np.fromiter(imap(itemgetter(x_idx),    data), count=num_snps, dtype=np.float32)*1000
    y_chunk[j]     = np.fromiter(imap(itemgetter(y_idx),    data), count=num_snps, dtype=np.float32)*1000
    x_raw_chunk[j] = np.fromiter(imap(itemgetter(x_raw_idx),data), count=num_snps, dtype=np.int16  )
    y_raw_chunk[j] = np.fromiter(imap(itemgetter(y_raw_idx),data), count=num_snps, dtype=np.int16  )
    baf_chunk[j]   = np.fromiter( (d[baf_idx] or '0' for d in data),  count=num_snps, dtype=np.float32)*10000
    lrr_chunk[j]   = np.fromiter( (d[lrr_idx] or '0' for d in data),  count=num_snps, dtype=np.float32)*10000
    print ', recode=%.2fs' % (time.time()-t1)

  if id_chunk:
    print 'Saving final chunk',
    t  = time.time()

    save_chunk(saved)
    saved += len(id_chunk)

    print ': %.2fs' % (time.time()-t)

  # Ensure no extra space has been allocated
  resize(saved)

  out.close()
  infile.close()


if __name__=='__main__':
  main()
