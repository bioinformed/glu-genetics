import sys
import h5py
import sqlite3
import numpy as np
import numpy.lib.recfunctions as rfn

from   itertools                 import izip

from   glu.lib.utils             import is_str,namedtuple
from   glu.lib.fileutils         import parse_augmented_filename, compressed_filename, namefile
from   glu.lib.glm               import Linear


IndexEntry = namedtuple('IndexEntry', 'id name gdat index manifest')


CHROM_LIST  = [ str(i) for i in range(1,23) ] + [ 'X', 'Y', 'XY', 'M' ]


BAF_TYPE   = np.int16
BAF_SCALE  =  10000
BAF_NAN    = np.iinfo(BAF_TYPE).min
BAF_MIN    = np.iinfo(BAF_TYPE).min+1
BAF_MAX    = np.iinfo(BAF_TYPE).max

LRR_TYPE   = np.int32
LRR_SCALE  = 100000
LRR_NAN    = np.iinfo(LRR_TYPE).min
LRR_MIN    = np.iinfo(LRR_TYPE).min+1
LRR_MAX    = np.iinfo(LRR_TYPE).max


def lazy_property(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazy_property


def gdat_encode_f(data, scale, minval, maxval, nanval):
  data = data*scale
  np.clip(data, minval, maxval, out=data)
  data[~np.isfinite(data)] = nanval
  return data


def gdat_decode(data, scale, nanval):
  nans       = data==nanval
  data       = data.astype(float)
  data[nans] = np.nan
  data      /= scale
  return data


class GDATFile(object):
  def __init__(self,filename,mode='r',extra_args=None,**kwargs):
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    if not is_str(filename):
      raise ValueError('Invalid filename')

    if compressed_filename(filename):
      raise ValueError('Binary genotype files must not have a compressed extension')

    if mode not in ('r','r+'):
      raise ValueError('Unsupported GDAT file mode: %s' % mode)

    self.gdat          = h5py.File(filename,mode)
    self.filename      = filename
    gdat               = self.gdat
    self.attrs         = gdat.attrs
    self.format_found  = gdat.attrs.get('GLU_FORMAT')
    self.gdat_version  = gdat.attrs.get('GLU_VERSION')
    self.snp_count     = gdat.attrs.get('SNPCount')
    self.sample_count  = gdat.attrs.get('SampleCount')

    if self.format_found!='gdat':
      raise ValueError('Input file "%s" does not appear to be in %s format.  Found %s.' \
                          % (namefile(filename),format,self.format_found))

    if self.gdat_version!=1:
      raise ValueError('Unknown gdat file version: %s' % self.gdat_version)

    if self.snp_count!=len(gdat['SNPs']):
      raise ValueError('Inconsistent gdat SNP metadata. gdat file may be corrupted.')

    if self.sample_count!=len(gdat['Genotype']):
      raise ValueError('Inconsistent gdat sample metadata. gdat file may be corrupted.')

  @lazy_property
  def snps(self):
    return self.gdat['SNPs'][:]

  @lazy_property
  def samples(self):
    return self.gdat['Samples'][:]

  @lazy_property
  def sample_index(self):
    return dict( (name,i) for i,name in enumerate(self.samples) )

  @lazy_property
  def chromosome_index(self):
    snps     = self.snps
    n        = len(snps)
    index    = np.arange(n).astype([('index',int)])
    chroms   = snps['chromosome'].astype('S10')
    chromset = set(chroms)
    chroms   = chroms.astype([('chromosome','S10')])
    locs     = snps['location'].astype(int).astype([('location',int)])
    snps     = rfn.merge_arrays([chroms,locs,index])

    snps.sort()

    index = {}
    for name in chromset:
      mask        = snps['chromosome']==name
      indices     = snps['index'][mask]
      pos         = snps['location'][mask]

      if name.startswith('chr'):
        name = name[3:]
      if name.upper()=='MT':
        name = 'M'

      index[name] = pos,indices

    return index

  def cnv_data(self,index):
    if isinstance(index,str):
      index = self.sample_index[index]

    if not isinstance(index,int):
      raise ValueError('Unsupported addressing mode: %s' % index)

    gdat       = self.gdat
    genos      = gdat['Genotype']

    if 'LRR_QN' in gdat:
      lrr      = gdat['LRR_QN']
      baf      = gdat['BAF_QN']
    else:
      lrr      = gdat['LRR']
      baf      = gdat['BAF']

    lrr_scale  = lrr.attrs['SCALE']
    lrr_nan    = lrr.attrs['NAN']

    baf_scale  = baf.attrs['SCALE']
    baf_nan    = baf.attrs['NAN']

    sample     = self.samples[index]
    sample_lrr = gdat_decode(lrr[index], lrr_scale, lrr_nan)
    sample_baf = gdat_decode(baf[index], baf_scale, baf_nan)

    return sample,genos[index],sample_lrr,sample_baf

  def cnv_iter(self,transform=None):
    gdat       = self.gdat
    snps       = self.snps
    samples    = self.samples

    genos      = gdat['Genotype']

    if 'LRR_QN' in gdat:
      lrr      = gdat['LRR_QN']
      baf      = gdat['BAF_QN']
    else:
      lrr      = gdat['LRR']
      baf      = gdat['BAF']

    lrr_scale  = lrr.attrs['SCALE']
    lrr_nan    = lrr.attrs['NAN']

    baf_scale  = baf.attrs['SCALE']
    baf_nan    = baf.attrs['NAN']

    if transform is not None:
      include  = transform.samples.include
      exclude  = transform.samples.exclude
      rename   = transform.samples.rename
    else:
      include  = exclude = rename = None

    for i,sample in enumerate(samples):
      if include is not None and sample not in include:
        continue

      if exclude is not None and sample in exclude:
        continue

      if rename is not None:
        sample = rename.get(sample,sample)

      if not sample:
        sys.stderr.write('Invalid null sample name... skipping.\n')
        continue

      sample_lrr = gdat_decode(lrr[i], lrr_scale, lrr_nan)
      sample_baf = gdat_decode(baf[i], baf_scale, baf_nan)

      yield sample,snps,genos[i],sample_lrr,sample_baf

  def __getitem__(self,item):
    return self.gdat[item]

  def __len__(self):
    return self.sample_count

  def close(self):
    self.gdat.close()


def gdat_decoder(table,rows):
  if 'SCALE' not in table.attrs:
    return rows

  def _decoder(table,rows):
    scale = table.attrs['SCALE']
    nan   = table.attrs['NAN']
    for row in rows:
      yield gdat_decode(row, scale, nan)

  return _decoder(table,rows)


def block_iter(table):
  chunksize = table.chunks[0]
  start     = 0
  last      = len(table)

  while start<last:
    stop  = min(last,start+chunksize)
    chunk = table[start:stop]
    for row in chunk:
      yield row
    start = stop


def parallel_gdat_iter(*tables):
  iters = [ gdat_decoder(table,block_iter(table)) for table in tables ]
  return izip(*iters)


class BatchTableWriter(object):
  def __init__(self,table):
    self.table     = table
    self.batchsize = table.chunks[0]
    self.scale     = table.attrs['SCALE']
    self.nan       = table.attrs['NAN']
    self.min       = table.attrs['MIN']
    self.max       = table.attrs['MAX']
    self.batch     = []
    self.index     = 0

  def write(self, value):
    batch = self.batch
    batch.append(value.reshape(-1))

    if len(batch)>=self.batchsize:
      self.flush()

  def flush(self):
    batch = self.batch

    if not batch:
      return

    table = self.table

    if len(batch)==1:
      chunk             = gdat_encode_f(batch[0], self.scale, self.min, self.max, self.nan)
      table[self.index] = chunk
      self.index       += 1
    else:
      chunk             = [ gdat_encode_f(c, self.scale, self.min, self.max, self.nan) for c in batch ]
      chunk             = np.array(chunk, dtype=table.dtype)

      start             = self.index
      end               = start+len(batch)

      table[start:end]  = chunk

      self.index        = end

    batch[:] = []

  def close(self):
    self.flush()
    self.table = None


def create_gdat_qn(gdatobject,s,n):
  comp         = dict(compression='gzip',compression_opts=5)
  chunks       = (1,s)
  shape        = (n,s)
  shuffle      = False

  gdat         = gdatobject.gdat
  BAF_QN       = gdat.require_dataset('BAF_QN', shape, BAF_TYPE,
                                      maxshape=shape,chunks=chunks,shuffle=shuffle,
                                      fillvalue=BAF_NAN,**comp)
  BAF_QN.attrs['SCALE'] = BAF_SCALE
  BAF_QN.attrs['NAN']   = BAF_NAN
  BAF_QN.attrs['MIN']   = BAF_MIN
  BAF_QN.attrs['MAX']   = BAF_MAX

  LRR_QN       = gdat.require_dataset('LRR_QN', shape, LRR_TYPE,
                                      maxshape=shape,chunks=chunks,shuffle=shuffle,
                                      fillvalue=LRR_NAN,**comp)

  LRR_QN.attrs['SCALE'] = LRR_SCALE
  LRR_QN.attrs['NAN']   = LRR_NAN
  LRR_QN.attrs['MIN']   = LRR_MIN
  LRR_QN.attrs['MAX']   = LRR_MAX


def sqlite_magic(con):
  con.execute('PRAGMA synchronous=OFF;')
  con.execute('PRAGMA journal_mode=OFF;')
  con.execute('PRAGMA count_changes=OFF;')
  con.execute('PRAGMA cache_size=200000;')
  con.execute('PRAGMA default_cache_size=200000;')


class GDATIndex(object):
  def __init__(self,filename):
    self.filename = filename
    self.con      = sqlite3.connect(filename)
    self.gdats    = {}

    sqlite_magic(self.con)

    try:
      sql = '''
      CREATE TABLE gdatindex (id            INTEGER PRIMARY KEY,
                              sample        TEXT,
                              gdat          TEXT,
                              idx           INTEGER,
                              manifest      TEXT);'''

      cur = self.con.cursor()
      cur.execute(sql)
    except sqlite3.OperationalError:
      pass

    try:
      sql = 'CREATE UNIQUE INDEX idx_gdatindex_id ON gdatindex (id)';
      cur.execute(sql)
    except sqlite3.OperationalError:
      pass

    try:
      sql = 'CREATE INDEX idx_gdatindex_sample ON gdatindex (sample)';
      cur.execute(sql)
    except sqlite3.OperationalError:
      pass


  def query(self,name):
    sql = '''SELECT *
             FROM   gdatindex
             WHERE  sample=?'''

    make = IndexEntry._make

    cur = self.con.cursor()
    cur.execute(sql, (name,))
    return [ make(r) for r in cur ]


  def get(self,names):
    results = self.query(names)

    for i,sample,filename,index,manifest in results:
      gdat = self.gdats.get(filename)

      if gdat is None:
        try:
          self.gdats[filename] = gdat = GDATFile(filename)
        except:
          sys.stderr.write('[ERROR] Cannot open GDAT file: %s\n' % filename)
          raise

      yield gdat,index


  def clear_index(self,gdat):
    sql = '''DELETE
             FROM   gdatindex
             WHERE  gdat = ?'''

    cur = self.con.cursor()
    cur.execute(sql, (gdat.filename,))


  def index(self,gdat):
    self.clear_index(gdat)

    filename = gdat.filename
    manifest = gdat.attrs['ManifestName']
    rows = [ (name,filename,i,manifest) for i,name in enumerate(gdat.samples) ]

    sql = 'INSERT INTO gdatindex VALUES (NULL,?,?,?,?);'
    cur = self.con.cursor()
    cur.executemany(sql, rows)
    self.con.commit()


def get_gcmodel(filename,chrom_indices=None,extra_terms=0):
  gcdata     = h5py.File(filename,'r')

  try:
    gc       = gcdata['GC'][:].T
    gcmeans  = gc.sum(axis=1)
    gc      -= gcmeans.reshape(-1,1)
    n        = len(gc)

    if chrom_indices is None:
      means  = np.ones((n,1+extra_terms), dtype=float)
    else:
      m      = len(CHROM_LIST)+extra_terms
      means  = np.zeros((n,m), dtype=float)
      for i,chrom in enumerate(CHROM_LIST):
        if chrom in chrom_indices:
          pos,index      = chrom_indices[chrom]
          means[index,i] = 1

    gcdesign = np.hstack( [means,gc] )
    gcmask   = np.isfinite(gcdesign.sum(axis=1))

    return gcdesign,gcmask

  finally:
    gcdata.close()


def gc_correct(lrr,gcdesign,gcmask,minval=None,maxval=None,thin=None):
  mask  = np.isfinite(lrr)
  mask &= gcmask

  if minval is not None:
    mask &= lrr>=minval

  if maxval is not None:
    mask &= lrr<=maxval

  lrr_masked      = lrr[mask]
  lrr_masked     -= lrr_masked.mean()
  gcdesign_masked = gcdesign[mask]

  if not len(lrr_masked):
    return lrr

  if thin is None:
    lm = Linear(lrr_masked, gcdesign_masked)
  else:
    lm = Linear(lrr_masked[::thin], gcdesign_masked[::thin])

  lm.fit()

  print '  GC correct: r2=%.2f, p=%s' % (lm.r2(),lm.p_values(phred=True))

  beta = lm.beta.reshape(-1,1)

  n         = len(CHROM_LIST)
  weights   = gcdesign[:,:n].sum(axis=0)
  weights  /= weights.sum()
  alphas    = (beta[:n].reshape(-1)*weights)
  alpha     = alphas.sum()

  if 0:
    print '!!!   alpha=',alpha
    print '      betas=',beta[:n].reshape(-1)
    print '    weights=',weights
    print '     alphas=',alphas

  beta[:n]  = 0
  lrr_adj   = lrr-np.dot(gcdesign,beta).reshape(-1)
  lrr_adj  -= alpha

  return lrr_adj
