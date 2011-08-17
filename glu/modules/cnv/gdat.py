import h5py
import sqlite3
import numpy as np
import numpy.lib.recfunctions as rfn

from   operator                  import itemgetter
from   itertools                 import groupby, izip, count
from   collections               import defaultdict

import scipy.stats

from   glu.lib.utils             import is_str,namedtuple
from   glu.lib.fileutils         import parse_augmented_filename,get_arg, \
                                        compressed_filename,namefile
from   glu.lib.glm               import Linear

IndexEntry = namedtuple('IndexEntry', 'id name gdat index manifest')


CHROM_LIST  = [ str(i) for i in range(1,23) ] + [ 'X', 'Y', 'XY', 'M' ]


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

    self.gdat     = h5py.File(filename,mode)
    self.filename = filename
    gdat          = self.gdat
    self.attrs    = gdat.attrs
    attrs         = gdat.attrs
    format_found  = attrs.get('GLU_FORMAT')
    gdat_version  = attrs.get('GLU_VERSION')
    snp_count     = attrs.get('SNPCount')
    sample_count  = attrs.get('SampleCount')

    if format_found!='gdat':
      raise ValueError('Input file "%s" does not appear to be in %s format.  Found %s.' \
                          % (namefile(filename),format,format_found))

    if gdat_version!=1:
      raise ValueError('Unknown gdat file version: %s' % gdat_version)

    if snp_count!=len(gdat['SNPs']):
      raise ValueError('Inconsistent gdat SNP metadata. gdat file may be corrupted.')

    if sample_count!=len(gdat['Genotype']):
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
        self.gdats[filename] = gdat = GDATFile(filename)

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


def print_regression_results(sample,scheme,linear_model):
  y     = linear_model.y
  b     = linear_model.beta.reshape(-1)
  stde  = (linear_model.ss*linear_model.W.diagonal())**0.5
  t     = b/stde
  n,m   = linear_model.X.shape
  p     = 2*scipy.stats.distributions.t.cdf(-abs(t),n-m)

  ss_t  = np.var(linear_model.y,ddof=1)
  r2    = 1 - linear_model.ss/ss_t

  print '  GC correct: ss_r=%.2f ss_t=%.2f r2=%.2f' % (linear_model.ss**0.5,ss_t**0.5,r2)
  #tvs   = ['%.2f' % tv for tv in t[1:]]
  #tvs   = []
  #out.writerow([sample,scheme]+tvs+[linear_model.ss**0.5,ss_t**0.5,r2])


def get_gcmodel(filename,chrom_indices,chrom_means=True):
  gcdata     = h5py.File(filename,'r')

  try:
    gc       = gcdata['GC'][:].T
    gcmeans  = gc.sum(axis=1)
    gc      -= gcmeans.reshape(-1,1)
    n        = len(gc)

    if not chrom_means:
      means  = np.ones((n,1), dtype=float)
    else:
      m      = len(CHROM_LIST)
      means  = np.zeros((n,m), dtype=float)
      for i,chrom in enumerate(CHROM_LIST):
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

  if thin is None:
    lm = Linear(lrr_masked, gcdesign_masked)
  else:
    lm = Linear(lrr_masked[::thin], gcdesign_masked[::thin])

  lm.fit()

  print_regression_results('','LRR', lm)

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
