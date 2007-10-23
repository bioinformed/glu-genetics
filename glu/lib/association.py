# -*- coding: utf-8 -*-
'''
File:          association.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:      Support library for association testing

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   math              import ceil,log
from   itertools         import izip
from   collections       import defaultdict
from   operator          import itemgetter

from   numpy             import array,matrix,asarray,asanyarray,zeros,zeros_like, \
                                exp,nan,abs,arange,median,sqrt,inf
from   scipy             import stats

from   glu.lib.utils     import tally,pick
from   glu.lib.fileutils import autofile,namefile,load_list,load_map,load_table
from   glu.lib.genolib   import load_genostream, get_genorepr


LOGE_10 = log(10)

def log10(n):
  if n in (0.0,-0.0,0):
    return -inf
  return log(n)/LOGE_10


def format_pvalue(p,scale=8):
  '''
  Format pvalues sensibly

  >>> format_pvalue(1)
  '1.0000000'
  >>> format_pvalue(1./3)
  '0.3333333'
  >>> format_pvalue(0.001256)
  '0.0012560'
  >>> format_pvalue(0.0001256)
  '1.256e-04'
  >>> format_pvalue(0.00001256)
  '1.256e-05'
  '''
  if -2*log10(p)>(scale-1):
    return '%%#0.0%de' % (scale-5) % p
  else:
    return '%%#0.0%df' % (scale-1) % p


def contingency_table_indicators(ys,gs):
  '''
  Build contingency table from categorical data

    ys = y0,y1,y2,..,yn: vector of phenotype indicators
    gs = g0,g1,g2,..,gm: vector of genotype indicators

  >>> ys = [(1,1,1,0,0,0,0,0,0),
  ...       (0,0,0,1,1,1,0,0,0),
  ...       (0,0,0,0,0,0,1,1,1)]
  >>> gs = [(1,1,0,0,0,0,0,1,1),
  ...       (0,0,1,0,1,0,0,0,0),
  ...       (0,0,0,1,0,1,1,0,0)]
  >>> contingency_table_indicators(ys,gs)
  array([[2, 1, 0],
         [0, 1, 2],
         [2, 0, 1]])
  '''
  ys = asanyarray(ys,dtype=int)
  gs = asanyarray(gs,dtype=int)

  assert ys.shape == gs.shape

  c = zeros( (len(ys),len(gs)), dtype=int )
  for i,y in enumerate(ys):
    for j,g in enumerate(gs):
      c[i,j] = (y & g).sum()
  return c


def contingency_table_categorical(ys,gs):
  '''
  Build contingency table from categorical data

    ys: vector of phenotype categories (0..n)
    gs: vector of genotype categories (0..m)

  >>> ys = [0,0,0,1,1,1,2,2,2]
  >>> gs = [0,0,1,2,1,2,2,0,0]
  >>> contingency_table_categorical(ys,gs)
  array([[2, 1, 0],
         [0, 1, 2],
         [2, 0, 1]])
  '''
  ys = asanyarray(ys,dtype=int)
  gs = asanyarray(gs,dtype=int)

  assert ys.shape == gs.shape

  n,m = max(ys)+1,max(gs)+1
  c = zeros( (n,m), dtype=int )

  for y,g in izip(ys,gs):
    c[y,g] += 1

  return c


contingency_table = contingency_table_categorical


def contingency_analysis(c,mincell=5,minmargin=15):
  if not c.size:
    return nan,0,array([],dtype=float)

  if minmargin>0:
    # Remove columns with less that the required margins
    c = c[:,c.sum(axis=0)>=minmargin]

    # Remove rows with less that the required margins
    if c.size:
      c = c[c.sum(axis=1)>=minmargin,:]

  # Remove empty columns
  if c.size:
    c = c[:,c.min(axis=0)>0]

  # Remove empty rows
  if c.size:
    c = c[c.min(axis=1)>0,:]

  m,n = c.shape
  df = (m-1)*(n-1)

  if df <= 0:
    return nan,0,array([],dtype=float)

  if c.min() <= mincell:
    import rpy
    cc = rpy.array(c.tolist())
    p = rpy.r.fisher_test(cc, conf_int=False, workspace=500000)['p.value']
    t = stats.distributions.chi2.ppf(1-p,df)

  else:
    m1 = c.sum(axis=0)
    m2 = c.sum(axis=1)
    e  = matrix(m2).T*matrix(m1)/float(sum(m1))
    t  = float(((c-e).A**2/e).sum())

  ors = zeros( (m-1,n-1), dtype=float )

  for i in range(m-1):
    for j in range(n-1):
      ors[i,j] = (float(c[0  ,0] * c[i+1,j+1])
               /       (c[i+1,0] * c[  0,j+1]))

  return t,df,ors


# Cochrane-Armitage Trend Test for 2xc tables
def table_trend(x):
  '''
  >>> x = [[26,26,23,18,9],[6,7,9,14,23]]
  >>> table_trend(x)
  22.9614677639
  '''
  x = asanyarray(x, dtype=int)

  if len(x.shape) != 2 or x.shape[0] != 2:
    raise ValueError,'tabletrend requires a 2xc table'

  n_i   = x.sum(axis=0)
  n     = n_i.sum()
  r_i   = arange(len(n_i))
  r_bar = float((n_i*r_i).sum())/n
  s2    = (n_i*(r_i-r_bar)**2).sum()
  p1    = float(x[0].sum())/n

  t = (x[0]*(r_i-r_bar)/(p1*(1-p1)*s2)**0.5).sum()**2

  return t


def permutations_needed(p_hat, w, g):
  '''
  Estimate the number of permutations required to estimate an empirical
  p-value of a test statistic using a Monte Carlo permutation procedure with N
  replicate permutations.  We choose N such that the estimated empirical
  p-value, p_hat, is within a proportion w (the width parameter) of its true
  p-value, p, with predetermined confidence probability g. That is, we want
  the standard deviation sigma_p_hat of p_hat to be proportional to p_hat.
  This permutation process can be viewed as a set of N independent Bernoulli
  trials each with success probability p.

  Using a Normal approximation for the distribution of p_hat, we obtain

       1 - p_hat
  N = -----------  Phi^-2[ (g+1)/2 ]
      w^2 * p_hat

  where Phi^-2 is the squared inverse of the standard normal cumulative
  distribution function.

  For example, to estimate an empirical p-value that is within 20% of its
  true value with 95% confidence (w=.2, g=.95), then N is approximately
  N ~= 100*(1-p_hat)/p_hat.  For suitably small values of p_hat, this
  further simplifies to N ~= 100/p_hat.

  >>> permutations_needed(0.05, 0.5, 0.90)
  206
  >>> permutations_needed(0.05, 0.2, 0.95)
  1825
  >>> permutations_needed(0.01, 0.5, 0.90)
  1072
  >>> permutations_needed(0.01, 0.2, 0.95)
  9508
  >>> permutations_needed(0.001, 0.2, 0.95)
  95941
  >>> permutations_needed(0.0001, 0.2, 0.95)
  960269
  >>> permutations_needed(0.00001, 0.2, 0.95)
  9603552
  '''
  return int(ceil((1-p_hat)/(w*w*p_hat) * stats.distributions.norm.ppf( (g+1)/2. )**2))


def normalize(row,n):
  '''
  Normalize a list to a fixed size by cutting off extra elements and adding
  empty strings, as necessary such that len(normalize(row,n)) == n.

  >>> normalize([],2)
  ['', '']
  >>> normalize(['  '],2)
  ['', '']
  >>> normalize(['1','2','3'],2)
  ['1', '2']
  '''
  return map(str.strip,row[:n]) + ['']*(n-len(row))


def strip_trailing_empty(row):
  '''
  Strip trailing empty elements in a list.

  >>> strip_trailing_empty([1,2,3,'','',''])
  [1, 2, 3]
  >>> strip_trailing_empty(['','',''])
  []
  >>> strip_trailing_empty(['','','1'])
  ['', '', '1']
  '''
  i = len(row)-1
  while i>=0 and not row[i]:
    i -= 1
  return row[:i+1]


# FIXME: Needs docs+tests
def load_phenos(filename,deptype=int,allowdups=False,verbose=1,errs=sys.stderr):
  '''
  Load phenotypes from a tab-delimited file

  Format:
    Row  1  : Column headers
    Rows 2..: Phenotype records

    Column  1  : Subject ID that cooresponds to genotype subject ID
    Column  2  : Phenotype, either categorical or continious
    Columns 3..: Additional integer or real valued covariates

  NOTE: Subjects may be repeated to allow for incidence-density and other
        schemes where subject genotypes may be included in a model multiple
        times with potentially different phenotypes.
  '''
  phenos = load_table(filename,want_header=True)

  try:
    header = strip_trailing_empty(phenos.next())
  except StopIteration:
    raise ValueError('Empty phenotype file')

  if len(header) < 2:
    header = normalize(header,2)

  for i,h in enumerate(header):
    if not h:
      if i == 0:
        header[i] = 'ID'
      elif i == 1:
        header[i] = 'PHENO'
      else:
        header[i] = 'Covariate_%02d' % (i-1)

  dups = [ (h,n) for h,n in tally(header).iteritems() if n>1 ]
  if dups:
    errs.write('[ERROR] Duplicate headers detected (n=%d)' % len(dups))
    if verbose > 1:
      for h,n in sorted(dups, key=itemgetter(1,0), reverse=True):
        errs.write('         %-12s : %2d\n' % (h,n))
    raise ValueError('Invalid repeated heading in phenotype header')

  def _phenos():
    warn_msg = '[WARNING] Subject "%s" dropped on line #%04d due to missing data in column "%s"\n'
    note_msg = '[NOTICE] Read phenotypes from %s: ' \
               '%d covariates, %d valid records, %d records dropped, %d distinct subjects\n'

    dropped  = 0
    subjects = defaultdict(int)
    for i,row in enumerate(phenos):
      row = normalize(strip_trailing_empty(row),len(header))

      if '' in row:
        dropped += 1
        if verbose > 1:
          j = row.index('')
          errs.write(warn_msg % (row[0],i+2,header[j]))
        continue

      subjects[row[0]] += 1
      row[1]  = deptype(row[1])
      row[2:] = map(float, row[2:])

      yield row

    if verbose:
      errs.write(note_msg % (namefile(filename),len(header)-2,i+1-dropped,dropped,len(subjects)))

    # Check for duplicate subjects
    dups = [ (pid,n) for pid,n in subjects.iteritems() if n>1 ]
    if dups:
      if allowdups:
        errs.write('[NOTICE] Duplicate subjects (n=%d)\n' % len(dups))
      else:
        errs.write('[ERROR] Duplicate subjects (n=%d)\n' % len(dups))

      if verbose > 1:
          for pid,n in sorted(dups, key=itemgetter(1,0), reverse=True):
            errs.write('         %-12s : %2d\n' % (pid,n))

      if not allowdups:
        raise ValueError('Duplicate subjects detected')

  return header,_phenos()


def build_models(phenofile, genofile, options,deptype=int):
  header,phenos = load_phenos(phenofile,deptype=deptype,allowdups=options.allowdups)
  phenos        = list(phenos)
  phenocount1   = len(phenos)
  genorepr      = get_genorepr(options.genorepr)
  loci          = load_genostream(genofile,format=options.format,genorepr=genorepr).as_ldat()

  if options.includeloci or options.excludeloci:
    loci = loci.transformed(include_loci=options.includeloci,
                            exclude_loci=options.excludeloci)

  if options.includesamples:
    keep   = set(load_list(options.includesamples))
    phenos = (p for p in phenos if p[0] in keep)

  if options.excludesamples:
    drop   = set(load_list(options.excludesamples))
    phenos = (p for p in phenos if p[0] not in drop)

  if loci.samples:
    samples = set(loci.samples)
    phenos  = (p for p in phenos if p[0] in samples)

  phenos      = list(phenos)
  phenocount2 = len(phenos)

  if phenocount1 != phenocount2:
    print >> sys.stderr, '[NOTICE] After exclusions, %d subjects remain, %d subjects excluded' % (phenocount2,phenocount1-phenocount2)

  reference_alleles = load_map(options.refalleles) if options.refalleles else None

  models = LocusModelBuilder(loci.samples,header,phenos,
                             reference_alleles=reference_alleles,
                             minmaf=options.minmaf,mingenos=options.mingenos)

  return loci,models


def estimate_maf(genocounts):
  '''
  Estimate minor allele frequency (MAF) from a dictionary of
  genotypes->counts.

  Genotypes must be two-item sequences of alleles, except the missing
  genotype which must be the only value to evaluate to False.

  >>> estimate_maf( {'AA':1,    'AB':2,    'BB':1   } )
  0.5
  >>> estimate_maf( {'AA':1000, 'AB':2000, 'BB':1000} )
  0.5
  >>> estimate_maf( {'AA':1000, 'AB':2000, 'BB':1000, None:10000 } )
  0.5
  >>> estimate_maf( {('A','A'):2000, ('A','B'):2000, None:10000 } )
  0.25
  >>> estimate_maf( {('A','A'):1000} )
  0.0
  >>> estimate_maf( {None:1000} )
  0.0
  '''
  hom1 = hom2 = het = 0

  for g,n in genocounts.iteritems():
    if not g:
      continue
    elif g[0] != g[1]:
      if het:
        raise ValueError('Locus may only have two alleles')
      het = n
    elif hom1:
      hom2 = n
    elif not hom1:
      hom1 = n
    else:
      raise ValueError('Locus may only have two alleles')

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)

  n = (hom1+het+hom2)*2

  if not n:
    return 0.0

  return (2.0*hom1 + het)/n


###################################################################


class TERM(object):
  def __init__(self, lname):
    self.lname = lname
    self.index = None

  def loci(self):
    if self.lname is None:
      return []
    return [self.lname]

  def terms(self):
    return [self]

  def __mul__(self, other):
    t1 =  self.terms if isinstance(self,  INTERACTION) else [self]
    t2 = other.terms if isinstance(other, INTERACTION) else [other]
    return INTERACTION(t1 + t2)

  def __add__(self, other):
    t1 =  self.terms if isinstance(self,  COMBINATION) else [self]
    t2 = other.terms if isinstance(other, COMBINATION) else [other]
    return COMBINATION(t1 + t2)


class NULL(TERM):
  def __init__(self, lname=None):
    self.lname = lname

  def effects(self, loci, i):
    if self.lname is None:
      return []

    lmodel = loci[self.lname]
    if lmodel.genos[i] not in lmodel.genomap:
      return [None]
    return []

  def names(self, loci):
    return []

  def __len__(self):
    return 0


class GENO(TERM):
  def effects(self, loci, i):
    lmodel  = loci[self.lname]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]

    if geno not in genomap:
      return [None]

    genos = [0,0]
    j = genomap[geno]
    if j:
      genos[j-1] = 1

    return genos

  def names(self, loci):
    lmodel = loci[self.lname]
    if lmodel.genocount < 3:
      return []
    return [ '%s:%s' % (self.lname,''.join(g)) for g in lmodel.tests[1:] ]

  def estimates(self,p):
    return p[self.index:self.index+2,0].A.flatten()

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i],c[i+1,i+1]]

  def se(self,c):
    return sqrt(self.var(c))

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1,p2 = self.estimates(p)
    e1,e2 = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1)),
            (exp(p2-a*e2),exp(p2+a*e2))]

  def __len__(self):
    return 2


class ADOM(TERM):
  def effects(self, loci, i):
    lmodel  = loci[self.lname]
    genomap = loci[self.lname].genomap
    geno    = lmodel.genos[i]
    if geno not in lmodel.genomap:
      return [None]
    return ([0,0],[1,1],[2,0])[lmodel.genomap[geno]]

  def names(self, loci):
    lmodel = loci[self.lname]
    if lmodel.genocount < 3:
      return []
    return [ '%s:trend:-%s+%s' % (self.lname,lmodel.alleles[0],lmodel.alleles[1]),
             '%s:domdev:%s%s'  % (self.lname,lmodel.alleles[0],lmodel.alleles[1]) ]

  def estimates(self,p):
    return p[self.index:self.index+2,0].A.flatten()

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i],c[i+1,i+1]]

  def se(self,c):
    return sqrt(self.var(c))

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1,p2 = self.estimates(p)
    e1,e2 = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1)),
            (exp(p2-a*e2),exp(p2+a*e2))]

  def __len__(self):
    return 2


class TREND(TERM):
  def effects(self, loci, i):
    lmodel = loci[self.lname]
    return [ lmodel.genomap.get(lmodel.genos[i]) ]

  def names(self, loci):
    lmodel = loci[self.lname]
    if lmodel.genocount < 2:
      return []
    return ['%s:trend:-%s+%s' % (self.lname,lmodel.alleles[0],lmodel.alleles[1])]

  def estimates(self,p):
    return [p[self.index,0],p[self.index,0]*2]

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i],4*c[i,i]]

  def se(self,c):
    return sqrt(self.var(c))

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1,p2 = self.estimates(p)
    e1,e2 = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1)),
            (exp(p2-a*e2),exp(p2+a*e2))]

  def __len__(self):
    return 1


class DOM(TERM):
  def effects(self, loci, i):
    lmodel  = loci[self.lname]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]
    if geno not in genomap:
      return [None]
    return [ int(genomap.get(geno)>0) ]

  def names(self, loci):
    lmodel = loci[self.lname]
    if lmodel.genocount < 2:
      return []
    return ['%s:dom:%s<%s' % (self.lname,lmodel.alleles[0],lmodel.alleles[1]) ]

  def estimates(self,p):
    return [p[self.index,0],p[self.index,0]]

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i],c[i,i]]

  def se(self,c):
    return sqrt(self.var(c))

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1,p2 = self.estimates(p)
    e1,e2 = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1)),
            (exp(p2-a*e2),exp(p2+a*e2))]

  def __len__(self):
    return 1


class REC(TERM):
  def effects(self, loci, i):
    lmodel  = loci[self.lname]
    genomap = lmodel.genomap
    geno    = lmodel.genos[i]
    if geno not in genomap:
      return [None]
    return [ int(genomap.get(geno)>1) ]

  def names(self, loci):
    lmodel = loci[self.lname]
    if lmodel.genocount < 2:
      return []
    return ['%s:rec:%s>%s' % (self.lname,lmodel.alleles[0],lmodel.alleles[1]) ]

  def estimates(self,p):
    return [0,p[self.index,0]]

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [0,c[i,i]]

  def se(self,c):
    return sqrt(self.var(c))

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1,p2 = self.estimates(p)
    e1,e2 = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1)),
            (exp(p2-a*e2),exp(p2+a*e2))]

  def __len__(self):
    return 1


class MISSING(TERM):
  def effects(self, loci, i):
    if loci[self.lname].genos[i]:
      return [1]
    else:
      return [0]

  def names(self, loci):
    if not loci[self.lname].genocount:
      return []
    return ['%s:missing' % self.lname]

  def estimates(self,p):
    return [p[self.index,0]]

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i]]

  def se(self,c):
    return [sqrt(self.var(c))]

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1, = self.estimates(p)
    e1, = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1))]

  def __len__(self):
    return 1


class NOT_MISSING(TERM):
  def effects(self, loci, i):
    if loci[self.lname].genos[i]:
      return [0]
    else:
      return [1]

  def names(self, loci):
    if not loci[self.lname].genocount:
      return []
    return ['%s:not_missing' % self.lname]

  def estimates(self,p):
    return [p[self.index,0]]

  def odds_ratios(self,p):
    return exp(self.estimates(p))

  def var(self,c):
    i = self.index
    return [c[i,i]]

  def se(self,c):
    return [sqrt(self.var(c))]

  def odds_ratio_ci(self,p,c,alpha=0.95):
    a = stats.distributions.norm.ppf( (1+alpha)/2 )
    p1, = self.estimates(p)
    e1, = self.se(c)
    return [(exp(p1-a*e1),exp(p1+a*e1))]

  def __len__(self):
    return 1


class INTERACTION(TERM):
  def __init__(self, terms):
    self.terms = terms

  def loci(self):
    results = set()
    for term in self.terms:
      results.update(term.loci())
    return sorted(results)

  def terms(self):
    results = []
    for term in self.terms:
      results.extend(term.terms())
    return results

  def effects(self, loci, i):
    results = []
    for term in self.terms:
      effects = term.effects(loci,i)
      if None in effects:
        return [None]
      if not results:
        results = list(effects)
      else:
        results = [ t1*t2 for t1 in results for t2 in effects ]
    return results

  def names(self, loci):
    results = []
    for term in self.terms:
      names = term.names(loci)

      if not results:
        results.extend(names)
      else:
        results = [ '%s x %s' % (t1,t2) for t1 in results for t2 in names ]

    return results

  def __len__(self):
    if not self.terms:
      return 0

    l = 1
    for term in self.terms:
      l *= len(term)
    return l


class COMBINATION(TERM):
  def __init__(self, terms):
    self.terms = terms

  def loci(self):
    results = set()
    for term in self.terms:
      results.update(term.loci())
    return sorted(results)

  def terms(self):
    results = []
    for term in self.terms:
      results.extend(term.terms())
    return results

  def effects(self, loci, i):
    results = []
    for term in self.terms:
      effects = term.effects(loci,i)
      if None in effects:
        return [None]
      results.extend(effects)
    return results

  def names(self, loci):
    results = []
    for term in self.terms:
      results.extend(term.names(loci))
    return results

  def __len__(self):
    return sum(len(term) for term in self.terms)


def geno_terms(model_term):
  return [ t for t in model_term.terms() if not isinstance(t, (NULL,MISSING,NOT_MISSING,NULL)) ]


termmap = { 'GENO'           : GENO,
            'GENOTYPE'       : GENO,
            'ADDDOM'         : ADOM,
            'ADOM'           : ADOM,
            'TREND'          : TREND,
            'MULTIPLICATIVE' : TREND,
            'MULT'           : TREND,
            'DOMINANT'       : DOM,
            'DOM'            : DOM,
            'RECESSIVE'      : REC,
            'REC'            : REC,
            'MISSING'        : MISSING,
            'MISS'           : MISSING,
            'NOT_MISSING'    : NOT_MISSING,
            'NOT_MISS'       : NOT_MISSING,
            'NULL'           : NULL,
          }

def get_term(name):
  try:
    return termmap[name.upper()]
  except KeyError:
    raise KeyError('Unknown genetic model term: %s' % name)


# FIXME: Needs docs+tests
class BiallelicLocusModel(object):
  def __init__(self, lname, genos, geno_indices, reference_allele=None):
    self.lname        = lname
    self.genos        = genos
    genos             = pick(genos[:], geno_indices.itervalues())
    self.allelecounts = tally(a for g in genos if g for a in g if a)
    self.genocounts   = tally(genos)
    self.genocount    = len([ 1 for g,n in self.genocounts.iteritems() if g and n ])
    self.maf          = estimate_maf(self.genocounts)
    self.alleles      = map(itemgetter(0),sorted(self.allelecounts.iteritems(),key=itemgetter(1), reverse=True))

    if len(self.alleles) != 2:
      return

    a1,a2 = self.alleles

    if reference_allele is not None:
      if reference_allele not in self.alleles:
        msg = 'Invalid reference allele %s specified for locus %s.  Alleles found: %s'
        raise ValueError(msg % (reference_allele,lname,', '.join(self.alleles)))

      if a1 != reference_allele:
        self.alleles = a1,a2 = a2,a1

    model        = genos[0].model
    self.tests   = [ model.add_genotype( (a1,a1) ),
                     model.add_genotype( (a1,a2) ),
                     model.add_genotype( (a2,a2) ) ]
    self.genomap = dict( (g,i) for i,g in enumerate(self.tests) )


class LocusModel(object):
  def __init__(self, y, X, pheno, vars, loci, model_loci):
    self.y          = y
    self.X          = X
    self.pheno      = pheno
    self.vars       = vars
    self.loci      = loci
    self.model_loci = model_loci

# FIXME: Needs docs+tests
class LocusModelBuilder(object):
  def __init__(self,locus_header,pheno_header,phenos,
                    reference_alleles=None,minmaf=0.01,mingenos=10):

    self.locus_header      = locus_header
    self.pheno_header      = pheno_header
    self.reference_alleles = reference_alleles or {}
    self.minmaf            = minmaf
    self.mingenos          = mingenos

    # Ensure phenotypes are materialized
    try:
      len(phenos)
    except TypeError:
      phenos = list(phenos)

    pidset                 = set(p[0] for p in phenos)
    self.geno_indices      = dict( (pid,i) for i,pid in enumerate(locus_header) if pid in pidset )
    self.phenos            = [ p for p in phenos if p[0] in self.geno_indices ]

  def build_model(self,term,loci):
    model_names = []
    model_terms = []
    model_loci  = {}

    for lname in set(term.loci()):
      ref    = self.reference_alleles.get(lname)
      lmodel = BiallelicLocusModel(lname,loci[lname],self.geno_indices,ref)
      if len(lmodel.alleles) != 2 or lmodel.maf < self.minmaf:
        return None
      model_loci[lname] = lmodel

    model_names = term.names(model_loci)
    k = len(model_names)

    # Reject models without all terms
    if k != len(term):
      return None

    index = 1
    for t in term.terms():
      t.index = index
      index += len(t)

    X = []
    y = []
    for row in self.phenos:
      pid  = row[0]
      stat = row[1]
      covs = row[2:]

      geno_effects = term.effects(model_loci, self.geno_indices[pid])
      if None in geno_effects:
        continue

      y.append([stat])
      X.append( [1.0] + geno_effects + covs )

    y = matrix(y, dtype=float)
    X = matrix(X, dtype=float)

    colcounts = (X.A[:,1:k+1]!=0).sum(axis=0)
    if k and colcounts.min() < self.mingenos:
      return None

    vars = ['_mean'] + model_names + self.pheno_header[2:]

    return LocusModel(y,X,self.pheno_header[1],vars,loci,model_loci)


def variable_summary(out, x, categorical_limit=5, verbose=1):
  '''
  Generate a text summary of a random variable based on the observed
  distribution.

  >>> variable_summary(sys.stdout, [1]*100, verbose=True)
  Constant = 1
  >>> variable_summary(sys.stdout, [1]*100+[0]*50, verbose=True)
  Binary variable
               0:    50   33.33%
               1:   100   66.67%
  >>> variable_summary(sys.stdout, [2]*10+[1]*100+[0]*50, verbose=True)
  Categorical variable with 3 distinct values
               0:    50   31.25%
               1:   100   62.50%
               2:    10    6.25%
  >>> variable_summary(sys.stdout, range(100)*100, verbose=True)
  Continuous variable
        min=0, mean=49.5, median=49.5, max=99, sd=28.8661
  '''
  x    = asarray(x).ravel()
  dist = tally(x)
  n    = len(dist)

  if not n:
    raise ValueError('No values in phenotype distribution')
  elif n == 1:
    out.write('Constant = %g\n' % dist.popitem()[0])
  elif n <= categorical_limit:
    if n == 2:
      out.write('Binary variable\n')
    else:
      out.write('Categorical variable with %d distinct values\n' % n)

    if verbose:
      total = sum(dist.itervalues())
      for v,n in sorted(dist.iteritems()):
        out.write('         %5g: %5d  %6.2f%%\n' % (v,n,float(n)/total*100))
  else:
    out.write('Continuous variable\n')
    out.write('      min=%g, mean=%g, median=%g, max=%g, sd=%g\n' % \
                                 (x.min(),x.mean(),median(x),x.max(),x.std()))


# FIXME: Needs a better name
# FIXME: Needs docs+tests
def print_results(out,locus_model,linear_model,verbose=1):
  b     = linear_model.beta.T
  stde  = linear_model.W.diagonal().A**.5
  z     = b.A/stde
  oddsr = exp(b)
  p     = 2*stats.distributions.norm.sf(abs(z))

  vars = locus_model.vars or linear_model.vars or \
         [ 'Covariate_%02d' % i for i in xrange(model.X.shape[1]) ]

  n = len(vars)

  if len(linear_model.categories) > 2:
    out.write('Multinomial logistic regression\n')
  else:
    out.write('Logistic regression\n')

  out.write('Observations = %d\n' % len(linear_model.X))
  for c in linear_model.categories:
    out.write('  CATEGORY %2d: %5d\n' % (c,(linear_model.y_ord==c).sum()))
  out.write('Covariate summary:\n')
  for name,x in izip(vars,linear_model.X.T):
    out.write('  %-12s: ' % name)
    variable_summary(out, x, verbose=verbose)
  out.write('\n')
  out.write('Log likelihood = %f\n' % linear_model.L)
  out.write('\n')

  if len(linear_model.categories) <= 2:
    out.write('Variable                       Coef.   Std.Err.     OR        Z      Pr > |Z| \n')
    out.write('------------------------ ----------- ---------- ---------- ------- -----------\n')
    for k,var in enumerate(vars):
      pv = format_pvalue(p[0,k],10)
      out.write('%-25s %10.6f %10.6f %10.6f %7.2f %s\n' % (var,b[0,k],stde[0,k],oddsr[0,k],z[0,k],pv))
    out.write('\n')
  else:
    out.write('TYPE       Variable                  Coef.   Std.Err.     OR        Z      Pr > |Z| \n')
    out.write('---------- ------------------- ----------- ---------- ---------- ------- -----------\n')
    for i,cat in enumerate(linear_model.categories[1:]):
      for j,var in enumerate(vars):
        k = i*n+j
        pv = format_pvalue(p[0,k],10)
        out.write('%-10s %-20s %10.6f %10.6f %10.6f %7.2f %s\n' % (cat,var,b[0,k],stde[0,k],oddsr[0,k],z[0,k],pv))
      out.write('\n')

  out.write('\n')


# FIXME: Needs a better name
# FIXME: Needs docs+tests
def print_results_linear(out,locus_model,linear_model,verbose=1):
  y     = linear_model.y
  b     = linear_model.beta.T
  stde  = (linear_model.ss*linear_model.W.diagonal()).A**0.5
  t     = b.A/stde
  oddsr = exp(b)
  n,m   = linear_model.X.shape
  p     = 2*stats.distributions.t.sf(abs(t),n-m)

  vars = locus_model.vars or linear_model.vars or \
         [ 'Covariate_%02d' % i for i in xrange(model.X.shape[1]) ]

  out.write('Linear (Gaussian) regression\n')
  out.write('Observations = %d\n' % n)
  out.write('Phenotype distribution:\n')
  out.write('  min=%g, mean=%g, median=%g, max=%g, sd=%g\n' % \
                     (y.min(),y.mean(),median(y),y.max(),y.std()))
  out.write('Covariate summary:\n')
  for name,x in izip(vars,linear_model.X.T):
    out.write('  %-12s: ' % name)
    variable_summary(out,x, verbose=verbose)
  out.write('Residual variance = %f\n' % linear_model.ss)
  out.write('\n')
  out.write('Variable                  Coef.   Std.Err.    T     Pr > |T| \n')
  out.write('------------------- ----------- ---------- ------- ----------\n')
  for i,var in enumerate(vars):
    pv = format_pvalue(p[0,i])
    out.write('%-20s %10.6f %10.6f %7.2f %s\n' % (var,b[0,i],stde[0,i],t[0,i],pv))
  out.write('\n')


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
