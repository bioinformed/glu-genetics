# -*- coding: utf-8 -*-

__abstract__  = 'library for association testing'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   math               import ceil,log
from   itertools          import izip,chain
from   collections        import defaultdict
from   operator           import itemgetter

from   numpy              import array,matrix,asarray,asanyarray,zeros, \
                                 exp,nan,abs,arange,median,inf,minimum
from   scipy              import stats

from   glu.lib.utils      import tally,as_set,is_str
from   glu.lib.fileutils  import namefile,list_reader,map_reader,table_reader,table_columns,resolve_column_headers,tryint1,\
                                 parse_augmented_name,resolve_column_header_atom
from   glu.lib.genolib    import load_genostream,pick
from   glu.lib.genolib.transform import _union_options, _intersect_options
from   glu.lib.formula    import INTERCEPT,NO_INTERCEPT,GENOTERM,PHENOTERM,COMBINATION, \
                                 GENO,TREND,FormulaParser


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
  '''
  >>> t,df,ors = contingency_analysis( [[1,2,3],[3,3,1],[1,4,1]],mincell=0,minmargin=1)
  >>> t
  3.7798941798941801
  >>> df
  4
  >>> ors
  array([[ 0.5       ,  0.11111111],
         [ 2.        ,  0.33333333]])
  '''
  c = asanyarray(c,dtype=int)

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


def table_trend(x):
  '''
  Cochrane-Armitage Trend Test for 2xc tables

  @param  x: 2xc table of counts
  @type   x: 2d array or list of lists
  @return  : Cochrane-Armitage trend test statistic (chi-squared w/ 1 df)
  @rtype   : float

  >>> x = [[26,26,23,18,9],[6,7,9,14,23]]
  >>> table_trend(x)
  22.961467763890997
  '''
  x = asanyarray(x, dtype=int)

  if len(x.shape) != 2 or x.shape[0] != 2:
    raise ValueError('tabletrend requires a 2xc table')

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
def load_phenos(filename,pid=0,pheno=1,columns=None,deptype=int,categorical=None,includevar=None,
                         excludevar=None,allowdups=False,verbose=1,errs=sys.stderr):
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
  phenos = table_reader(filename,want_header=True)

  try:
    header = strip_trailing_empty(phenos.next())
  except StopIteration:
    raise ValueError('Empty phenotype file')

  if categorical:
    header,phenos = create_all_categorical(header,phenos,categorical)

  if includevar or excludevar:
    header,phenos = subset_all_variables(header,phenos,includevar,excludevar)

  indices = resolve_column_headers(header,[pid,pheno])
  if columns is not None:
    indices.extend(resolve_column_headers(header,columns))
  else:
    used = set(indices)
    indices.extend(i for i in range(len(header)) if i not in used)

  if indices != range(len(header)):
    phenos = table_columns(phenos, columns=indices, header=header, want_header=True)
    header = phenos.next()

  # Assign default names to missing headers to ensure that output code
  # doesn't get cranky
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
    if verbose > 0:
      errs.write('[ERROR] Duplicate headers detected (n=%d)' % len(dups))
    if verbose > 1:
      for h,n in sorted(dups, key=itemgetter(1,0), reverse=True):
        errs.write('         %-12s : %2d\n' % (h,n))
    raise ValueError('Invalid repeated heading in phenotype header')

  def _phenos():
    warn_msg = '[WARNING] Subject "%s" dropped on line #%04d due to missing data in column "%s"\n'
    note_msg = '[NOTICE] Read phenotypes from %s: ' \
               '%d covariates, %d valid records, %d records dropped due to missing data, %d distinct subjects\n'

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
      if verbose > 0:
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


def subset_variable(header,data,variable,include=None,exclude=None):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = subset_variable(header,data,'a',include='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']

  >>> h,d = subset_variable(header,data,'b',exclude='?')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_variable(header,data,'c',include='')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  '''
  index = resolve_column_header_atom(header,variable)

  if include is not None and exclude is not None:
    include = as_set(include)-as_set(exclude)
    exclude = None

  if include is not None:
    def _subset():
      for row in data:
        if row[index] in include:
          yield row

  elif exclude is not None:
    def _subset():
      for row in data:
        if row[index] not in exclude:
          yield row

  return header,_subset()


def subset_all_variables(header,data,include=None,exclude=None):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = subset_all_variables(header,data,include='a=1,2')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_all_variables(header,data,exclude='b=?')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_all_variables(header,data,include='c')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']

  >>> h,d = subset_all_variables(header,data,include='c=')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_all_variables(header,data,exclude='c')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']

  >>> h,d = subset_all_variables(header,data,exclude='c=')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['3', '?', '1']
  '''
  if is_str(include):
    include = [include]
  else:
    include = include or []

  if is_str(exclude):
    exclude = [exclude]
  else:
    exclude = exclude or []

  for invar in include:
    if '=' not in invar:
      header,data = subset_variable(header,data,invar,exclude='')
    else:
      var,values = invar.split('=',1)
      values = set(v.strip() for v in values.split(','))
      header,data = subset_variable(header,data,var,include=values)

  for exvar in exclude:
    if '=' not in exvar:
      header,data = subset_variable(header,data,exvar,include='')
    else:
      var,values = exvar.split('=',1)
      values = set(v.strip() for v in values.split(','))
      header,data = subset_variable(header,data,var,exclude=values)

  return header,data


def create_categorical(header,data,variable,prefix=None,ref=None,include=None,exclude=None,
                              yes='1',no='0',missing=''):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = create_categorical(header,data,'a')
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0']
  ['2', 'F', '', '0', '1', '0']
  ['3', '?', '1', '0', '0', '1']

  >>> h,d = create_categorical(header,data,'a',ref=['1'])
  >>> h
  ['a', 'b', 'c', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '0', '0']
  ['2', 'F', '', '1', '0']
  ['3', '?', '1', '0', '1']

  >>> h,d = create_categorical(header,data,'b',prefix='',exclude=['?'],yes='Y',no='N')
  >>> h
  ['a', 'b', 'c', 'F', 'M']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', 'N', 'Y']
  ['2', 'F', '', 'Y', 'N']
  ['3', '?', '1', '', '']

  >>> h,d = create_categorical(header,data,'c',ref='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  ['3', '?', '1']

  >>> h,d = create_categorical(header,data,'c',exclude='1')
  >>> h
  ['a', 'b', 'c']
  >>> for row in d:
  ...   print row
  ['1', 'M', '']
  ['2', 'F', '']
  ['3', '?', '1']
  '''
  index = resolve_column_header_atom(header,variable)
  if include is not None and exclude is not None:
    include = as_set(include)-as_set(exclude)
    exclude = None

  if not isinstance(data,(tuple,list)):
    data = list(data)

  values = set(row[index] for row in data if len(row)>=index)
  values.discard('')

  if include is not None:
    values &= as_set(include)
  if exclude is not None:
    values -= as_set(exclude)

  if prefix is None:
    prefix = '%s_' % variable

  if ref is None:
    ref = set()
  else:
    ref = as_set(ref)

  values = sorted(values-ref)

  if not values:
    return header,data

  header = header + [ '%s%s' % (prefix,value) for value in sorted(values) ]
  values = dict( (v,i) for i,v in enumerate(values) )

  def _make_category():
    n = len(values)
    missingrow = [missing]*n
    for row in data:
      if index>=len(row):
        yield row+missingrow
        continue

      val = row[index]
      cats = [no]*n
      if val in values:
        cats[ values[val] ] = yes
      elif val not in ref:
        cats = missingrow

      yield row+cats

  return header,_make_category()


def create_all_categorical(header,phenos,categorical):
  '''
  >>> header =  ['a','b','c']
  >>> data   = [['1','M','' ],
  ...           ['2','F','' ],
  ...           ['3','?','1']]

  >>> h,d = create_all_categorical(header,data,['a'])
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0']
  ['2', 'F', '', '0', '1', '0']
  ['3', '?', '1', '0', '0', '1']

  >>> h,d = create_all_categorical(header,data,['a','b:ref=M:prefix=:exclude=?','c:exclude=1:yes=Y:no=N'])
  >>> h
  ['a', 'b', 'c', 'a_1', 'a_2', 'a_3', 'F']
  >>> for row in d:
  ...   print row
  ['1', 'M', '', '1', '0', '0', '0']
  ['2', 'F', '', '0', '1', '0', '1']
  ['3', '?', '1', '0', '0', '1', '']
  '''
  allowed_args = set(['prefix','ref','include','exclude','missing','yes','no'])

  for cat in categorical:
    opts = {}
    var  = parse_augmented_name(cat,opts)

    illegal = set(opts) - allowed_args
    if illegal:
      raise ValueError('Illegal argument(s) to categorical: %s' % ','.join(sorted(illegal)))

    for arg,val in opts.items():
      if arg in ('ref','include','exclude'):
        opts[arg] = set(v.strip() for v in val.split(','))

    header,phenos = create_categorical(header,phenos,var,**opts)

  return header,phenos


def _load_loci(filename,options,keep):
  if options.includesamples:
    keep &= _union_options(options.includesamples)

  if options.excludesamples:
    keep -= _intersect_options(options.excludesamples)

  loci = None
  if filename is not None:
    loci = load_genostream(filename,format=options.informat,genorepr=options.ingenorepr,
                           genome=options.loci,phenome=options.pedigree,
                           transform=options).as_ldat()
    keep &= set(loci.samples)

  fixedloci = None
  if options.fixedloci:
    fixedloci = load_genostream(options.fixedloci,format=options.informat,genorepr=options.ingenorepr,
                                genome=options.loci,phenome=options.pedigree,
                                transform=options).as_ldat()
    keep     &= set(fixedloci.samples)
    samples   = loci.samples if loci is not None else None
    fixedloci = fixedloci.transformed(include_samples=keep,order_samples=samples)

  if fixedloci and loci:
    assert fixedloci.samples == loci.samples

  if loci is None:
    return None,fixedloci,list(fixedloci.samples)

  loci        = loci.transformed(include_samples=keep)
  samples     = loci.samples

  return loci,fixedloci,loci.samples


def parse_formulae(options,models):
  scan = options.scan = options.scan or 'locus'

  if options.test is not None:
    _,options.test = FormulaParser().parse(options.test)

  if options.display is not None:
    _,options.display = FormulaParser().parse(options.display)

  if options.model and scan not in options.model.loci():
    raise ValueError('Model does not contain any genotype scan terms')

  if options.model and not options.test:
    options.test = COMBINATION(t for t in options.model.terms() if scan in t.loci())

  if not options.test:
    options.test = GENO(scan)

  if scan not in options.test.loci():
    raise ValueError('Test does not contain any genotype scan terms')

  if not options.model:
    phenos = COMBINATION( PHENOTERM(pheno) for pheno in models.pheno_header[2:] )
    options.model = options.test + phenos

  if not options.display:
    options.display = options.test

  try:
    options.test = options.model.find(options.test)
  except KeyError:
    raise ValueError('Formula does not contain all terms to be tested')

  try:
    options.display = options.model.find(options.display)
  except KeyError:
    raise ValueError('Formula does not contain all terms to display')

  options.null = COMBINATION(t for t in options.model.terms()
                                if t not in options.test and scan not in t.loci())

  options.stats = set(t.strip().lower() for t in options.stats.split(','))
  options.stats.discard('')
  extra = options.stats - set(['score','wald','lrt'])
  if extra:
    raise ValueError('Unknown test(s) specified: %s' % ','.join(sorted(extra)))


def build_models(phenofile, genofile, options, deptype=int, errs=sys.stderr):
  warn_msg = '[WARNING] Subject "%s" excluded from analysis\n'

  if options.model is not None:
    pheno,options.model = FormulaParser().parse(options.model)
    if pheno is not None:
      options.pheno = pheno
    covs = set(t.name for t in options.model.expand_terms() if isinstance(t,PHENOTERM))
  else:
    pheno  = covs = None
    covs = None

  verbose       = options.verbose
  header,phenos = load_phenos(phenofile,pid=options.pid,pheno=options.pheno,columns=covs,deptype=deptype,
                                        categorical=options.categorical,includevar=options.includevar,
                                        excludevar=options.excludevar,allowdups=options.allowdups,
                                        verbose=verbose,errs=errs)
  phenos        = list(phenos)
  subjects      = set(p[0] for p in phenos)
  keep          = subjects.copy()
  phenocount1   = len(phenos)

  loci,fixedloci,samples = _load_loci(genofile,options,keep)

  if subjects != keep:
    phenos = [ p for p in phenos if p[0] in keep ]

    if verbose > 0:
      errs.write('[NOTICE] After exclusions, %d subjects remain, %d subjects excluded\n' % (len(phenos),len(subjects)-len(keep)))
    if verbose > 1:
      for pid in sorted(subjects-keep):
        errs.write(warn_msg % pid)

  reference_alleles = map_reader(options.refalleles) if options.refalleles else None

  models = LocusModelBuilder(samples,header,phenos,
                             reference_alleles=reference_alleles,
                             minmaf=options.minmaf,mingenos=options.mingenos)

  parse_formulae(options,models)

  fixedloci = dict(fixedloci) if fixedloci is not None else {}

  # Collect unbound genotype terms and check fixed loci
  fixed = set()
  gterms = []
  for t in options.model.expand_terms():
    if isinstance(t, GENOTERM):
      if t.name==options.scan:
        gterms.append(t)
      elif t.name not in fixedloci:
        raise ValueError('Locus %s used in model, but not found in fixed loci' % t.name)
      else:
        fixed.add(t.name)

  # Remove unused terms in fixed loci
  for lname in list(fixedloci):
    if lname not in fixed:
      del fixedloci[lname]

  if not gterms:
    raise ValueError('No genotype terms to scan')

  return loci,fixedloci,gterms,models


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

# FIXME: Needs docs+tests
class BiallelicLocusModel(object):
  def __init__(self, lname, genos, geno_indices, reference_allele=None):
    self.lname        = lname
    self.genos        = genos
    genos             = pick(genos, geno_indices.itervalues())
    self.allelecounts = tally(a for g in genos if g for a in g if a)
    self.genocounts   = tally(genos)
    self.genocount    = len([ 1 for g,n in self.genocounts.iteritems() if g and n ])
    self.maf          = estimate_maf(self.genocounts)
    self.alleles      = map(itemgetter(0),sorted(self.allelecounts.iteritems(),key=itemgetter(1),reverse=True))

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
    self.counts  = [ self.genocounts.get(g,0) for g in self.tests ]
    self.genomap = dict( (g,i) for i,g in enumerate(self.tests) )


class LocusModel(object):
  def __init__(self, formula, y, X, pheno, vars, loci, model_loci, pids):
    self.formula    = formula
    self.y          = y
    self.X          = X
    self.pheno      = pheno
    self.vars       = vars
    self.loci       = loci
    self.model_loci = model_loci
    self.pids       = pids


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

    pidset            = set(p[0] for p in phenos)
    self.geno_indices = dict( (pid,i) for i,pid in enumerate(locus_header) if pid in pidset )
    self.phenos       = [ p for p in phenos if p[0] in self.geno_indices ]

  def build_model(self,term,loci,mingenos=None):
    genoterms    = []
    phenoterms   = []
    interactions = []
    intercept    = 0
    no_intercept = 0

    for t in term.terms():
      if isinstance(t, INTERCEPT):
        intercept += 1
      elif isinstance(t, NO_INTERCEPT):
        no_intercept += 1
      elif isinstance(t, GENOTERM):
        genoterms.append(t)
      elif isinstance(t, PHENOTERM):
        phenoterms.append(t)
      else:
        interactions.append(t)

    term = COMBINATION()
    if not no_intercept or intercept:
      term += INTERCEPT()
    for t in chain(genoterms,phenoterms,interactions):
      term += t

    index = 0
    geno_indices = []
    for t in term.terms():
      n = len(t)
      if n:
        t.index = index
        if t.loci():
          geno_indices.extend( range(index,index+n) )
        index += n

    for t in term.expand_terms():
      if isinstance(t,PHENOTERM):
        try:
          t.pindex = self.pheno_header.index(t.name)
        except (ValueError,IndexError):
          # FIXME: Index refers to file header, not pheno header.  More so,
          # you can't specify 0 or 1 because those are taken as intercept
          # terms.
          try:
            t.pindex = tryint1(t.name)
            t.name   = self.pheno_header[t.pindex]
          except (ValueError,TypeError):
            raise ValueError("Cannot find phenotype column '%s'" % t.name)

    model_names = []
    model_terms = []
    model_loci  = {}

    for lname in set(term.loci()):
      ref    = self.reference_alleles.get(lname)
      lmodel = BiallelicLocusModel(lname,loci[lname],self.geno_indices,ref)
      if len(lmodel.alleles) != 2 or lmodel.maf < self.minmaf:
        return None
      model_loci[lname] = lmodel

    model_names = term.term_names(model_loci)
    k = len(model_names)

    # Reject models without all terms
    if k != len(term):
      return None

    X = []
    y = []
    pids = []
    for row in self.phenos:
      pid  = row[0]
      stat = row[1]

      effects = term.effects(model_loci, row, self.geno_indices[pid])
      if None in effects:
        continue

      pids.append(pid)
      y.append([stat])
      X.append(effects)

    if not y:
      return None

    # Do not fit models with no contrast
    if len(set(i[0] for i in y)) < 2:
      return None

    y = matrix(y, dtype=float)
    X = matrix(X, dtype=float)

    if mingenos is None:
      mingenos = self.mingenos

    # FIXME: What does mingenos mean in a world with arbitrary formulae?
    if geno_indices and mingenos:
      colcounts0 = (X.A[:,geno_indices]==0).sum(axis=0)
      colcounts1 = (X.A[:,geno_indices]!=0).sum(axis=0)
      colcounts  = minimum(colcounts0, colcounts1)
      if colcounts.min() < self.mingenos:
        return None

    return LocusModel(term,y,X,self.pheno_header[1],model_names,loci,model_loci,pids)


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
  p     = 2*stats.distributions.norm.cdf(-abs(z))

  vars = locus_model.vars or linear_model.vars or \
         [ 'Covariate_%02d' % i for i in xrange(linear_model.X.shape[1]) ]

  n = len(vars)

  if len(linear_model.categories) > 2:
    out.write('Multinomial logistic regression\n')
  else:
    out.write('Logistic regression\n')

  out.write('Observations = %d\n' % len(linear_model.X))
  for c in linear_model.categories:
    out.write('  CATEGORY %2d: %5d\n' % (c,(linear_model.y==c).sum()))
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
  p     = 2*stats.distributions.t.cdf(-abs(t),n-m)

  vars = locus_model.vars or linear_model.vars or \
         [ 'Covariate_%02d' % i for i in xrange(linear_model.X.shape[1]) ]

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
