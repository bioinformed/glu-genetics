import os
import time

from   operator            import itemgetter
from   itertools           import islice,izip,groupby,imap
from   collections         import defaultdict

import numpy
import tables

from   biozilla.utils      import tally
from   biozilla.genoarray  import snp_acgt,snp_acgt2
from   biozilla.genoarray2 import UnphasedMarkerModel,GenotypeArrayDescriptor,GenotypeArray,Genotype


def ilen(s):
  i = 0
  for x in s:
    i+=1
  return i


def xlen(s):
  return s[0],ilen(s[1])


def save_genotriples_binary(filename,triples,chunksize=232960):
  '''
  Write the genotype triple data to file.

  @param filename: a file name or file object
  @type  filename: str or file object
  @param      triples: genotype triple data
  @type       triples: sequence
  '''

  gfile = tables.openFile(filename,mode='w')
  gfile.root._v_attrs.format = 'genotriple'

  class TripleDesc(tables.IsDescription):
    sample = tables.Int32Col(pos=0)
    locus  = tables.Int32Col(pos=1)
    geno   = tables.UInt8Col(pos=2)

  filters = tables.Filters(complevel=5, complib='zlib',shuffle=True,fletcher32=True)
  #filters = None # tables.Filters(complevel=5, complib='zlib',shuffle=True)
  genotypes = gfile.createTable(gfile.root, 'genotypes', TripleDesc,
                            'Sequence of encoded sample, locus, genotype triples',
                            filters=filters, chunkshape=(chunksize//4,),expectedrows=5000000)

  samplemap = {}
  locusmap  = {}

  sd = samplemap.setdefault
  sl = samplemap.__len__
  ld = locusmap.setdefault
  ll = locusmap.__len__

  row = genotypes.row
  for sample,locus,geno in triples:
    row['sample'] = sd(sample,sl())
    row['locus']  = ld(locus,ll())
    row['geno']   = geno.index
    row.append()

  genotypes.flush()

  samples = map(itemgetter(0), sorted(samplemap.iteritems(), key=itemgetter(1)))
  loci    = map(itemgetter(0), sorted(locusmap.iteritems(),  key=itemgetter(1)))

  gfile.createArray(gfile.root, 'samples', samples).flush()
  gfile.createArray(gfile.root, 'loci',    loci).flush()

  gfile.close()


def load_genotriples_binary(filename,limit=None):
  '''
  Load genotype triples from file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of genotypes loaded
  @type         limit: int or None
  @return: sequence of tuples of sample name, locus name, and genotype representation
  @rtype:  generator
  '''
  gfile = tables.openFile(filename,mode='r')
  format = gfile.root._v_attrs.format

  if format != 'genotriple':
    raise ValieError, 'Unknown format: %s' % format

  samples = map(str,gfile.root.samples[:])
  loci    = map(str,gfile.root.loci[:])
  genos   = snp_acgt2.model.genotypes

  for x in gfile.root.genotypes:
    yield samples[x[0]],loci[x[1]],genos[x[2]]

  gfile.close()


def save_strings(gfile,name,data,filters=None):
  n = max(len(s) for s in data)
  a = gfile.createCArray(gfile.root, name, tables.StringAtom(itemsize=n),
                         (len(data),), filters=filters)
  a[:] = data
  a.flush()


def save_models(gfile, models):
  allelemap = {}
  ad = allelemap.setdefault
  al = allelemap.__len__

  modelmap = {}

  filters = tables.Filters(complevel=5, complib='zlib', shuffle=True, fletcher32=True)

  #####
  # WRITE LOCUS MODELS: vector of locus -> model index
  #
  # Collect and collapse redundant models and write index array
  # also, collect an index array of alleles to write later
  class LocusModelDesc(tables.IsDescription):
    model = tables.Int32Col(pos=0)

  locus_models = gfile.createTable(gfile.root, 'locus_models', LocusModelDesc, 'locus models',
                                       filters=filters, expectedrows=len(models))

  locus_row = locus_models.row
  for model in models:
    genotypes = (model.max_alleles,model.allow_hemizygote)+tuple(model.genotypes[1:])
    index = modelmap.get(genotypes)
    if index is None:
      index = modelmap[genotypes] = len(modelmap)
      for allele in model.alleles:
        ad(allele,al())

    locus_row['model'] = index
    locus_row.append()

  locus_models.flush()

  # Smash modelmap down to an ordered list of tuples
  models = map(itemgetter(0), sorted(modelmap.iteritems(),  key=itemgetter(1)))

  #####
  # WRITE MODELS: sequence of model max_alleles and allow_hemizygote parameters
  #
  # Used to re-construct model objects
  class ModelDesc(tables.IsDescription):
    max_alleles      = tables.UInt16Col(pos=0)
    allow_hemizygote = tables.UInt16Col(pos=1)

  mods = gfile.createTable(gfile.root, 'models', ModelDesc, 'models',
                                       filters=filters, expectedrows=len(models))

  model_row = mods.row
  for model in models:
    model_row['max_alleles']      = model[0]
    model_row['allow_hemizygote'] = model[1]
    model_row.append()

  mods.flush()

  #####
  # WRITE MODEL_GENOTYPES: model -> allele1/allele2
  #
  # Used to re-construct model objects.  Ordered list of genotypes per model

  class GenotypeDesc(tables.IsDescription):
    model   = tables.Int32Col(pos=0)
    allele1 = tables.Int32Col(pos=1)
    allele2 = tables.Int32Col(pos=2)

  genos = gfile.createTable(gfile.root, 'model_genotypes', GenotypeDesc, 'genotypes in each model',
                                        filters=filters)

  geno_row = genos.row
  for i,model in enumerate(models):
    for allele1,allele2 in model[2:]:
      geno_row['model']   = i
      geno_row['allele1'] = allelemap[allele1]
      geno_row['allele2'] = allelemap[allele2]
      geno_row.append()

  genos.flush()

  #####
  # WRITE MODEL_ALLELES: sequence of allele strings
  #
  # Used to re-construct model objects.  Ordered list of all possible alleles
  alleles = map(itemgetter(0), sorted(allelemap.iteritems(),  key=itemgetter(1)))
  alleles[0] = ''
  save_strings(gfile,'model_alleles',alleles,filters=filters)


def load_models(gfile):
  alleles         = map(str,gfile.root.model_alleles[:])
  alleles[0]      = None
  mods            = list(gfile.root.models[:])
  model_genotypes = list(gfile.root.model_genotypes[:])
  model_genotypes = groupby(model_genotypes, itemgetter(0))

  models = []
  for mod,(i,mgenos) in izip(mods,model_genotypes):
    model = UnphasedMarkerModel(mod[1],mod[0])
    for j,allele1,allele2 in mgenos:
      allele1 = alleles[allele1] or None
      allele2 = alleles[allele2] or None
      model.add_genotype( (allele1,allele2) )
    models.append(model)

  locus_models = [ m[0] for m in gfile.root.locus_models[:] ]
  return [ models[i] for i in locus_models ]


def save_genomatrix_binary(filename,matrix,format,compress=True,scratch=16*1024*1024):
  '''
  Write the genotype matrix data to file.

  @param     filename: a file name or file object
  @type      filename: str or file object
  @param       matrix: the genotype matrix data
  @type        matrix: sequence
  @param       format: text string output in the first header field to
                       indicate data format (default is blank)
  @type        format: string
  '''
  if format not in ('ldat','sdat'):
    raise ValueError('format must be either ldat or sdat')

  matrix  = iter(matrix)
  header  = matrix.next()
  l1,row1 = matrix.next()
  row1    = row1
  n       = len(row1.data)

  gfile = tables.openFile(filename,mode='w')
  gfile.root._v_attrs.format = format

  crows = min(max(8, int(scratch//n)),8192)
  ccols = min(n,8192)

  if compress:
    filters = tables.Filters(complevel=5, complib='zlib', shuffle=(format=='sdat'),fletcher32=True)
  else:
    filters = tables.Filters(fletcher32=True)

  genos = gfile.createEArray(gfile.root, 'genotypes', tables.UInt8Atom(), (0,n),
                             'Matrix of binary encoded genotypes values',
                             chunkshape=(crows,ccols), filters=filters, expectedrows=50000)

  rowlabels = [l1]
  chunk     = [row1.data]

  if format == 'sdat':
    models = row1.descriptor.models
    for label,row in matrix:
      rowlabels.append(label)
      chunk.append(row.data)
      if len(chunk) >= crows:
        genos.append(chunk)
        chunk = []

  elif format == 'ldat':
    models = [row1.descriptor.models[0]]
    for label,row in matrix:
      rowlabels.append(label)
      chunk.append(row.data)
      models.append(row.descriptor.models[0])
      if len(chunk) >= crows:
        genos.append(chunk)
        chunk = []

  if chunk:
    genos.append(chunk)
    chunk = []

  genos.flush()

  save_strings(gfile, 'rows', rowlabels, filters)
  save_strings(gfile, 'cols', header,    filters)
  save_models(gfile, models)

  gfile.close()


def load_genomatrix_binary(filename,format,limit=None,chunksize=4096,scratch=32*1024*1024):
  '''
  Load the genotype matrix data from file.
  Note that the first row is header and the rest rows are genotypes,
  and the file is tab delimited.

  @param     filename: a file name or file object
  @type      filename: str or file object
  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param        limit: limit the number of columms loaded
  @type         limit: int or None
  @param       unique: verify that rows and columns are uniquely labeled
                       (default is True)
  @type        unique: bool
  @return:             format and sequence of column names followed by
                       tuples of row label and row data
  @rtype:              tuple of string and generator
  '''
  gfile = tables.openFile(filename,mode='r')
  format_found = gfile.root._v_attrs.format

  # 'key' and blank are always allowed for backward compatibility
  if format_found == 'key':
    format_found = ''

  if format not in ('ldat','sdat'):
    raise ValieError, 'Unknown format: %s' % format

  cols   = map(intern,map(str,gfile.root.cols[:]))
  rows   = map(intern,map(str,gfile.root.rows[:]))
  models = list(load_models(gfile))

  if format == format_found == 'sdat':
    def _gen_load_genomatrix():
      yield cols

      assert len(models) == len(cols)
      descr = GenotypeArrayDescriptor(models)

      chunksize = max(2, int(scratch//gfile.root.genotypes.rowsize))
      chunks    = int(len(rows)+chunksize-1)//chunksize

      stop = 0
      for i in range(chunks):
        start,stop = stop,stop+chunksize
        labels = rows[start:stop]
        chunk  = gfile.root.genotypes[start:stop,:]
        for j,label in enumerate(labels):
          g = GenotypeArray(descr)
          g.data = chunk[j,:]
          yield label,g

      gfile.close()

  elif format == format_found == 'ldat':
    def _gen_load_genomatrix():
      yield cols

      assert len(models) == len(rows)
      chunksize = max(2, int(scratch//gfile.root.genotypes.rowsize))
      chunks    = int(len(rows)+chunksize-1)//chunksize

      stop = 0
      mods = iter(models)
      for i in range(chunks):
        start,stop = stop,stop+chunksize
        labels = rows[start:stop]
        chunk  = gfile.root.genotypes[start:stop,:]
        for j,label in enumerate(labels):
          descr = GenotypeArrayDescriptor( [mods.next()]*len(cols) )
          g = GenotypeArray(descr)
          g.data = chunk[j,:]
          yield label,g

      gfile.close()

  if format == 'ldat' and format_found == 'sdat':
    def _gen_load_genomatrix():
      yield rows

      descr = GenotypeArrayDescriptor(models)

      chunkrows,chunkcols = gfile.root.genotypes.chunkshape
      chunksize = max(1,int(scratch/(chunkcols*len(rows))))*chunkcols
      chunkbits = chunksize*8
      chunks    = int((gfile.root.genotypes.rowsize+chunksize-1)//chunksize)

      stopbit = 0
      stop    = 0
      mods = iter(models)
      for i in range(chunks):
        start    = stop
        startbit = stopbit

        # Note: O(N) sequential search.  This could be done via binary search
        while (stopbit-startbit) < chunkbits and stop < len(cols):
          stop   += 1
          stopbit = descr.offsets[stop]

        labels     = cols[start:stop]
        startbyte  = int(startbit//8)
        stopbyte   = int((stopbit+7)//8)
        offset     = int(startbit%8)
        chunk      = gfile.root.genotypes[:,startbyte:stopbyte]
        chunkdescr = GenotypeArrayDescriptor(models[start:stop],initial_offset=offset)

        chunkgenos = []
        for j in xrange(len(rows)):
          g = GenotypeArray(chunkdescr)
          g.data = chunk[j,:]
          chunkgenos.append(g[:])

        for j,label in enumerate(labels):
          coldescr = GenotypeArrayDescriptor( [mods.next()]*len(rows) )
          g = GenotypeArray(coldescr, imap(itemgetter(j), chunkgenos))
          yield label,g

      gfile.close()

  elif format == 'sdat' and format_found == 'ldat':
    def _gen_load_genomatrix():
      yield rows

      assert len(models) == len(rows)
      descr = GenotypeArrayDescriptor(models)

      chunkrows,chunkcols = gfile.root.genotypes.chunkshape
      chunksize = max(1,int(scratch/(chunkcols*len(rows))))*chunkcols
      chunkbits = chunksize*8
      chunks    = int((gfile.root.genotypes.rowsize+chunksize-1)//chunksize)

      stopbit = 0
      stop    = 0
      for i in range(chunks):
        start    = stop
        startbit = stopbit

        # Note: O(N) sequential search.  This could be done via binary search
        while (stopbit-startbit) < chunkbits and stop < len(cols):
          stop   += 1
          stopbit = descr.offsets[stop]

        labels     = cols[start:stop]
        startbyte  = int(startbit//8)
        stopbyte   = int((stopbit+7)//8)
        offset     = int(startbit%8)
        chunk      = gfile.root.genotypes[:,startbyte:stopbyte]

        chunkgenos = []
        for j in xrange(len(rows)):
          chunkdescr = GenotypeArrayDescriptor([models[j]]*(stop-start),initial_offset=offset)
          g = GenotypeArray(chunkdescr)
          g.data = chunk[j,:]
          chunkgenos.append(g[:])

        for j,label in enumerate(labels):
          g = GenotypeArray(descr, imap(itemgetter(j), chunkgenos))
          yield label,g

      gfile.close()

  return format_found or format, _gen_load_genomatrix()


def recode_genomatrix_bitpacked(genos, format, old_genorepr):
  '''
  Returns a new genomatrix with the genotypes recoded to a new internal representation

  @param        genos: genomatrix
  @type         genos: genomatrix generator
  @param old_genorepr: internal representation of genotypes to be transformed from
  @type  old_genorepr: UnphasedMarkerRepresentation or similar object
  @param     new_repr: internal representation of genotypes to be transformed to
  @type      new_repr: UnphasedMarkerRepresentation or similar object
  @return            : a genomatrix with a packed internal format
  @rtype             : genomatrix generator

  >>> genos = [('s1','s2','s3'),
  ...          ('l1',[17,0,51]),
  ...          ('l2',[0,0,0]),
  ...          ('l3',[17,0,0]),
  ...          ('l4',[52,0,68])]
  >>> for row in recode_genomatrix(genos,snp_acgt,snp_marker):
  ...   print row
  ('s1', 's2', 's3')
  ('l1', [('A', 'A'), None, ('G', 'G')])
  ('l2', [None, None, None])
  ('l3', [('A', 'A'), None, None])
  ('l4', [('G', 'T'), None, ('T', 'T')])
  '''
  if format == 'ldat':
    genos  = iter(genos)
    header = genos.next()
    yield header

    modelmap = {}

    for i,(label,row) in enumerate(genos):
      geno = [ g or (None,None) for g in old_genorepr.genos_from_reps(row) ]
      genocounts = tally(geno)
      genocounts.pop( (None,None), None )
      key = tuple(imap(itemgetter(0), sorted(genocounts.iteritems(), key=itemgetter(1), reverse=True)))

      models  = modelmap.get(key)
      if models is None:
        model = UnphasedMarkerModel()
        for g in key:
          model.add_genotype(g)
        models = [model]*len(header)
        modelmap[key] = models

      descr = GenotypeArrayDescriptor(models)
      row = GenotypeArray(descr, geno)
      yield label,row

  elif format == 'sdat':
    if not isinstance(genos, (list,tuple)):
      genos = list(genos)

    header = genos[0]
    yield header

    counts = [ defaultdict(int) for i in range(len(header)) ]

    for label,row in islice(genos,1,None):
      geno = [ g or (None,None) for g in old_genorepr.genos_from_reps(row) ]
      for i,g in enumerate(geno):
        counts[i][g] += 1

    modelmap = {}
    models   = []
    for genocounts in counts:
      genocounts.pop( (None,None), None )
      key = tuple(imap(itemgetter(0), sorted(genocounts.iteritems(), key=itemgetter(1), reverse=True)))
      model = modelmap.get(key)

      if model is None:
        model = UnphasedMarkerModel()
        for g in key:
          model.add_genotype(g)
        modelmap[key] = model

      models.append(model)

    descr = GenotypeArrayDescriptor(models)

    for label,row in islice(genos,1,None):
      geno = [ g or (None,None) for g in old_genorepr.genos_from_reps(row) ]
      yield label,GenotypeArray(descr, geno)

  else:
    raise ValueError('Unsupported format')


def test(descr,filename,command,genotypes):
  t = time.time()
  command(filename)
  t = time.time()-t
  s = os.stat(filename).st_size
  gps,bpg = genotypes/t,8.0*s/genotypes
  print '%-38s  %6.2fs  %10d  %10d  %6.2f'  % (descr,t,s,gps,bpg)


def main():
  from   random import shuffle

  from biozilla.genodata import (load_genomatrixstream, save_genomatrix,
                                 build_genotriples_by_locus, recode_genomatrix,
                                 save_genotriples, load_genotriples)

  if 0:
    f      = '/usr/local/share/hapmap/build21/fwd_strand/non-redundant/genotypes_chr2_CEU_r21a_nr_fwd.txt.gz'
    #f     = '/usr/local/share/hapmap/build21/fwd_strand/non-redundant/genotypes_chr2_YRI_r21a_nr_fwd.txt.gz'
    matrix = load_genomatrixstream(f,'hapmap').materialize()
    format = 'hapmap'
  else:
    f = '/home/jacobske/projects/CGEMS/Scans/Breast/1/current/genotypes/STUDY/subjects_STUDY_CASE.ldat.gz'
    matrix = load_genomatrixstream(f).materialize()
    format = None

  bmatrix  = list(recode_genomatrix_bitpacked(matrix,  'ldat', snp_acgt))

  if 0:
    matrix2  = matrix.transposed().materialize()
    bmatrix2 = list(recode_genomatrix_bitpacked(matrix2, 'sdat', snp_acgt))
    matrix2  = list(matrix2)

  matrix   = list(matrix)

  loci     = len(matrix)-1
  subjects = len(matrix[0])
  g = loci*subjects
  print 'DATA: %d loci, %d subjects, %d genotypes' % (loci,subjects,g)
  print

  test('Load   compressed HapMap file', f,
         lambda f: ilen(load_genomatrixstream(f,format)), g)

  if 0:
    test('Save uncompressed   triple file ldat', 'data/g2.trip',
           lambda f: save_genotriples(f,build_genotriples_by_locus(matrix)), g)
    test('Save   compressed   triple file ldat', 'data/g2.trip.gz',
           lambda f: save_genotriples(f,build_genotriples_by_locus(matrix)), g)
    test('Save       binary   triple file ldat', 'data/g2.tdat',
           lambda f: save_genotriples_binary(f,build_genotriples_by_locus(bmatrix)), g)
    test('Load uncompressed   triple file ldat', 'data/g2.trip',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load   compressed   triple file ldat', 'data/g2.trip.gz',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load       binary   triple file ldat', 'data/g2.tdat',
           lambda f: ilen(load_genotriples_binary(f)), g)

  if 1:
    test('Save   compressed   ldat file',        'data/g2.ldat.gz',
           lambda f: save_genomatrix(f,matrix,format='ldat'), g)

  if 1:
    test('Save uncompressed   ldat file',        'data/g2.ldat',
           lambda f: save_genomatrix(f,matrix,format='ldat'), g)
    test('Load   compressed   ldat file',        'data/g2.ldat.gz',
           lambda f: ilen(load_genomatrixstream(f,'ldat')), g)
    test('Load uncompressed   ldat file',        'data/g2.ldat',
           lambda f: ilen(load_genomatrixstream(f,'ldat')), g)

  if 1:
    test('Save       binary   ldat file',        'data/g2.gdat',
           lambda f: save_genomatrix_binary(f,bmatrix,format='ldat'), g)

  if 1:
    test('Load       binary   ldat file',        'data/g2.gdat',
           lambda f: xlen(load_genomatrix_binary(f,'ldat')), g)

  test('Save      ubinary   ldat file',        'data/u2.gdat',
         lambda f: save_genomatrix_binary(f,bmatrix,format='ldat',compress=False), g)
  test('Load      ubinary   ldat file',        'data/u2.gdat',
         lambda f: xlen(load_genomatrix_binary(f,'ldat')), g)

  test('Load       binary   ldat file as sdat', 'data/g2.gdat',
         lambda f: xlen(load_genomatrix_binary(f,'sdat')), g)

  test('Load   compressed   ldat file as sdat',        'data/g2.ldat.gz',
         lambda f: ilen(load_genomatrixstream(f,'ldat').as_sdat()), g)

  matrix  = None
  bmatrix = None

  # Materialize for use later (but don't time)
  if 1:
    format,bmatrix2 = load_genomatrix_binary('data/g2.gdat','sdat')
    bmatrix2 = list(bmatrix2)
    matrix2 = [None]*len(bmatrix2)
    matrix2[0] = bmatrix2[0]
    for i in range(1,len(bmatrix2)):
      label,genos = bmatrix2[i]
      matrix2[i] = label,snp_acgt.pack_genos(map(Genotype.alleles,genos))

  test('Save       binary   sdat file',        'data/g22.gdat',
         lambda f: save_genomatrix_binary(f,bmatrix2,format='sdat'), g)
  test('Load       binary   sdat file',        'data/g22.gdat',
         lambda f: xlen(load_genomatrix_binary(f,'sdat')), g)
  test('Save      ubinary   sdat file',        'data/u22.gdat',
         lambda f: save_genomatrix_binary(f,bmatrix2,format='sdat',compress=False), g)
  test('Load      ubinary   sdat file',        'data/u22.gdat',
         lambda f: xlen(load_genomatrix_binary(f,'sdat')), g)

  test('Load       binary   sdat file as ldat', 'data/g22.gdat',
         lambda f: xlen(load_genomatrix_binary(f,'ldat')), g)

  test('Save   compressed   sdat file',        'data/g2.sdat.gz',
         lambda f: save_genomatrix(f,matrix2,format='sdat'), g)
  test('Save uncompressed   sdat file',        'data/g2.sdat',
         lambda f: save_genomatrix(f,matrix2,format='sdat'), g)
  test('Load   compressed   sdat file',        'data/g2.sdat.gz',
         lambda f: ilen(load_genomatrixstream(f,'sdat')), g)
  test('Load uncompressed   sdat file',        'data/g2.sdat',
         lambda f: ilen(load_genomatrixstream(f,'sdat')), g)

  test('Load   compressed   sdat file as ldat','data/g2.sdat.gz',
         lambda f: ilen(load_genomatrixstream(f,'sdat').as_ldat()), g)

  if 0:
    test('Save uncompressed   triple file sdat', 'data/g22.trip',
           lambda f: save_genotriples(f,build_genotriples_by_locus(matrix2)), g)
    test('Save   compressed   triple file sdat', 'data/g22.trip.gz',
           lambda f: save_genotriples(f,build_genotriples_by_locus(matrix2)), g)
    test('Save       binary   triple file sdat', 'data/g23.tdat',
           lambda f: save_genotriples_binary(f,build_genotriples_by_locus(bmatrix2)), g)
    test('Load uncompressed   triple file sdat', 'data/g22.trip',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load   compressed   triple file sdat', 'data/g22.trip.gz',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load       binary   triple file sdat', 'data/g23.tdat',
           lambda f: ilen(load_genotriples_binary(f)), g)

  if 0: # VERY VERY VERY VERY (VERY!) EXPENSIVE
    triples = list(build_genotriples_by_locus(matrix2))
    shuffle(triples)

    test('Save uncompressed   triple file random', 'data/g32.trip',
           lambda f: save_genotriples(f,triples), g)
    test('Save   compressed   triple file random', 'data/g32.trip.gz',
           lambda f: save_genotriples(f,triples), g)

    triples = list(build_genotriples_by_locus(bmatrix2))
    shuffle(triples)

    test('Save       binary   triple file random', 'data/g33.tdat',
           lambda f: save_genotriples_binary(f,triples), g)
    test('Load uncompressed   triple file random', 'data/g32.trip',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load   compressed   triple file random', 'data/g32.trip.gz',
           lambda f: ilen(load_genotriples(f)), g)
    test('Load       binary   triple file random', 'data/g33.tdat',
           lambda f: ilen(load_genotriples_binary(f)), g)


if __name__ == '__main__':
  main()
