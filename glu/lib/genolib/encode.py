# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GLU genotype data encoding functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   operator                  import getitem
from   itertools                 import izip,imap,repeat,count

from   glu.lib.utils             import is_str

from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,            \
                                        GenotypeLookupError, GenotypeRepresentationError, \
                                        build_model, build_descr


def _sample_encoding_error(loci,models,genos,warn=False):
  '''
  Handle sample genotype encoding errors by either producing an informative
  exception or a warning message.
  '''
  for i in xrange(len(genos)):
    model = models[i]
    geno  = genos[i]
    try:
      model[geno]
    except GenotypeLookupError:
      _encoding_error(loci[i],geno,model,warn)
      genos[i] = model[None,None]


def _encoding_error(locus,item,model,warn=False):
  '''
  Handle genotype encoding error by either producing an informative exception
  or a warning message.
  '''
  if is_str(item) or len(item)==1:
    item = 'allele %s' % (''.join(item))
  else:
    item = 'genotype %s' % (','.join(item))

  msg = 'Locus model %s cannot accommodate %s (max_alleles=%d,alleles=%s)' \
                      % (locus,item,model.max_alleles,','.join(model.alleles[1:]))

  if warn:
    sys.stderr.write('[WARNING] %s\n' % msg)
  else:
    raise GenotypeRepresentationError(msg)


def merge_locus(loc,new_model):
  old_model = loc.model

  if old_model is None:
    loc.model = new_model
  elif new_model is None or new_model.replaceable_by(old_model):
    pass
  elif new_model is old_model or old_model.replaceable_by(new_model):
    loc.model = new_model
  else:
    # Recoding needed
    loc.model = build_model(alleles=new_model.alleles[1:],base=loc.model)

    return True

  return False


def merge_locus2(old_model,new_model):
  if old_model is None:
    return False,new_model
  elif new_model is None or new_model.replaceable_by(old_model):
    return False,old_model
  elif new_model is old_model or old_model.replaceable_by(new_model):
    return False,new_model

  # Recoding needed
  return True,build_model(alleles=new_model.alleles[1:],base=old_model)


def update_model(old_model, new_model):
  if old_model is new_model or old_model.replaceable_by(new_model):
    return new_model

  # i.e., an unnecessary model update
  if not new_model.replaceable_by(old_model):
    raise GenotypeRepresentationError('Invalid model update')

  # Return the old model
  return old_model


def pack_genomatrixstream(genos):
  '''
  Transform a genomatrix into an internal packed representation

  @param      genos: genomatrix stream
  @type       genos: sequence
  @param     genome: genome descriptor
  @type      genome: Genome instance
  @return          : genomatrix with a packed internal format
  @rtype           : genomatrix generator

  >>> from glu.lib.genolib.streams import GenomatrixStream
  >>> samples = ('s1','s2','s3')
  >>> rows = [('l1',[ ('A', 'A'), (None, None), ('G', 'G') ]),
  ...          ('l2',[(None, None),(None, None),(None, None)]),
  ...          ('l3',[ ('A', 'A'), (None, None),(None, None)]),
  ...          ('l4',[ ('G', 'T'), (None, None), ('T', 'T') ])]
  >>> genos = GenomatrixStream.from_tuples(rows,'ldat',samples=samples)
  >>> genos = pack_genomatrixstream(genos)
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> for row in genos:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])
  '''
  if not genos.columns:
    return genos.clone([],packed=True,samples=(),loci=(),models=[])

  if genos.format=='sdat':
    updates = []

    def _pack(genos):
      locs   = [ genos.genome.get_locus_model(lname) for lname in genos.loci ]
      models = []

      for loc,model in izip(locs,genos.models):
        new_model = loc.model = update_model(loc.model,model)
        models.append(new_model)

      descr  = GenotypeArrayDescriptor(models)

      for label,row in genos:
        if genos.updates:
          for i,model in genos.updates:
            loc       = locs[i]
            new_model = update_model(loc.model,model)
            new_model = update_model(models[i],new_model)
            if new_model is not loc.model or new_model is not models[i]:
              loc.model = descr[i] = models[i] = new_model
              updates.append( (i,new_model) )

          genos.updates[:] = []

        yield label,GenotypeArray(descr,row)

  else:
    updates = genos.updates

    def _pack(genos):
      n = len(genos.columns)

      for lname,row in genos:
        model = genos.genome.get_model(lname)
        descr = build_descr(model,n)
        yield lname,GenotypeArray(descr,row)


  return genos.clone(_pack(genos),packed=True,materialized=False,updates=updates)


def recode_genomatrixstream(genos, genome, warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded with representations
  defined by the supplied genome object.  Locus metadata other than models
  are merged and discrepancies raise errors, If genotype models change, then
  all genotypes are recoded to use the same representation provided the
  models are compatible.

  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : new genomatrixstream with encoding identical to the
                       supplied genome
  @rtype             : GenomatrixStream

  >>> from glu.lib.genolib.streams   import GenomatrixStream
  >>> from glu.lib.genolib.genoarray import Genotype
  >>> defmodel = build_model(alleles='ACGT',allow_hemizygote=True)
  >>> genome = Genome()
  >>> genome.set_locus('l1',defmodel)
  >>> genome.set_locus('l2',defmodel)
  >>> genome.set_locus('l3',defmodel)
  >>> genome.set_locus('l4',defmodel)

  Test ldat "unknown" models remapped to a default model

  >>> samples = ('s1', 's2', 's3')
  >>> rows = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...         ('l2', [(None, None), (None, None),  (None, None)]),
  ...         ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...         ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples)
  >>> len(genos.models)
  0
  >>> genos = recode_genomatrixstream(genos,genome).materialize()
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> all(model is defmodel for model in genos.models)
  True
  >>> for (label,row),model in izip(genos,genos.models):
  ...   print label,row,model is defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test sdat known models

  >>> samples = ('s1', 's2', 's3')
  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).as_sdat()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,genome).materialize()
  >>> genos.loci
  ('l1', 'l2', 'l3', 'l4')
  >>> genos.samples
  ('s1', 's2', 's3')
  >>> all(model is defmodel for model in genos.models)
  True
  >>> for label,row in genos:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  Test ldat fastpath for no recoding

  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).materialize()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,Genome())
  >>> for (label,row),model in izip(genos,genos.models):
  ...   print label,row,model is not defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test ldat fastpath for known models

  >>> genos = GenomatrixStream.from_tuples(rows, 'ldat', samples=samples).materialize()
  >>> len(genos.models)
  4
  >>> genos = recode_genomatrixstream(genos,genome)
  >>> for (label,row),model in izip(genos,genos.models):
  ...   print label,row,model is defmodel,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True True
  l2 [(None, None), (None, None), (None, None)] True True
  l3 [('A', 'A'), (None, None), (None, None)] True True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True True

  Test sequential recoding to ensure that the same model is used throughout

  >>> samples =          ('s1',         's2',        's3')
  >>> rows1 = [('l1',[ ('G', 'G'),   ('G', 'T'),  ('T', 'T')]),
  ...          ('l2',[ ('A', 'A'),   ('T', 'T'),  ('A', 'T')])]
  >>> rows2 = [('l1',[(None, None),  ('T', 'T'),  ('G', 'G')]),
  ...          ('l3',[ ('A', 'A'),  (None, None), ('A', 'T')])]
  >>> genos1 = GenomatrixStream.from_tuples(rows1,'ldat',samples=samples)
  >>> genos2 = GenomatrixStream.from_tuples(rows2,'ldat',samples=samples)
  >>> genome = Genome()
  >>> genos1 = recode_genomatrixstream(genos1, genome).materialize()
  >>> sorted(genome.loci)
  ['l1', 'l2']
  >>> for locus,row in genos1:
  ...   print locus,row
  l1 [('G', 'G'), ('G', 'T'), ('T', 'T')]
  l2 [('A', 'A'), ('T', 'T'), ('A', 'T')]

  >>> genos2 = recode_genomatrixstream(genos2, genome).materialize()
  >>> for locus,row in genos2:
  ...   print locus,row
  l1 [(None, None), ('T', 'T'), ('G', 'G')]
  l3 [('A', 'A'), (None, None), ('A', 'T')]

  >>> sorted(genome.loci)
  ['l1', 'l2', 'l3']
  >>> for locus,model in genos1.model_pairs:
  ...   assert genome.get_model(locus) is model
  >>> for locus,model in genos2.model_pairs:
  ...   assert genome.get_model(locus) is model
  '''
  # Fastpath for null recoding -- DISABLED due to some operations leaving
  # streams with inconsistent encoding (like renaming ldat rows)
  #if genos.genome is genome:
  #  return genos

  # Data for slowpath
  models  = []
  updates = []

  if genos.format=='ldat':
    def _recode_genomatrixstream():
      n = len(genos.samples)

      for (lname,row),old_model in izip(genos,genos.models):
        old_locus = genos.genome.loci[lname]

        if lname not in genome.loci:
          loc = genome.loci[lname] = old_locus
        else:
          genome.merge_locus(lname, None, old_locus.chromosome, old_locus.location,
                                          old_locus.strand, warn)
          loc = genome.loci[lname]

        try:
          recode,loc.model = merge_locus2(loc.model, old_model)
        except GenotypeRepresentationError:
          # FIXME: The third argument may be incorrect elsewhere too
          _encoding_error(lname,set(old_model.alleles)-set(loc.model.alleles),loc.model,warn)
          row = [loc.model.genotypes[0]]*len(row)
          packed = False

        model = loc.model

        # If recoding or packing is required
        if recode or not genos.packed:
          descr = build_descr(model,n)
          row   = GenotypeArray(descr,row)

        # N.B. Does not aggressively recode model to loc.model unless
        # encoding changed or stream is not packed

        models.append(model)
        yield lname,row

  # sdat format
  elif genos.format=='sdat':
    assert genos.loci is not None and len(genos.loci) == len(genos.models)

    for i,lname in enumerate(genos.loci):
      old_locus = genos.genome.loci[lname]
      old_model = old_locus.model or genos.models[i]

      # Fastpath since no merging is necessary-- just fire and forget
      if lname not in genome.loci:
        loc = genome.loci[lname] = old_locus
      else:
        genome.merge_locus(lname, None, old_locus.chromosome, old_locus.location,
                                        old_locus.strand, warn)

        loc = genome.get_locus(lname)

        try:
          recode,loc.model = merge_locus2(loc.model, old_model)
        except GenotypeRepresentationError:
          _encoding_error(lname,set(old_model.alleles)-set(loc.model.alleles),loc.model,warn)

      models.append(loc.model)

    # No FASTPATH is generally possible, since other streams may be updating our genome

    def _recode_genomatrixstream():
      descr = GenotypeArrayDescriptor(models)

      for sample,row in genos:
        if genos.updates:
          for i,model in genos.updates:
            loc = genome.get_locus(genos.loci[i])
            old_model = loc.model

            try:
              recode,loc.model = merge_locus2(loc.model, model)
            except GenotypeRepresentationError:
              _encoding_error(genos.loci[i],set(old_model.alleles)-set(model.alleles),model,warn)

            if loc.model is not models[i]:
              descr[i] = models[i] = loc.model
              updates.append( (i,loc.model) )

          genos.updates[:] = []

        row = GenotypeArray(descr,row)

        yield sample,row

  else:
    raise ValueError('Unknown format')

  return genos.clone(_recode_genomatrixstream(),models=models,updates=updates,
                     genome=genome,packed=True,materialized=False)


def encode_genomatrixstream_from_tuples(columns, genos, format, genome=None,
                                                 unique=False, warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       format: format of input genomatrix, either 'ldat' or 'sdat'
  @type        format: str
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
  @type        unique: bool
  @return            : tuple of columns and a genomatrix generator in packed format
  @rtype             : 2-tuple of list of str and genomatrix generator

  >>> from glu.lib.genolib.genoarray import Genotype
  >>> defmodel = build_model(alleles='ACGT',allow_hemizygote=True)
  >>> genome  = Genome(default_model=defmodel, default_fixed=True)

  With genome:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
  >>> samples
  ('s1', 's2', 's3')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True
  l2 [(None, None), (None, None), (None, None)] True
  l3 [('A', 'A'), (None, None), (None, None)] True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')]),
  ...          ('s2', [(None, None), (None, None), (None, None), (None, None)]),
  ...          ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])]
  >>> loci,models,updates,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  No genome:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat')
  >>> samples
  ('s1', 's2', 's3')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A'), (None, None), ('G', 'G')] True
  l2 [(None, None), (None, None), (None, None)] True
  l3 [('A', 'A'), (None, None), (None, None)] True
  l4 [('G', 'T'), (None, None), ('T', 'T')] True

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')]),
  ...          ('s2', [(None, None), (None, None), (None, None), (None, None)]),
  ...          ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])]
  >>> loci,models,updates,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  s1 [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')] True
  s2 [(None, None), (None, None), (None, None), (None, None)] True
  s3 [('G', 'G'), (None, None), (None, None), ('T', 'T')] True

  See if we can provide a subtle cache bug when models are cached too
  aggressively for non-unique loci:

  >>> samples = ('s1',)
  >>> genos = [('l1', [ ('A','A') ]),
  ...          ('l2', [ ('A','A') ]),
  ...          ('l1', [ ('A','T') ]),
  ...          ('l2', [ ('A','G') ])]
  >>> genome = Genome()
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
  >>> new_rows = list(new_rows)
  >>> samples
  ('s1',)
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A')] True
  l2 [('A', 'A')] True
  l1 [('A', 'T')] True
  l2 [('A', 'G')] True
  '''
  if genome is None:
    genome = Genome()

  if not columns:
    return columns,[],[],genome,[]

  models  = []
  updates = []

  if format=='ldat':
    def _encode():
      n = len(columns)

      def get_alleles(genos):
        alleles = set(a for g in set(genos) for a in g)
        alleles.discard(None)
        return tuple(sorted(alleles))

      for lname,row in genos:
        loc = genome.get_locus(lname)

        if loc.model is None:
          loc.model = build_model(alleles=get_alleles(row))

        model = loc.model
        descr = build_descr(model,n)

        try:
          row = GenotypeArray(descr,row)
        except GenotypeLookupError:
          try:
            new_alleles = get_alleles(row)
            new_model = build_model(alleles=new_alleles,base=model)
          except GenotypeRepresentationError:
            _encoding_error(lname,set(new_alleles)-set(model.alleles),model,warn)

          loc.model = model = new_model

          descr = build_descr(model,n)
          row = GenotypeArray(descr,row)

        models.append(model)
        yield lname,row

  elif format=='sdat':
    loci = []
    for i,lname in enumerate(columns):
      loc = genome.get_locus_model(lname)
      loci.append(loc)
      models.append(loc.model)

    def _encode():
      descr = GenotypeArrayDescriptor(models)

      for sample,row in genos:
        try:
          row = GenotypeArray(descr,row)
        except GenotypeLookupError:
          new_row = [None]*len(columns)
          for i,geno,model in izip(count(),row,models):
            if geno not in model:
              try:
                model = build_model(alleles=geno,base=model)
              except GenotypeRepresentationError:
                _encoding_error(columns[i],set(geno)-set(model.alleles),model,warn)

              loci[i].model = descr[i] = models[i] = model
              updates.append( (i,model) )

            new_row[i] = model[geno]

          row = GenotypeArray(descr,new_row)

        yield sample,row
  else:
    raise ValueError('Unknown format')

  return columns,models,updates,genome,_encode()


def encode_genomatrixstream_from_strings(columns,genos,format,genorepr,genome=None,
                                                 unique=False,warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded to a new internal representation

  @param      columns: matrix column names
  @type       columns: sequence of strs
  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       format: format of input genomatrix, either 'ldat' or 'sdat'
  @type        format: str
  @param     genorepr: internal representation of genotypes for the input/output
  @type      genorepr: UnphasedMarkerRepresentation or similar object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: flag indicating if repeated elements do not exist within the stream. Default is 'False'
  @type        unique: bool
  @return            : tuple of columns and a genomatrix generator in packed format
  @rtype             : 2-tuple of list of str and genomatrix generator

  >>> from glu.lib.genolib.reprs import snp
  >>> from glu.lib.genolib.genoarray import Genotype
  >>> defmodel  = build_model(alleles='ACGT',allow_hemizygote=True)
  >>> genome = Genome(default_model=defmodel, default_fixed=True)

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', ['AA', '  ', 'AA', 'GT']),
  ...          ('s2', ['  ', '',   '  ', '  ']),
  ...          ('s3', ['GG', '  ', '  ', 'TT'])]
  >>> loci,models,updates,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp,genome)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in new_rows:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp)
  >>> samples
  ('s1', 's2', 's3')
  >>> for row in new_rows:
  ...   print row
  ('l1', [('A', 'A'), (None, None), ('G', 'G')])
  ('l2', [(None, None), (None, None), (None, None)])
  ('l3', [('A', 'A'), (None, None), (None, None)])
  ('l4', [('G', 'T'), (None, None), ('T', 'T')])

  >>> loci = ('l1','l2','l3', 'l4')
  >>> genos = [('s1', ['AA', '  ', 'AA', 'GT']),
  ...          ('s2', ['  ', '',   '  ', '  ']),
  ...          ('s3', ['GG', '  ', '  ', 'TT'])]
  >>> loci,models,updates,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp)
  >>> loci
  ('l1', 'l2', 'l3', 'l4')
  >>> for row in new_rows:
  ...   print row
  ('s1', [('A', 'A'), (None, None), ('A', 'A'), ('G', 'T')])
  ('s2', [(None, None), (None, None), (None, None), (None, None)])
  ('s3', [('G', 'G'), (None, None), (None, None), ('T', 'T')])

  See if we can provide a subtle cache bug when models are cached too
  aggressively for non-unique loci:

  >>> samples = ('s1',)
  >>> genos = [('l1', [ ('AA') ]),
  ...          ('l2', [ ('AA') ]),
  ...          ('l1', [ ('AT') ]),
  ...          ('l2', [ ('AG') ])]
  >>> genome = Genome()
  >>> samples,models,updates,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
  >>> new_rows = list(new_rows)
  >>> samples
  ('s1',)
  >>> for label,row in new_rows:
  ...   print label,row,all(isinstance(g,Genotype) for g in row)
  l1 [('A', 'A')] True
  l2 [('A', 'A')] True
  l1 [('A', 'T')] True
  l2 [('A', 'G')] True
  '''
  if genome is None:
    genome = Genome()

  models  = []
  updates = []

  if format=='ldat':
    def _encode():
      n = len(columns)

      cachemap     = {}
      from_strings = genorepr.from_strings

      for lname,row in genos:
        loc   = genome.get_locus(lname)
        cache = cachemap.get(loc.model)

        # FAST PATH: known model, existing string cache
        if loc.model is not None and cache is not None:
          try:
            descr = build_descr(loc.model,n)
            row   = GenotypeArray(descr,imap(getitem, repeat(cache), row))

            models.append(loc.model)
            yield lname,row
            continue

          except (GenotypeLookupError,KeyError):
            pass

        # SLOW PATH: no model, new alleles, no string cache, or new representation
        gstrs   = tuple(set(row))
        gtups   = from_strings(gstrs)
        alleles = set(a for g in gtups for a in g)
        alleles.discard(None)

        try:
          loc.model = build_model(alleles,base=loc.model)
        except GenotypeRepresentationError:
          _encoding_error(lname,set(alleles)-set(loc.model.alleles),loc.model,warn)
          # FIXME: Need recovery code
          continue

        cache = cachemap.get(loc.model)
        if cache is None:
          cache = cachemap[loc.model] = {}
        cache.update( (gs,loc.model[gt]) for gs,gt in izip(gstrs,gtups) )

        descr = build_descr(loc.model,n)
        row   = GenotypeArray(descr,imap(getitem, repeat(cache), row))

        models.append(loc.model)
        yield lname,row

  elif format=='sdat':
    loci      = []
    cachemap  = {}

    to_string    = genorepr.to_string

    def mk_cache(model):
      cache = dict( (to_string(g),g) for g in model.genotypes )
      cache.update( (genorepr.to_string_from_alleles((g[1],g[0])),g) for g in model.genotypes )
      for g in genorepr.missing_geno_strs:
        cache[g] = model[None,None]
      return cache

    for lname in columns:
      loc = genome.get_locus_model(lname)
      loci.append(loc)

      model = loc.model
      models.append(model)

    def _encode():
      repr   = genorepr.from_string
      descr  = GenotypeArrayDescriptor(models)
      errors = (KeyError,ValueError,GenotypeRepresentationError)

      cachelist = [{}]*len(columns)

      for sample,row in genos:
        try:
          row = GenotypeArray(descr,imap(getitem, cachelist, row) )
        except errors:
          new_row = [None]*len(loci)

          for i,gstr,cache in izip(count(),row,cachelist):
            g = cache.get(gstr)
            if g is None:
              loc  = loci[i]
              geno = repr(gstr)

              if geno not in loc.model:
                try:
                  loc.model = build_model(alleles=geno,base=loc.model)
                except GenotypeRepresentationError:
                  _encoding_error(columns[i],set(geno)-set(loc.model.alleles),loc.model,warn)
                  # FIXME: Need recovery code

              model = loc.model
              if models[i] is not model:
                # Update model in various locations
                descr[i] = models[i] = model
                updates.append( (i,model) )

              # Update string to genotype cache
              cache = cachemap.get(model)
              if cache is None:
                cache = cachemap[model] = mk_cache(model)
              cachelist[i] = cache

              # Update cache
              g = loc.model[geno]
              cache[gstr] = g

              assert g is loc.model[geno] and cache[gstr] is g

            # Set genotype
            new_row[i] = g

          if 0: # DEBUG
            for i,g in enumerate(new_row):
              assert g in descr[i], 'Genotype %s not in model %s at locus %s' % (type(g),descr[i].genotypes[1:],columns[i])
              assert g.model.replaceable_by(models[i]),'model with alleles "%s" is not replaceable by model with alleles "%s"' % (''.join(g.model.alleles[1:]),''.join(models[i].alleles[1:]))

          row = GenotypeArray(descr,new_row)

        yield sample,row

  else:
    raise ValueError('Unknown format')

  return columns,models,updates,genome,_encode()


def recode_genotriples(triples,genome,warn=False):
  '''
  Returns a new genotriples with the genotypes encoded to a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation). e.g.
                       ('s1','l1','AA'),...
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> from glu.lib.genolib.streams import GenotripleStream
  >>> triples = [('s3', 'l1', ('G', 'G')),('s3', 'l2', ('A', 'A')),
  ...            ('s2', 'l3', ('G', 'T')),('s1', 'l1', ('T', 'T')),
  ...            ('s1', 'l1', ('G', 'G')),('s2', 'l2', ('A', 'A'))]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> genome = Genome()
  >>> for row in recode_genotriples(triples,genome):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  >>> sorted(genome.loci)
  ['l1', 'l2', 'l3']
  '''
  updates = []
  def _recode():
    updated = {}
    for sample,lname,geno in triples:
      model = updated.get(lname)

      if model is None:
        old_model = geno.model
        old_locus = triples.genome.loci[lname]
        #assert old_model is triples.genome.loci[lname].model
        #assert old_locus.model is old_model

        genome.merge_locus(lname, chromosome=old_locus.chromosome, location=old_locus.location)

        loc = genome.get_locus(lname)
        if loc.model is None:
          loc.model = model = geno.model
        elif geno.model is not loc.model:
          try:
            loc.model = model = build_model(alleles=old_model.alleles[1:],base=loc.model)
            updates.append( (lname,model) )
          except GenotypeRepresentationError:
            _encoding_error(lname,set(old_model.alleles)-set(loc.model.alleles),loc.model,warn)
        else:
          model = loc.model

        updated[lname] = model

      assert model is not None

      try:
        geno = model[geno]
      except GenotypeLookupError:
        loc = genome.get_locus(lname)
        try:
          loc.model = model = build_model(alleles=geno.model.alleles[1:],base=loc.model)
          updates.append( (lname,model) )
        except GenotypeRepresentationError:
          _encoding_error(lname,set(geno.model.alleles)-set(model.alleles),model,warn)

        updated[lname] = model
        geno = model[geno]

      yield sample,lname,geno

  return triples.clone(_recode(),updates=updates,genome=genome,materialized=False)


def encode_genotriples_from_tuples(triples,genome,warn=False):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation)
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> triples = [('s3', 'l1', ('G', 'G')),('s3', 'l2', ('A', 'A')),
  ...            ('s2', 'l3', ('G', 'T')),('s1', 'l1', ('T', 'T')),
  ...            ('s1', 'l1', ('G', 'G')),('s2', 'l2', ('A', 'A'))]
  >>> triples,updates = encode_genotriples_from_tuples(triples,Genome())
  >>> for row in triples:
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  updates = []

  def _encode():
    for sample,lname,geno in triples:
      loc = genome.get_locus_model(lname)

      try:
        geno = loc.model[geno]
      except GenotypeLookupError:
        try:
          new_model = loc.model = build_model(alleles=geno,base=loc.model)
        except GenotypeRepresentationError:
          _encoding_error(lname,geno,loc.model,warn)

        updates.append( (lname,new_model) )

        geno = new_model[geno]

      yield sample,lname,geno

  return _encode(),updates


def encode_genotriples_from_strings(triples,genorepr,genome,warn=False):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation)
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> from glu.lib.genolib.reprs import snp
  >>> triples = [('s3', 'l1', 'GG'),('s3', 'l2', 'AA'),
  ...            ('s2', 'l3', 'GT'),('s1', 'l1', 'TT'),
  ...            ('s1', 'l1', 'GG'),('s2', 'l2', 'AA')]
  >>> triples,updates = encode_genotriples_from_strings(triples,snp,Genome())
  >>> for row in triples:
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  updates = []
  def _encode():
    from_string = genorepr.from_string
    to_string   = genorepr.to_string

    def mk_cache(model):
      cache = dict( (to_string(g),g) for g in model.genotypes )
      cache.update( (genorepr.to_string_from_alleles((g[1],g[0])),g) for g in model.genotypes )
      for g in genorepr.missing_geno_strs:
        cache[g] = model[None,None]
      return cache

    cachemap = {}
    for sample,lname,geno in triples:
      loc = genome.get_locus_model(lname)
      model = loc.model

      cache = cachemap.get(model)
      if cache is None:
        cache = cachemap[model] = mk_cache(model)

      try:
        geno = cache[geno]
      except KeyError:
        try:
          loc.model = model = build_model(alleles=geno,base=loc.model)
        except GenotypeRepresentationError:
          _encoding_error(lname,geno,loc.model,warn)

        updates.append( (lname,model) )

        cache = cachemap.get(model)
        if cache is None:
          cache = cachemap[model] = mk_cache(model)

        geno = cache[geno] = model[from_string(geno)]

      yield sample,lname,geno

  return _encode(),updates


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
