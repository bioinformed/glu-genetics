# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GLU genotype data encoding functions'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   operator                  import getitem
from   itertools                 import izip,imap,chain,repeat

from   glu.lib.utils             import izip_exact,is_str

from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.genoarray import GenotypeArrayDescriptor,GenotypeArray,Genotype,   \
                                        GenotypeLookupError, GenotypeRepresentationError, \
                                        model_from_alleles


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
  if is_str(item):
    item = 'allele %s' % item
  else:
    item = 'genotype %s' % (','.join(item))

  msg = 'Locus model %s cannot accommodate %s (max_alleles=%d,alleles=%s)' \
                      % (locus,item,model.max_alleles,','.join(model.alleles[1:]))

  if warn:
    sys.stderr.write('[WARNING] %s\n' % msg)
  else:
    raise GenotypeRepresentationError(msg)


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
  def _pack(genos):
    if genos.format=='sdat':
      descr = GenotypeArrayDescriptor(genos.models)
      for label,row in genos:
        yield label,GenotypeArray(descr,row)
    else:
      n = len(genos.columns)
      descrcache = {}

      for label,row in genos:
        model = row[0].model
        descr = descrcache.get(model)
        if descr is None:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        yield label,GenotypeArray(descr,row)

  if not genos.columns:
    return genos

  return genos.clone(_pack(genos),packed=True)


# FIXME: Implement COW on genotype models
def recode_genomatrixstream(genos, genome, warn=False):
  '''
  Returns a new genomatrix with the genotypes encoded with representations
  defined by the supplied genome object.  Locus metadata other than models
  are merged and discrepencies raise errors, If genotype models change, then
  all genotypes are recoded to use the same representation provided the
  models are compatible.

  @param        genos: genomatrix stream
  @type         genos: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : new genomatrixstream with encoding identical to the
                       supplied genome
  @rtype             : GenomatrixStream

  >>> from glu.lib.genolib.streams import GenomatrixStream
  >>> defmodel = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome = Genome()
  >>> genome.set_locus('l1',model=defmodel)
  >>> genome.set_locus('l2',model=defmodel)
  >>> genome.set_locus('l3',model=defmodel)
  >>> genome.set_locus('l4',model=defmodel)

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
  >>> for (label,row),model in izip_exact(genos,genos.models):
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
  >>> for (label,row),model in izip_exact(genos,genos.models):
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
  >>> for (label,row),model in izip_exact(genos,genos.models):
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
  # streams with inconsistant encoding (like renaming ldat rows)
  #if genos.genome is genome:
  #  return genos

  # Data for slowpath
  models = []

  # All loci and models are known
  if genos.loci is not None and len(genos.models) == len(genos.loci):
    recode = False
    for lname,old_model in izip(genos.loci,genos.models):
      # Get the new model or fix the old model
      old_locus = genos.genome.loci[lname]
      #assert old_locus.model is old_model or None in (old_model,old_locus.model)

      if lname not in genome.loci:
        loc = genome.loci[lname] = old_locus
      else:
        genome.merge_locus(lname, None, old_locus.fixed,    old_locus.chromosome,
                                        old_locus.location, old_locus.strand, warn)

        loc = genome.get_locus(lname)

      if loc.model is None:
        loc.model = old_model

      model = loc.model

      # Check to see if a full recoding or update is necessary
      if model is not old_model:
        recode = True
        try:
          for g in old_model.genotypes[1:]:
            model.add_genotype(g)
        except GenotypeRepresentationError:
          _encoding_error(lname,g,model,warn)

      models.append(model)

    # FASTPATH: No models change, so return with the updated genome
    if not recode:
      assert genos.models == models
      return genos.clone(genos.use_stream(), genome=genome, materialized=False)

  # MEDIUM PATH: Ldat, recoding needed, known models
  if models and genos.format=='ldat':
    def _recode_genomatrixstream():
      n = len(genos.samples)

      descrcache = {}
      packed = genos.packed

      for (lname,row),old_model,model in izip(genos,genos.models,models):
        # Cache the descriptor for this model, since we're likely to see it again
        if packed:
          assert old_model is row.descriptor.models[0]
          descr = descrcache[old_model] = row.descriptor

        # Get or build the new descriptor
        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        # If the model changed, recode by adding all genotypes and packing
        if old_model is not model:
          loc = genome.get_locus(lname)
          # FIXME: The semantics of the fixed flag are broken
          if not loc.fixed:
            # Unpack to speed updates and repacking
            row = row[:]
            try:
              for g in set(row):
                model.add_genotype(g)
            except GenotypeRepresentationError:
              _encoding_error(lname,g,model,warn)
              row = None

          row = GenotypeArray(descr,row)

        # Otherwise, only repack if necessary
        elif not packed:
          row = GenotypeArray(descr,row)

        yield lname,row

  # SLOWPATH: Ldat without model information
  elif genos.format=='ldat':
    def _recode_genomatrixstream():
      n = len(genos.samples)

      descrcache = {}
      packed = genos.packed

      for (lname,row),old_model in izip_exact(genos,genos.models):
        # Get the new model or fix the old model
        old_locus = genos.genome.loci[lname]
        #assert old_locus.model is old_model or None in (old_model,old_locus.model)

        if lname not in genome.loci:
          loc = genome.loci[lname] = old_locus
        else:
          genome.merge_locus(lname, None, old_locus.fixed,    old_locus.chromosome,
                                          old_locus.location, old_locus.strand, warn)

          loc = genome.get_locus(lname)

        if loc.model is None:
          loc.model = old_model
        else:
          model = loc.model

        # Cache the descriptor for this model, since we're likely to see it again
        if packed:
          assert old_model is row.descriptor.models[0]
          descr = descrcache[old_model] = row.descriptor

        # Get or build the new descriptor
        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        # If the model changed, recode by adding all genotypes and packing
        if old_model is not model:
          # FIXME: The semantics of the fixed flag are broken
          if not loc.fixed:
            # Unpack to speed updates and repacking
            row = row[:]
            try:
              for g in set(row):
                model.add_genotype(g)
            except GenotypeRepresentationError:
              _encoding_error(lname,g,model,warn)
              row = None

          row = GenotypeArray(descr,row)

        # Otherwise, only repack if necessary
        elif not packed:
          row = GenotypeArray(descr,row)

        models.append(model)
        yield lname,row

  # sdat format
  elif genos.format=='sdat':
    #assert genome is not genos.genome
    assert genos.loci is not None and len(genos.loci) == len(genos.models) == len(models)

    # Find all models that must be updated
    updates = []
    for i,(lname,old_model,model) in enumerate(izip(genos.loci,genos.models,models)):
      if model is not old_model:
        loc = genome.get_locus(lname)
        # FIXME: The semantics of the fixed flag are broken
        if not loc.fixed:
          updates.append( (i,model.add_genotype) )

    # FASTERPATH: If all models are fixed, recoding is straightforward
    if not updates:
      def _recode_genomatrixstream():
        descr = GenotypeArrayDescriptor(models)
        for sample,row in genos:
          try:
            # Recode and yield new row
            row = GenotypeArray(descr,row)
          except GenotypeLookupError:
            _sample_encoding_error(genos.loci,models,row,warn)

          yield sample,row

    # SLOWPATH: Otherwise, recode by adding genotypes from all changed
    # models and packing.
    else:
      def _recode_genomatrixstream():
        descr = GenotypeArrayDescriptor(models)
        for sample,row in genos:
          # Unpack row to speed access -- both updates and GenotypeArray will
          # need to unpack
          if updates:
            row = row[:]

          # Try to yield updated genotype array, hoping that all alleles are represented
          try:
            yield sample,GenotypeArray(descr,row)

          except GenotypeRepresentationError:
            # Update all changed models to ensure they contain the needed alleles
            for i,add in updates:
              try:
                add(row[i])
              except GenotypeRepresentationError:
                _encoding_error(genos.loci[i],row[i],models[i],warn)
                row[i] = models[i][None,None]

            # Recode and yield new row
            yield sample,GenotypeArray(descr,row)

  else:
    raise ValueError('Uknown format')

  return genos.clone(_recode_genomatrixstream(),models=models,genome=genome,packed=True,materialized=False)


# FIXME: Implement COW on genotype models
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

  >>> defmodel = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome  = Genome(default_model=defmodel, default_fixed=True)

  With genome:

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', [ ('A', 'A'),  (None, None),   ('G', 'G') ]),
  ...          ('l2', [(None, None), (None, None),  (None, None)]),
  ...          ('l3', [ ('A', 'A'),  (None, None),  (None, None)]),
  ...          ('l4', [ ('G', 'T'),  (None, None),   ('T', 'T')] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
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
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
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
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat')
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
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_tuples(loci,genos,'sdat')
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
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_tuples(samples,genos,'ldat',genome)
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

  models = []

  if format=='ldat':
    def _encode():
      n = len(columns)
      m = genome.max_alleles+1

      # FIXME: add support for hemizygote models
      modelcache = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                               if locus.model is not None and len(locus.model.alleles)==m )
      descrcache = {}

      def _genokey(genos):
        gset = set(a for g in set(genos) for a in g)
        gset.discard(None)
        return tuple(sorted(gset))

      for lname,row in genos:
        key = None

        loc = genome.get_locus(lname)

        if loc.model is None:
          key = _genokey(row)
          cachable = genome.default_model is None and len(key) == genome.max_alleles

          if cachable:
            # Aggressively reuse models with fully compatible alleles
            loc.model = modelcache.get(key)

          if loc.model is None:
            genome.get_locus_model(lname)

          if cachable:
            modelcache[key] = loc.model

        model = loc.model

        assert model is not None

        # FIXME: The semantics of the fixed flag are broken
        if not loc.fixed:
          key = key or _genokey(row)
          try:
            for a in key:
              model.add_allele(a)
          except GenotypeRepresentationError:
            _encoding_error(lname,a,model,warn)
            row = None

        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        assert descr.models[0] is model
        models.append(model)
        yield lname,GenotypeArray(descr,row)

  elif format=='sdat':

    genos = sdat_model_lookahead_from_tuples(columns,genos,genome)

    updates = []

    for i,lname in enumerate(columns):
      loc = genome.get_locus_model(lname)
      models.append(loc.model)
      # FIXME: The semantics of the fixed flag are broken
      if not loc.fixed:
        updates.append( (i,loc.model.add_genotype) )

    def _encode():
      descr = GenotypeArrayDescriptor(models)

      if not updates:
        for sample,row in genos:
          yield sample,GenotypeArray(descr,row)
      else:
        for sample,row in genos:
          for i,add in updates:
            try:
              add(row[i])
            except GenotypeRepresentationError:
              _encoding_error(columns[i],row[i],models[i],warn)
              row[i] = models[i][None,None]

          yield sample,GenotypeArray(descr,row)
  else:
    raise ValueError('Uknown format')

  return columns,models,genome,_encode()


# FIXME: Implement COW on genotype models
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

  >>> from reprs import snp
  >>> defmodel  = model_from_alleles('ACGT',allow_hemizygote=True)
  >>> genome = Genome(default_model=defmodel, default_fixed=True)

  >>> samples = ('s1', 's2', 's3')
  >>> genos = [('l1', ['AA', '  ', 'GG']),
  ...          ('l2', ['  ', '',   '  ']),
  ...          ('l3', ['AA', '  ', '  ']),
  ...          ('l4', ['GT', '  ', 'TT'] )]
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
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
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp,genome)
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
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp)
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
  >>> loci,models,genome,new_rows = encode_genomatrixstream_from_strings(loci,genos,'sdat',snp)
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
  >>> samples,models,genome,new_rows = encode_genomatrixstream_from_strings(samples,genos,'ldat',snp,genome)
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

  models = []

  if format=='ldat':
    def _encode():
      n = len(columns)
      m = genome.max_alleles+1

      # FIXME: add support for hemizygote models
      modelcache   = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                               if locus.model is not None and len(locus.model.alleles)==m )
      descrcache   = {}
      strcache     = {}
      to_string    = genorepr.to_string
      from_strings = genorepr.from_strings

      def _genokey(genos):
        gset = set(a for g in set(from_strings(genos)) for a in g)
        gset.discard(None)
        return tuple(sorted(gset))

      for lname,row in genos:
        key = None
        loc = genome.get_locus(lname)

        if loc.model is None:
          key = _genokey(row)
          cachable = genome.default_model is None and len(key) == genome.max_alleles

          if cachable:
            # Aggressively reuse models with fully compatible alleles
            loc.model = modelcache.get(key)

          if loc.model is None:
            genome.get_locus_model(lname)

          if cachable:
            modelcache[key] = loc.model

        model = loc.model

        # FIXME: The semantics of the fixed flag are broken
        if not loc.fixed:
          key = key or _genokey(row)
          try:
            for a in key:
              model.add_allele(a)
          except GenotypeRepresentationError:
            _encoding_error(lname,a,model,warn)
            row = ['']*len(row)

        descr = descrcache.get(model)
        if not descr:
          descr = descrcache[model] = GenotypeArrayDescriptor( [model]*n )

        cache = strcache.get(model)
        if cache is None:
          cache = strcache[model] = dict( (to_string(g),g) for g in model.genotypes )
          for g in genorepr.missing_geno_strs:
            cache[g] = model[None,None]

        models.append(model)

        try:
          yield lname,GenotypeArray(descr,imap(getitem, repeat(cache), row))
        except KeyError:
          gset = set(row)
          try:
            cache.update( (g,model[r]) for g,r in izip(gset,from_strings(gset)) )
            yield lname,GenotypeArray(descr,imap(getitem, repeat(cache), row))
          except KeyError,g:
            _encoding_error(lname,g,model,warn)
            yield lname,GenotypeArray(descr)

  elif format=='sdat':

    genos = sdat_model_lookahead_from_strings(columns,genos,genome,genorepr)

    n = len(columns)
    updates   = []
    cachemap  = {}
    cachelist = []

    to_string    = genorepr.to_string
    from_strings = genorepr.from_strings

    for lname in columns:
      loc = genome.get_locus_model(lname)
      model = loc.model
      models.append(model)
      # FIXME: The semantics of the fixed flag are broken
      if loc.fixed:
        update = model.get_genotype
      else:
        update = model.add_genotype

      cache = cachemap.get(model)
      if cache is None:
        cache = cachemap[model] = dict( (to_string(g),g) for g in model.genotypes )
        cache.update( (genorepr.to_string_from_alleles((g[1],g[0])),g) for g in model.genotypes )
        for g in genorepr.missing_geno_strs:
          cache[g] = model[None,None]

      cachelist.append(cache)
      updates.append( (lname,update,cache) )

    def _encode():
      repr   = genorepr.from_string
      descr  = GenotypeArrayDescriptor(models)
      errors = (KeyError,ValueError,GenotypeRepresentationError)

      for sample,row in genos:
        try:
          row = GenotypeArray(descr,imap(getitem, cachelist, row) )
        except errors:
          geno_tuples = from_strings(row)

          for (lname,update,cache),gstr,g in izip(updates,row,geno_tuples):
            # Aggressively form homozygote genotypes and cache them.  Thus
            # we only see cache misses when we encounter previously
            # unobserved alleles or when genotype formatting is off.
            g1,g2 = g
            if g1!=g2:
              gh = g1,g1
              try:
                cache[genorepr.to_string_from_alleles(gh)] = update(gh)
              except errors:
                pass

              gh = g2,g2
              try:
                cache[genorepr.to_string_from_alleles(gh)] = update(gh)
              except errors:
                pass

            try:
              cache[gstr] = update(g)
            except errors:
              pass

          try:
            row = GenotypeArray(descr,imap(getitem, cachelist, row) )
          except errors:
            _sample_encoding_error(columns,models,geno_tuples,warn)

        yield sample,row
  else:
    raise ValueError('Uknown format')

  return columns,models,genome,_encode()


# FIXME: Implement COW on genotype models
def recode_genotriples(triples,genome):
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
  def _recode():
    updates = {}
    try:
      for sample,lname,geno in triples:
        ud = updates.get(lname)

        if ud is None:
          old_model = geno.model
          old_locus = triples.genome.loci[lname]
          assert old_model is triples.genome.loci[lname].model
          assert old_locus.model is old_model

          genome.merge_locus(lname, fixed=old_locus.fixed, chromosome=old_locus.chromosome,
                                          location=old_locus.location)

          loc = genome.get_locus(lname)
          if loc.model is None:
            loc.model = old_model

          # FIXME: The semantics of the fixed flag are broken
          if old_model is loc.model or loc.fixed:
            ud = loc.model.get_genotype
          else:
            ud = loc.model.add_genotype

          updates[lname] = ud

        yield sample,lname,ud(geno)

    except GenotypeRepresentationError:
      model = triples.genome.get_model(lname)
      _encoding_error(lname,geno,model)

  return triples.clone(_recode(),genome=genome,materialized=False)


# FIXME: Implement COW on genotype models
def encode_genotriples_from_tuples(triples,genome):
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
  >>> for row in encode_genotriples_from_tuples(triples,Genome()):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  try:
    for sample,lname,geno in triples:
      loc = genome.get_locus_model(lname)

      # FIXME: The semantics of the fixed flag are broken
      if loc.fixed:
        geno = loc.model.get_genotype(geno)
      else:
        geno = loc.model.add_genotype(geno)

      yield sample,lname,geno
  except GenotypeRepresentationError:
    _encoding_error(lname,geno,loc.model)


# FIXME: Implement COW on genotype models
def encode_genotriples_from_strings(triples,genorepr,genome):
  '''
  Returns a new genotriples with the genotypes encoded to the a new internal representation

  @param      triples: sequence of genotriples(str,str,genotype representation)
  @type       triples: sequence
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @return            : genotriple in bitpacked format
  @rtype             : genotriple generator

  >>> from reprs import snp
  >>> triples = [('s3', 'l1', 'GG'),('s3', 'l2', 'AA'),
  ...            ('s2', 'l3', 'GT'),('s1', 'l1', 'TT'),
  ...            ('s1', 'l1', 'GG'),('s2', 'l2', 'AA')]
  >>> for row in encode_genotriples_from_strings(triples,snp,Genome()):
  ...   print row
  ('s3', 'l1', ('G', 'G'))
  ('s3', 'l2', ('A', 'A'))
  ('s2', 'l3', ('G', 'T'))
  ('s1', 'l1', ('T', 'T'))
  ('s1', 'l1', ('G', 'G'))
  ('s2', 'l2', ('A', 'A'))
  '''
  local_repr = genorepr.from_string
  if genome is None:
    genome = Genome()

  updates = {}

  try:
    for sample,lname,geno in triples:
      ud = updates.get(lname)

      if ud is None:
        loc = genome.get_locus_model(lname)
        # FIXME: The semantics of the fixed flag are broken
        if loc.fixed:
          ud = lambda g,get=loc.model.get_genotype: get(local_repr(g))
        else:
          ud = lambda g,add=loc.model.add_genotype: add(local_repr(g))

        updates[lname] = ud

      yield sample,lname,ud(geno)

  except GenotypeRepresentationError:
    _encoding_error(lname,geno,loc.model)


def sdat_model_lookahead_from_strings(loci,genos,genome,genorepr,min_unknown=10,max_lookahead=50,warn=False):
  '''
  Lookahead in an sdat genotype stream to determine alleles and fixed
  models.  This is a major optimization that can usually avoid the majority
  of the overhead associated with creating len(loci) default models
  '''
  # Do not attempt to lookahead if a fixed model is used for unknown loci
  if genome.default_fixed:
    return genos

  # Do not attempt to lookahead if only a small number of loci are unknown
  unknown = sum(1 for lname in loci if genome.get_locus(lname).model is None)
  if unknown < min_unknown:
    return genos

  # Track the alleles seen
  alleles_seen = [ [] for i in range(len(loci)) ]

  # Start looking ahead up to max_lookahead rows
  lookahead_rows = []
  for sample,row in genos:
    lookahead_rows.append( (sample,row) )

    row = genorepr.from_strings(row)

    for g,seen in izip(row,alleles_seen):
      seen.append(g[0])
      seen.append(g[1])

    if len(lookahead_rows) == max_lookahead:
      break

  # If there are rows in the lookahead buffer
  if not lookahead_rows:
    return genos

  max_alleles = genome.max_alleles
  modelcache = {}

  try:
    # Review the alleles seen at each locus
    for lname,seen in izip(loci,alleles_seen):
      loc = genome.get_locus(lname)
      model = loc.model

      # If the model is unknown, then check the alleles
      if model is None:
        seen = set(seen)
        seen.discard(None)

        # Create or reuse a fixed model if all alleles have been seen
        # FIXME: add support for hemizygote models
        # FIXME: The semantics of the fixed flag are broken
        if len(seen) == max_alleles:
          seen  = tuple(sorted(seen))
          model = modelcache.get(seen)
          if model is None:
            model = modelcache[seen] = model_from_alleles(seen, max_alleles=max_alleles)
          loc.model = model
          loc.fixed = True
          continue

        # Otherwise create an empty default model
        model = genome.get_model(lname)

      # Populate the observed alleles when a fixed model cannot be used
      for allele in seen:
        model.add_allele(allele)

  except GenotypeRepresentationError:
    _encoding_error(lname,allele,model,warn)

  return chain(lookahead_rows,genos)


def sdat_model_lookahead_from_tuples(loci,genos,genome,min_unknown=10,max_lookahead=25,warn=False):
  '''
  Lookahead in an sdat genotype stream to determine alleles and fixed
  models.  This is a major optimization that can usually avoid the majority
  of the overhead associated with creating len(loci) default models
  '''
  # Do not attempt to lookahead if a fixed model is used for unknown loci
  if genome.default_fixed:
    return genos

  # Do not attempt to lookahead if only a small number of loci are unknown
  unknown = sum(1 for lname in loci if genome.get_locus(lname).model is None)
  if unknown < min_unknown:
    return genos

  # Track the alleles seen
  alleles_seen = [ [] for i in range(len(loci)) ]

  # Start looking ahead up to max_lookahead rows
  lookahead_rows = []
  for sample,row in genos:
    lookahead_rows.append( (sample,row) )

    for g,seen in izip(row,alleles_seen):
      seen.append(g[0])
      seen.append(g[1])

    if len(lookahead_rows) == max_lookahead:
      break

  # If there are rows in the lookahead buffer
  if not lookahead_rows:
    return genos

  max_alleles = genome.max_alleles
  modelcache = {}

  try:
    # Review the alleles seen at each locus
    for lname,seen in izip(loci,alleles_seen):
      loc = genome.get_locus(lname)
      model = loc.model

      # If the model is unknown, then check the alleles
      if model is None:
        seen = set(seen)
        seen.discard(None)

        # Create or reuse a fixed model if all alleles have been seen
        # FIXME: add support for hemizygote models
        # FIXME: The semantics of the fixed flag are broken
        if len(seen) == max_alleles:
          seen  = tuple(sorted(seen))
          model = modelcache.get(seen)
          if model is None:
            model = modelcache[seen] = model_from_alleles(seen, max_alleles=max_alleles)
          loc.model = model
          loc.fixed = True
          continue

        # Otherwise create an empty default model
        model = genome.get_model(lname)

      # Populate the observed alleles when a fixed model cannot be used
      for allele in seen:
        model.add_allele(allele)

  except GenotypeRepresentationError:
    _encoding_error(lname,allele,model,warn)

  return chain(lookahead_rows,genos)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
