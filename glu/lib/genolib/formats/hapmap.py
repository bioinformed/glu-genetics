# -*- coding: utf-8 -*-
'''
File:          hapmap.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU text genotype format input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   itertools                 import islice,dropwhile

from   glu.lib.fileutils         import autofile,namefile,tryint,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import model_from_alleles
from   glu.lib.genolib.reprs     import hapmap
from   glu.lib.genolib.locus     import Genome

__all__ = ['load_hapmap']


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']


def load_hapmap(filename,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a HapMap genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  unique = get_arg(args, ['unique'], True)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  gfile = autofile(filename)
  gfile = dropwhile(lambda s: s.startswith('#'), gfile)

  try:
    header = gfile.next()
  except StopIteration:
    header = ''

  if not any(header.startswith(h) for h in HAPMAP_HEADERS):
    raise ValueError("Input file '%s' does not appear to be in HapMap format." % namefile(filename))

  columns = [ intern(h) for h in islice(header.split(),11,None) ]

  if genome is None:
    genome = Genome()

  # FIXME: Add support for hemizygote models
  m = genome.max_alleles+1
  modelcache = dict( (tuple(sorted(locus.model.alleles[1:])),locus.model) for locus in genome.loci.itervalues()
                           if locus.model is not None and len(locus.model.alleles)==m )

  def _load_hapmap():
    n = len(columns)
    for i,line in enumerate(gfile):
      fields = [ f.strip() for f in line.split(' ') ]

      # FIXME: Issue warning or raise error
      if len(fields) < 11:
        sys.stderr.write('Invalid HapMap data row at %s:%d.  Expected >11 records, found %d.' \
                                      % (namefile(filename),i+1,len(fields)))
        continue

      locus      = intern(fields[0].strip())
      alleles    = tuple(sorted(fields[1].split('/')))
      chromosome = fields[2].strip()
      position   = tryint(fields[3].strip())
      strand     = intern(fields[4].strip())
      genos      = fields[11:]

      if len(genos) != n:
        raise ValueError('Invalid genotype length in %s:%d for locus %s.  Expected %d genotypes, found %d.' \
                              % (namefile(filename),i+1,locus,n,len(genos)))

      # Normalize 'chrXX' names to just 'XX'
      if chromosome.startswith('chr'):
        chromosome = chromosome[3:].strip()

      chromosome = intern(chromosome)

      if len(alleles) != 2 or any(a not in 'ACGT' for a in alleles):
        alleles = tuple(set(a for g in genos for a in g if a!='N'))

      # FIXME: Add error recovery and detection
      assert len(alleles)<=2

      model = genome.get_locus(locus).model

      if not model:
        model = modelcache.get(alleles)

      if not model:
        model = model_from_alleles(alleles,max_alleles=genome.max_alleles)

        if genome.default_model is None and len(alleles) == genome.max_alleles:
          modelcache[alleles] = model

      genome.merge_locus(locus, model, False, chromosome, position, strand)

      yield locus,genos

  genos = GenomatrixStream.from_strings(_load_hapmap(),'ldat',genorepr=hapmap,samples=columns,
                                                              genome=genome,phenome=phenome,unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
