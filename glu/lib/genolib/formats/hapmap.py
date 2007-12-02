# -*- coding: utf-8 -*-
'''
File:          hapmap.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2006-01-01

Abstract:      GLU text genotype format input/output objects

Requires:      Python 2.5

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


from   itertools                 import islice,dropwhile

from   glu.lib.fileutils         import autofile,namefile,tryint

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import model_from_alleles
from   glu.lib.genolib.reprs     import hapmap
from   glu.lib.genolib.locus     import Genome

__all__ = ['load_genomatrix_hapmap']


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']


def load_genomatrix_hapmap(filename,limit=None,genome=None):
  '''
  Load a HapMap genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param        limit: limit the number of samples loaded. Default is None
  @type         limit: int or None
  @rtype             : GenomatrixStream
  '''
  gfile = autofile(filename)
  gfile = dropwhile(lambda s: s.startswith('#'), gfile)

  try:
    header = gfile.next()
  except StopIteration:
    header = []

  if not any(header.startswith(h) for h in HAPMAP_HEADERS):
    raise ValueError("Input file '%s' does not appear to be in HapMap format." % namefile(filename))

  if limit is not None:
    limit += 11

  columns = [ intern(h.strip()) for h in islice(header.split(),11,limit) ]
  modelcache = {}

  if genome is None:
    genome = Genome()

  def _load():
    n = len(columns)
    for line in gfile:
      fields     = line.split()
      locus      = intern(fields[0].strip())
      alleles    = tuple(sorted(fields[1].split('/')))
      chromosome = fields[2].strip()
      position   = tryint(fields[3].strip())
      strand     = intern(fields[4].strip())
      genos      = fields[11:limit]

      # Normalize 'chrXX' names to just 'XX'
      if chromosome.startswith('chr'):
        chromosome = chromosome[3:].strip()

      chromosome = intern(chromosome)

      if len(alleles) != 2 or any(a not in 'ACGT' for a in alleles):
        alleles = tuple(set(a for g in genos for a in g if a!='N'))

      # FIXME: Add error recovery and detection
      assert len(alleles)<=2
      assert len(genos) == n

      model = genome.get_locus(locus).model

      if not model:
        model = modelcache.get(alleles)

      if not model:
        model = modelcache[alleles] = model_from_alleles(alleles,max_alleles=2)

      loc = genome.merge_locus(locus, model, True, chromosome, position, strand)

      yield locus,genos

  return GenomatrixStream.from_strings(_load(),'ldat',hapmap,samples=columns,genome=genome)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
