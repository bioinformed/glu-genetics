# -*- coding: utf-8 -*-

__abstract__  = 'HapMap genotype parser'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['load_hapmap']

__genoformats__ = [
  #    LOADER     SAVER   WRITER  PFORMAT    ALIAS    EXTS
  ('load_hapmap', None,    None,   'ldat', 'hapmap',  None) ]


import sys

from   itertools                 import islice,dropwhile

from   glu.lib.fileutils         import autofile,namefile,tryint,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import build_model
from   glu.lib.genolib.reprs     import hapmap
from   glu.lib.genolib.locus     import Genome


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']


def load_hapmap(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
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

  def _load_hapmap():
    modelcache = {}
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

      model = genome.get_model(locus)

      if not model:
        model = modelcache.get(alleles)

      if not model:
        model = modelcache[alleles] = build_model(alleles,max_alleles=genome.max_alleles)

      genome.merge_locus(locus, model, chromosome, position, strand)

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
