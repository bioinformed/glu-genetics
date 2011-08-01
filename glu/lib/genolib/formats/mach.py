# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'MACH imputed genotype format input objects'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['load_mach']

__genoformats__ = [
  #   LOADER       SAVER          WRITER       PFORMAT           ALIAS              EXTS
  ('load_mach',    None,          None,        'sdat',           'mach',            'geno') ]


from   itertools                 import islice

from   glu.lib.utils             import gcdisabled
from   glu.lib.fileutils         import autofile,namefile,table_reader,  \
                                        guess_related_file,related_file, \
                                        parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.locus     import Genome, load_locus_records, populate_genome
from   glu.lib.genolib.phenos    import Phenome,SEX_MALE,SEX_FEMALE,SEX_UNKNOWN
from   glu.lib.genolib.genoarray import GenotypeArray,build_model,build_descr


def load_mach_info(filename,genome):
  mfile = table_reader(filename)

  header = next(mfile)

  if header!=['SNP','Al1','Al2','Freq1','MAF','Quality','Rsq']:
    raise ValueError('Invalid MACH imputed genotype info header')

  modelcache = {}

  for i,row in enumerate(mfile):
    if len(row) != 7:
      raise ValueError('Invalid MACH info record %d' % (i+1))

    lname   = row[0]
    alleles = tuple(row[1:3])

    model = modelcache.get(alleles)
    if model is None:
      model = modelcache[alleles] = build_model(alleles=alleles,max_alleles=2)

    genome.merge_locus(lname,model)

    yield lname


def load_mach(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a MACH format genotype data file.

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
  unique   = get_arg(args, ['unique'], True)
  info     = get_arg(args, ['info']) or guess_related_file(filename,['info'])

  if info is None:
    raise ValueError('INFO file must be specified when loading MACH files')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if genome is None:
    genome = Genome()

  if phenome is None:
    phenome = Phenome()

  loci = list(load_mach_info(info,genome))

  gfile = autofile(filename)

  def _load_mach():
    n = len(loci)+2

    gmap = dict( ('%s/%s' % (a1,a2), (a1,a2)) for a1 in 'ACGT' for a2 in 'ACGT' )
    gmap['N/N'] = None

    for line_num,line in enumerate(gfile):
      with gcdisabled():
        fields = line.split()

        if len(fields) != n:
          raise ValueError('Invalid record on line %d of %s' % (line_num+1,namefile(filename)))

        name    = fields[0]
        rectype = fields[1]
        genos   = [ gmap[g] for g in islice(fields,2,None) ]

      yield name,genos

  genos = GenomatrixStream.from_tuples(_load_mach(),'sdat',loci=loci,genome=genome,phenome=phenome,unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
