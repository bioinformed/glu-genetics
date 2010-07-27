# -*- coding: utf-8 -*-

__abstract__  = 'VCF parser'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['load_vcf']

__genoformats__ = [
  #    LOADER     SAVER   WRITER  PFORMAT    ALIAS    EXTS
  ('load_vcf',    None,    None,   'vcf',    None,   'vcf') ]


import sys
import csv

from   itertools                 import takewhile

from   glu.lib.fileutils         import autofile,namefile,tryint,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.genoarray import build_model
from   glu.lib.genolib.locus     import Genome


VCF_HEADER = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']


def load_vcf(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
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

  gfile   = csv.reader(autofile(filename),dialect='excel-tab')
  info    = list(takewhile(lambda l: l and l[0].startswith('##'), gfile))
  headers = list(takewhile(lambda l: l and l[0].startswith('#'),  gfile))

  header = headers[-1]

  n = len(VCF_HEADER)

  if header[:n] != VCF_HEADER:
    raise ValueError("Input file '%s' does not appear to be in VCF format." % namefile(filename))

  samples = map(intern,header[n:])

  file_genome = Genome()

  models = []
  loci   = []

  def _load_vcf():
    modelcache = {}
    n = len(samples)
    strand = '+'

    for fields in gfile:
      if len(fields) < 9:
        sys.stderr.write('Invalid VCF data row at %s:%d.  Expected >8 records, found %d.' \
                                      % (namefile(filename),gfile.line_num+1,len(fields)))
        continue

      chromosome = fields[0]
      location   = tryint(fields[1])
      locus      = fields[2]
      ref_allele = intern(fields[3])
      alt_allele = intern(fields[4])

      if ref_allele=='.':
        ref_allele='-'

      if alt_allele=='.':
        alt_allele='-'

      alleles = ref_allele,alt_allele

      if locus=='.':
        locus = 'SNP%s-%s' % (chromosome,location)

      # Skip non-biallelic loci for now
      if ',' in alt_allele:
        continue

      genos = fields[9:]
      if len(genos) != n:
        raise ValueError('Invalid genotype length in %s:%d for locus %s.  Expected %d genotypes, found %d.' \
                              % (namefile(filename),gfile.line_num+1,locus,n,len(genos)))

      # Normalize 'chrXX' names to just 'XX'
      if chromosome.startswith('chr'):
        chromosome = chromosome[3:].strip()

      chromosome = intern(chromosome)

      if len(alleles) != 2 or any(a not in 'ACGT-' for a in alleles):
        alleles = tuple(set(a for g in genos for a in g if a!='N'))

      # FIXME: Add error recovery and detection
      assert len(alleles)<=2

      model = modelcache.get(alleles)

      if not model:
        model = modelcache[alleles] = build_model(alleles,max_alleles=file_genome.max_alleles)

      modelmap = {'0|0'  : model[ref_allele,ref_allele],
                  '0/0'  : model[ref_allele,ref_allele],
                  '0\\0' : model[ref_allele,ref_allele],
                  '0|1'  : model[ref_allele,alt_allele],
                  '0/1'  : model[ref_allele,alt_allele],
                  '0\\1' : model[ref_allele,alt_allele],
                  '1|0'  : model[alt_allele,ref_allele],
                  '1/0'  : model[alt_allele,ref_allele],
                  '1\\0' : model[alt_allele,ref_allele],
                  '1|1'  : model[alt_allele,alt_allele],
                  '1/1'  : model[alt_allele,alt_allele],
                  '1\\1' : model[alt_allele,alt_allele]}

      genos = [ modelmap[g.split(':')[0]] for g in genos ]
      file_genome.merge_locus(locus, model, chromosome, location, strand)

      loci.append(locus)
      models.append(model)

      yield locus,genos

  genos = GenomatrixStream(_load_vcf(),'ldat',samples=samples,loci=loci,models=models,
                                              genome=file_genome,phenome=phenome,unique=unique,
                                              packed=False)

  if genome:
    genos = genos.transformed(recode_models=genome)

  if unique:
    genos = genos.unique_checked()

  return genos


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
