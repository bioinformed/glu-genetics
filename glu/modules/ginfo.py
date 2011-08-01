# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Print information about a GLU genotype file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.lib.fileutils      import table_writer,autofile,namefile,hyphen

from   glu.lib.genolib.io     import load_genostream, geno_options
from   glu.lib.genolib.phenos import SEX_UNKNOWN,SEX_MALE,SEX_FEMALE, \
                                     PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


def emit(filename,rows):
  table_writer(filename,hyphen=sys.stdout).writerows(rows)


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('genofile', nargs='+', help='Genotype file(s)')

  geno_options(parser,input=True,filter=True)

  parser.add_argument('-z', '--lazy', action='store_true', default=False,
                      help='Be lazy and never materialize the genotypes.  Some results may come back unknown.')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                      help='Output results (default is "-" for standard out)')
  parser.add_argument('--outputloci',   metavar='FILE',
                       help='Output the list of loci to FILE')
  parser.add_argument('--outputsamples',metavar='FILE',
                       help='Output the list of samples to FILE')
  return parser


def ginfo(filename,options,out):
  genos = load_genostream(filename,format=options.informat,genorepr=options.ingenorepr,
                                  genome=options.loci,phenome=options.pedigree,
                                  transform=options,hyphen=sys.stdin)


  if None in (genos.samples,genos.loci) and not options.lazy:
    out.write('Materializing genotypes.\n')
    genos = genos.materialize()

  out.write('Filename    : %s\n' % namefile(filename))
  out.write('Format      : %s\n' % genos.format)

  if genos.samples is not None:
    out.write('sample count: %d\n' % len(genos.samples))

  if genos.loci is not None:
    out.write('locus  count: %d\n' % len(genos.loci))

  out.write('\n')

  if options.outputsamples:
    def _samples():
      yield ['SAMPLE','FAMILY','INDIVIDUAL','PARENT1','PARENT2','SEX','PHENOCLASS']

      phenome   = genos.phenome
      sex_map   = {SEX_UNKNOWN:'',SEX_MALE:'MALE',SEX_FEMALE:'FEMALE'}
      pheno_map = {PHENO_UNKNOWN:'',PHENO_UNAFFECTED:'UNAFFECTED',PHENO_AFFECTED:'AFFECTED'}

      for sample in genos.samples:
        phenos = phenome.get_phenos(sample)
        sex    = sex_map[phenos.sex]
        pheno  = pheno_map[phenos.phenoclass]

        yield [sample, phenos.family or '', phenos.individual,
                       phenos.parent1 or '', phenos.parent2 or '',
                       sex, pheno]

    emit(options.outputsamples,_samples())


  if options.outputloci:
    # FIXME: Add information to streams about known models
    if not genos.materialized: # not genos.known_models:
      for geno in genos:
        pass

    def _loci():
      yield ['LOCUS','MAX_ALLELES','ALLELES','CHROMOSOME','LOCATION','STRAND']
      for locus in genos.loci:
        loc = genos.genome.get_locus(locus)
        assert loc.model is not None
        yield [locus, str(loc.model.max_alleles), ','.join(sorted(loc.model.alleles[1:])),
                      loc.chromosome or '', loc.location or '', loc.strand or '']

    emit(options.outputloci,_loci())


def main():
  parser  = option_parser()
  options = parser.parse_args()

  out = autofile(hyphen(options.output,sys.stdout),'w')

  if len(options.genofile) > 1 and (options.outputsamples or options.outputloci):
    raise ValueError('Cannot output samples or loci when multiple input files are specified')

  for filename in options.genofile:
    ginfo(filename,options,out)


if __name__=='__main__':
  main()
