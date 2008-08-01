# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Print information about a GLU genotype file'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   itertools              import chain

from   glu.lib.fileutils      import table_writer,autofile,namefile,hyphen

from   glu.lib.genolib.io     import load_genostream, geno_options
from   glu.lib.genolib.phenos import SEX_UNKNOWN,SEX_MALE,SEX_FEMALE, \
                                     PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


def emit(filename,rows):
  table_writer(filename,hyphen=sys.stdout).writerows(rows)


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] genofile'

  parser = optparse.OptionParser(usage=usage)

  geno_options(parser,input=True)

  parser.add_option('-z', '--lazy', dest='lazy', action='store_true', default=False,
                    help='Be lazy and never materialize the genotypes.  Some results may come back unknown')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  parser.add_option('--outputloci',   dest='outputloci',    metavar='FILE',
                     help='Output the list of loci to FILE')
  parser.add_option('--outputsamples',dest='outputsamples', metavar='FILE',
                     help='Output the list of samples to FILE')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  genos = load_genostream(args[0],format=options.informat,genorepr=options.ingenorepr,
                                  genome=options.loci,phenome=options.pedigree,hyphen=sys.stdin)

  out = autofile(hyphen(options.output,sys.stdout),'w')

  if None in (genos.samples,genos.loci) and not options.lazy:
    out.write('Materializing genotypes.\n')
    genos = genos.materialize()

  out.write('Filename    : %s\n' % namefile(args[0]))
  out.write('Format      : %s\n' % genos.format)

  if genos.samples is not None:
    out.write('sample count: %d\n' % len(genos.samples))

  if genos.loci is not None:
    out.write('locus  count: %d\n' % len(genos.loci))

  if options.outputsamples:
    def _samples():
      yield ['SAMPLE','FAMILY','INDIVIDUAL','PARENT1','PARENT2','SEX','PHENOCLASS']

      phenome   = genos.phenome
      sex_map   = {SEX_UNKNOWN:'',SEX_MALE:'MALE',SEX_FEMALE:'FEMALE'}
      pheno_map = {PHENO_UNKNOWN:'',PHENO_UNAFFECTED:'UNAFFECTED',PHENO_AFFECTED:'AFFECTED'}

      for sample in genos.samples:
        phenos = genos.phenome.get_phenos(sample)
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


if __name__=='__main__':
  main()
