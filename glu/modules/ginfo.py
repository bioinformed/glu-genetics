# -*- coding: utf-8 -*-
'''
File:          ginfo.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       September 5, 2007

Abstract:      Extract available metadata from a GLU genotype file

Compatibility: Python 2.5 and above

Requires:      glu

Version:       0.99

Revision:      $Id: $
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


import sys

from   itertools              import chain

from   glu.lib.fileutils      import table_writer,autofile,namefile,hyphen

from   glu.lib.genolib.io     import load_genostream
from   glu.lib.genolib.phenos import SEX_UNKNOWN,SEX_MALE,SEX_FEMALE, \
                                     PHENO_UNKNOWN,PHENO_UNAFFECTED,PHENO_AFFECTED


def emit(filename,rows):
  table_writer(filename,hyphen=sys.stdout).writerows(rows)


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] genofile'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format',  dest='format', type='str',
                     help='The input genotype file format, possible values=hapmap,ldat,sdat,lbat,sbat,tbat,trip,genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP',
                    help='Input genotype representation')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output results (default is "-" for standard out)')
  parser.add_option('-z', '--lazy', dest='lazy', action='store_true', default=False,
                    help='Be lazy and never materialize the genotypes.  Some results may come back unknown')
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

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                  genome=options.loci,hyphen=sys.stdin)

  out = autofile(hyphen(options.output,sys.stdout),'w')

  if None in (genos.samples,genos.loci) and not options.lazy:
    out.write('Materializing genotypes.\n')
    genos = genos.materialize()

  out.write('Filename    : %s\n' % namefile(args[0]))
  out.write('Format      : %s\n' % genos.format)

  if genos.samples is not None:
    out.write('sample count: %d\n' % len(genos.samples))
    if options.outputsamples:

      def _samples():
        yield ['SAMPLE','FAMILY','INDIVIDUAL','PARENT1','PARENT2','SEX','PHENOCLASS']

        phenome = genos.phenome

        for sample in genos.samples:
          phenos = genos.phenome.get_phenos(sample)
          sex    = {SEX_UNKNOWN:'',SEX_MALE:'MALE',SEX_FEMALE:'FEMALE'}[phenos.sex]
          pheno  = {PHENO_UNKNOWN:'',PHENO_UNAFFECTED:'UNAFFECTED',PHENO_AFFECTED:'AFFECTED'}[phenos.phenoclass]

          yield [sample, phenos.family or '', phenos.individual,
                         phenos.parent1 or '', phenos.parent2 or '',
                         sex, pheno]

      emit(options.outputsamples,_samples())

  if genos.loci is not None:
    out.write('locus  count: %d\n' % len(genos.loci))

    if genos.format == 'genotriple':
      models = [ genos.genome.get_model(locus) for locus in genos.loci ]
    else:
      models = genos.models

    known_models = None not in (genos.loci,genos.samples) or all(m.fixed for m in models)
  else:
    models = None
    known_models = False

  if known_models:
    out.write('model  count: %d\n' % len(set(models)))

    if options.outputloci:
      def _loci():
        yield ['LOCUS','MAX_ALLELES','ALLELES','CHROMOSOME','LOCATION','STRAND']
        for locus in genos.loci:
          loc = genos.genome.get_locus(locus)

          yield [locus, str(loc.model.max_alleles), ','.join(sorted(loc.model.alleles[1:])),
                        loc.chromosome or '', loc.location or '', loc.strand or '']

      emit(options.outputloci,_loci())

  elif genos.loci is not None and options.outputloci:
    emit(options.outputloci, chain([['LOCUS']],([l] for l in genos.loci)))


if __name__=='__main__':
  main()
