# -*- coding: utf-8 -*-
'''
File:          snpstats.py

Authors:       Zhaoming Wang (wangzha@mail.nih.gov)

Created:       August 9, 2006

Abstract:      This utility script generates a file describing counts and
               characteristics of alleles and genotypes for each given SNP

Compatibility: Python 2.4 and above

Requires:      biozilla

Version:       0.99

Revision:      $Id: $

Copyright (c) 2006 BioInformed Consulting Services.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
'''

__version__ = '0.01'

import csv
import sys

from biozilla.genodata  import *
from biozilla.utils     import autofile,hyphen
from biozilla.hwp       import hwp_exact_biallelic
from biozilla.genoarray import snp_acgt
from tagzilla           import read_hapmap_nonfounders


HEADER = ['LOCUS_NAME','REFERENCE_ALLELE','REFERENCE_ALLELE_COUNT','REFERENCE_HOMOZYGOTE_COUNT',
          'OTHER_ALLELE', 'OTHER_ALLELE_COUNT', 'OTHER_HOMOZYGOTE_COUNT', 'HETEROZYGOTE_COUNT',
          'MISSING_ALLELE_COUNT','MISSING_ALLELE_FREQ','MISSING_GENOTYPE_COUNT','HWP_PVALUE',
          'GENOTYPE_COMPLETION_RATE','MAF']


COMPLEMENT = {'A':'T','T':'A','G':'C','C':'G'}


def option_parser():
  import optparse
  usage = 'Usage: %prog [options] ldatfile'

  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-a', '--allelicinfo',    metavar='FILE',
                    help='The file contains the reference allele and other allele for each SNP')
  parser.add_option('-A', '--autosomalsnps',  metavar='FILE',
                    help='The file lists the autosomal SNPs')
  parser.add_option('-S', '--samples',        metavar='FILE',
                    help='The file lists the included samples')
  parser.add_option('-s', '--snps',           metavar='FILE',
                    help='The file lists the included SNPs')
  parser.add_option('-o', '--outfile',        metavar='FILE', default='-',
                    help="The name of the output file(default='-' for standard out)")

  return parser


def load_allele_info(filename):
  rows = csv.reader(autofile(filename))
  allelicinfo = {}
  for row in rows:
    allelicinfo[row[0]] = row[1][1:-1].split('/')
  return allelicinfo


def count(locus,genos,allelicinfo):
  refa,othera = allelicinfo[locus]
  genos = list(genos)
  #complement the alleles if needed
  for geno in genos:
    if geno != '  ' and geno != '':
      genotype = geno
      break
  if refa not in genotype and othera not in genotype:
    refa,othera = COMPLEMENT[refa],COMPLEMENT[othera]

  missing_geno_cnt = 0
  het_cnt = 0
  ref_hom_cnt = 0
  other_hom_cnt = 0
  for geno in genos:
    if refa in geno and othera in geno:
      het_cnt += 1
    elif refa in geno:
      ref_hom_cnt += 1
    elif othera in geno:
      other_hom_cnt += 1
    elif geno == '  ' or geno == '':
      missing_geno_cnt += 1
    else:
      print >> sys.stderr, 'Bad geno',locus,refa,othera,geno,len(geno)

  return refa,othera,missing_geno_cnt,het_cnt,ref_hom_cnt,other_hom_cnt


def completion(missing_geno_cnt,het_cnt,ref_hom_cnt,other_hom_cnt):
  total_geno_cnt = het_cnt + ref_hom_cnt + other_hom_cnt
  total_attempted = total_geno_cnt + missing_geno_cnt
  completion_rate = float(total_geno_cnt) / total_attempted
  return completion_rate


def compute_snp_stats(allelicinfo,autosomalsnps,genomatrix):
  for locus,genos in genomatrix:
    refa,othera,missing_geno_cnt,het_cnt,ref_hom_cnt,other_hom_cnt = count(locus,genos,allelicinfo)
    missing_allele_cnt = 2 * missing_geno_cnt
    if autosomalsnps is not None and locus in autosomalsnps:
      hwp_pvalue = hwp_exact_biallelic(ref_hom_cnt,het_cnt,other_hom_cnt)
    else:
      hwp_pvalue = 1
    completion_rate = completion(missing_geno_cnt,het_cnt,ref_hom_cnt,other_hom_cnt)
    missing_allele_freq = 1 - completion_rate
    refacount = ref_hom_cnt * 2 + het_cnt
    otheracount = other_hom_cnt * 2 + het_cnt
    totala = refacount + otheracount
    maf = min(refacount,otheracount)/float(totala)
    row = [locus,refa,refacount,ref_hom_cnt,othera,otheracount,other_hom_cnt,
           het_cnt,missing_allele_cnt,missing_allele_freq,missing_geno_cnt,hwp_pvalue,completion_rate,maf]
    yield row


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) < 1 or options.allelicinfo is None:
    parser.print_help()
    return

  allelicinfo = load_allele_info(options.allelicinfo)

  genomatrix = load_genomatrix_by_locus(args[0],format='ldat',genorepr=list)

  if options.samples:
    includes = load_list(options.samples,skip=1)
    genomatrix = filter_genomatrix_by_column(genomatrix,includes,False)
  if options.snps:
    includes = load_list(options.snps,skip=1)
    genomatrix = filter_genomatrix_by_row(genomatrix,includes,False)

  genomatrix.next()

  autosomalsnps = None
  if options.autosomalsnps:
    autosomalsnps = load_list(options.autosomalsnps)
  outfile = csv.writer(autofile(hyphen(options.outfile,sys.stdout),'w'),dialect='excel-tab')
  outfile.writerow(HEADER)
  for row in compute_snp_stats(allelicinfo,autosomalsnps,genomatrix):
    outfile.writerow(row)


if __name__ == '__main__':
  main()
