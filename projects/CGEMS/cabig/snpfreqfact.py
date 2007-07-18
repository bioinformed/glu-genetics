'''
  1. freqfiles: the output files from snpstats.py
  2. snpidmap:  dbsnp id (rs number) to snp_anno_id(randomized sequence number) map
  3. specify the population id from the option
'''

__version__ = '0.01'

import csv
import sys

from itertools             import chain,islice,count
from biozilla.genodata     import load_genomatrix
from biozilla.fileutils    import autofile,hyphen,load_map
from tagzilla              import TagZillaOptionParser


HEADER = ['SNP_FREQ_ID','SNPANNO_ID','POPULATION_ID','REF_ALLELE','REF_ALLELE_COUNT','REF_HOMOZYGOTE_COUNT',
          'OTHER_ALLELE', 'OTHER_ALLELE_COUNT', 'OTHER_HOMOZYGOTE_COUNT', 'HETEROZYGOTE_COUNT',
          'MISSING_ALLELE_COUNT','MISSING_ALLELE_FREQ','MISSING_GENOTYPE_COUNT','HWP_PVALUE',
          'GENOTYPE_COMPLETION_RATE','MAF']

def option_parser():
  usage = 'Usage: %prog [options] freqfiles...'

  parser = TagZillaOptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-s', '--snpidmap', metavar='FILE',
                    help='The file containing the map from snp id to snp annotation id')
  parser.add_option('-o', '--outfile',  metavar='FILE', default='-',
                    help="The name of the output file(default='-' for standard out)")
  parser.add_option('-p', '--popid',    type='int')
  parser.add_option('-c', '--inicount', type='int')

  return parser

def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) < 1 or options.snpidmap is None:
    parser.print_help()
    return

  snpidmap = load_map(options.snpidmap)

  outfile = csv.writer(autofile(hyphen(options.outfile,sys.stdout),'w'))
  outfile.writerow(HEADER)

  freqid = count(options.inicount)
  for file_options,filename in args:
    snpstats = csv.reader(autofile(filename),dialect='excel-tab')
    snpstats.next()

    print len(snpidmap)
    for row in snpstats:
      #print row[0],snpidmap[row[0]]
      outfile.writerow(list(chain([freqid.next(),snpidmap[row[0]],file_options.popid],islice(row,1,None))))


if __name__ == '__main__':
  main()

