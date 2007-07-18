'''
  Split the single genotype file by chr
  Input:
    (1) SNP to Chr map file (snpidchr.map)
    (2) The directory for output files
    (3) The original genotype file
  Output:
    (1) Split by chr
'''
import csv
import sys

from biozilla.utils    import autofile
from biozilla.genodata import load_map


def main():
  chrs = map(str,range(1,23))
  chrs.extend(['X','XY','Y','M'])

  snpchr = load_map(sys.argv[1],skip=1)

  #open the file handles and store them in a dictionary
  filedict = {}
  for chr in chrs:
    filename = '%s/genotypes_chr%s.txt.gz' % (sys.argv[2],chr)
    w = csv.writer(autofile(filename,'w'),dialect='excel-tab')
    filedict[chr] = w

  genos = csv.reader(autofile(sys.argv[3]),dialect='excel-tab')

  header = genos.next()
  for w in filedict.itervalues():
    w.writerow(header)

  for row in genos:
    snpid = row[0]
    chr = snpchr[snpid]
    w = filedict[chr]
    w.writerow(row)


if __name__=='__main__':
  main()

