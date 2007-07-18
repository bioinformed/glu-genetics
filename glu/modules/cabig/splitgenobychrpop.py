# -*- coding: utf-8 -*-
'''
File:          splitgenobychrpop.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

'''
  Split the single genotype file by chr and or population
  Input:
    (1) SNP to Chr map file (snpidchr.map)
    (2) Sample(spcimen) to population map file (sample.def)
    (3) The directory for output files
    (4) The original genotype file
  Output:
    (1) Split by chr and pop
'''
import csv
import sys

from   glu.lib.utils import autofile
from   utils         import load_map


def main():
  chrs = map(str,range(1,23))
  chrs.extend(['X','XY','Y','M'])

  pops = ['CONTROL','CASE_NONAGGRESSIVE','CASE_AGGRESSIVE','CEPH']

  snpchr = load_map(sys.argv[1],0,1,1)
  sampop = load_map(sys.argv[2],12,8,1)
  for key,value in sampop.iteritems():
    if value == 'QC':
      sampop[key] = 'CEPH'

  #open the file handles and store them in a dictionary
  filedict = {}
  for chr in chrs:
    for pop in pops:
      filename = '%s/genotypes_%s_chr%s.txt.gz' % (sys.argv[3],pop,chr)
      w = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
      filedict[(pop,chr)] = w

  genos = csv.reader(autofile(sys.argv[4]),dialect='excel-tab')

  header = genos.next()
  for w in filedict.itervalues():
    w.writerow(header)

  for row in genos:
    specid = row[2]
    snpid = row[0]
    chr = snpchr[snpid]
    pop = sampop[specid]
    w = filedict[(pop,chr)]
    w.writerow(row)


if __name__=='__main__':
  main()

