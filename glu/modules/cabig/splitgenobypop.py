# -*- coding: utf-8 -*-
'''
File:          splitgenobypop.py.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

'''
  Split the single genotype file by population
  Input:
    (1) Sample(spcimen) to population map file (sample.def)
    (2) The directory for output files
    (3) The original genotype file
  Output:
    (2) Split by pop
'''
import csv
import sys

from   glu.lib.utils import autofile
from   utils         import load_map


def main():
  pops = ['CONTROL','CASE_NONAGGRESSIVE','CASE_AGGRESSIVE','CEPH']

  sampop = load_map(sys.argv[1],12,8,1)
  for key,value in sampop.iteritems():
    if value == 'QC':
      sampop[key] = 'CEPH'

  #open the file handles and store them in a dictionary
  filedict = {}
  for pop in pops:
    filename = '%s/genotypes_%s.txt.gz' % (sys.argv[2],pop)
    w = csv.writer(autofile(filename,'w'),dialect='excel-tab')
    filedict[pop] = w

  genos = csv.reader(autofile(sys.argv[3]),dialect='excel-tab')

  header = genos.next()
  for w in filedict.itervalues():
    w.writerow(header)

  for row in genos:
    specid = row[2]
    pop = sampop[specid]
    w = filedict[pop]
    w.writerow(row)


if __name__=='__main__':
  main()

