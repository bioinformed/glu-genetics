# -*- coding: utf-8 -*-
'''
File:          genotypefact_bulk.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   glu.lib.utils import autofile
from   utils         import load_map

HEADER = ['dbSNP ID','Study Participant de-identifier ID','Specimen de-identified ID',
          'Allele 1','Allele 2','Quality Score','QC Status',
          'Normalized X Intensity','Normalized Y Intensity','Raw X Intensity','Raw Y Intensity']


def main():
  smap = load_map(sys.argv[1],1,0,1,dialect='excel')
  snpmap = load_map(sys.argv[2],1,0,1,dialect='excel')
  genos = csv.reader(autofile(sys.argv[3]))
  genos.next()
  out = csv.writer(autofile(sys.argv[4],'w'),dialect='excel-tab')
  out.writerow(HEADER)
  for row in genos:
    a1,a2,qscore,nx,ny,rx,ry,assay,sid,annoid,studyname,scanid,status=row
    dbsnpid = snpmap[annoid]
    pid = smap[sid]
    out.writerow([dbsnpid,pid,sid,a1,a2,qscore,status,nx,nx,rx,ry])

if __name__ == '__main__':
  main()

