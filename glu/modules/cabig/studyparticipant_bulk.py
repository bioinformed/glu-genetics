# -*- coding: utf-8 -*-
'''
File:          studyparticipant_bulk.py

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
from   snpfreqfact_bulk  import load_study_pop


AGE={'0':'55-59','1':'60-64','2':'65-69','3':'70-74','':''}
HEADER = ['Study Participant de-identifier ID','Gender','Age','Affection Status','Family History','Population']


def main():
  popmap = load_study_pop(sys.argv[1])
  #read the phenotyic file
  r = csv.reader(autofile(sys.argv[2]))
  r.next()
  out = csv.writer(autofile(sys.argv[3],'w'),dialect='excel-tab')
  out.writerow(HEADER)

  for row in r:
    pid = row[0]
    age = AGE[row[1]]
    gender = row[3]
    status = row[7]
    famhist = row[6]
    pop = popmap[row[4]]
    out.writerow([pid,gender,age,status,famhist,pop])

if __name__=='__main__':
  main()
