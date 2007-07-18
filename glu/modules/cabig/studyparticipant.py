# -*- coding: utf-8 -*-
'''
File:          studyparticipant.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

'''
  Inputs:
  1. sample.def # list of samples for CGEMS Prostate Scan 1
  2. hapmap_ceph_id.out # CGEMS CEPH individuals also typed in Hapmap

  Ouputs:
  1. sample data file for the table of study_participant
  2. participant id map file
  3. scan1a individual id
  4. scan1b individual id
'''

import csv
import sys

from   itertools         import count

from   glu.lib.fileutils import autofile


HEADER = ['PARTICIPANT_PK','PARTICIPANT_DID','AGE_AT_ENROLL','ETHNIC_GROUP_CODE','GENDER','POPULATION_ID',
          'STUDY_NAME','FAMILY_HISTORY','CASE_CONTROL_STATUS']

POPMAP = {'CONTROL':'6','CASE_AGGRESSIVE':'7','CASE_NONAGGRESSIVE':'8','QC':'9'}


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-s', '--samdef',       dest='samdef',        metavar='FILE',   default= '-',
                    help='the generated samdef file')
  parser.add_option('-i', '--hapmapcephs',  dest='hapmapcephs',   metavar='FILE',   default= '-',
                    help='CGEMS CEPH individuals typed in Hapmap')
  parser.add_option('-o', '--outdatafile',  dest='outdatafile',   metavar='FILE',   default= '-',
                    help='Output file for table study_participant')
  parser.add_option('-m', '--outmapfile',   dest='outmapfile',    metavar='FILE',   default= '-',
                    help='Output a map file between pid and deidentified pid')
  parser.add_option('-a', '--scan1aout',    dest='scan1aout',     metavar='FILE',   default= '-',
                    help='Output all pids from scan1A')
  parser.add_option('-b', '--scan1bout',    dest='scan1bout',     metavar='FILE',   default= '-',
                    help='Output all pids from scan1B')
  parser.add_option('-c', '--inicount',     dest='inicount',      type = 'int')

  return parser



def main():

  parser = option_parser()
  options,args = parser.parse_args()

  outdatafile = csv.writer(autofile(options.outdatafile,'w'))
  outdatafile.writerow(HEADER)
  outmapfile  = csv.writer(autofile(options.outmapfile,'w'),dialect='excel-tab')
  hapmapcephs = set(row[0] for row in csv.reader(autofile(options.hapmapcephs)))

  samdef = csv.reader(autofile(options.samdef),dialect='excel-tab')
  samdef.next()
  #set global attributes
  studyname = 'CGEMS Prostate Cancer WGAS Phase 1'
  ethnic = 'CAUCASIAN'

  scan1aout = csv.writer(autofile(options.scan1aout,'w'),dialect='excel-tab')
  scan1bout = csv.writer(autofile(options.scan1bout,'w'),dialect='excel-tab')

  participantids = set()
  seqnum = count(options.inicount)
  for row in samdef:
    if row[13] != 'YES':
      continue
    if row[7] == 'CEPH' and row[3] not in hapmapcephs:
      continue
    if row[9] == 'HumanHap317K':
      scan1aout.writerow([row[3]])
    else:
      scan1bout.writerow([row[3]])
    participantid = row[12]
    if participantid in participantids:
      continue
    participantids.add(participantid)

    age      = row[17]
    gender   = row[6]
    popid    = POPMAP[row[8]]
    famhist  = row[10]
    ccstatus = row[8]

    outdatafile.writerow([seqnum.next(),participantid,age,ethnic,gender,popid,studyname,famhist,ccstatus])
    outmapfile.writerow([row[0],participantid])


if __name__=='__main__':
  main()

