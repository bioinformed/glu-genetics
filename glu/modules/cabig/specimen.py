# -*- coding: utf-8 -*-
'''
File:          specimen.py

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

from   glu.lib.fileutils import autofile,load_list


HEADER = ['PARTICIPANT_ID','SPECIMEN_ID','SPECIMEN_TYPE']

def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',        dest='outfile',        metavar='FILE',   default= '-',
                    help='Output file for table specimen')
  parser.add_option('-s', '--samdef',         dest='samdef',         metavar='FILE',   default= '-',
                    help='the generated sample def file')
  parser.add_option('-i', '--participants',   dest='participants',   metavar='FILE',   default= '-',
                    help='the generated csv file for table study_participant')

  return parser


def main():

  parser = option_parser()
  options,args = parser.parse_args()

  samdef = csv.reader(autofile(options.samdef),dialect='excel-tab')
  samdef.next()
  participants = load_list(options.participants)
  outfile = csv.writer(autofile(options.outfile,'w'))
  outfile.writerow(HEADER)

  specimentype = 'DNA'

  for row in samdef:
    if row[14] == 'TRUE' and  row[13] in participants:
      outfile.writerow([row[13],row[12],specimentype])


if __name__=='__main__':
  main()
