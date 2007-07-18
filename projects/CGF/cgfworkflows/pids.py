# -*- coding: utf-8 -*-
'''
File:

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      A utility for generating PID-based sample info.
               Used for generating analysis-ready-data set

Requires:      Python 2.4

Revision:
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv
from   itertools              import islice,imap,chain,groupby
from   operator               import itemgetter

from   biozilla.utils         import autofile,hyphen
from   biozilla.genoarray     import snp_acgt
from   biozilla.genodata      import  *
from   app_utils              import  mapper,mapper2
import re


def get_pids(filename):
  rows = csv.reader(autofile(filename),dialect='excel-tab')
  rows.next()

  adict = {}
  for r in rows:
    pid   = r[1].strip()
    status= r[2].strip()
    status= status.lower()
    if not pid:
      continue
    adict.setdefault(pid,{})[status] = None
    adict[pid].setdefault('sid',[]).append(r[0])

  for it in adict:
    if 'qc' in adict[it]:
      del adict[it]['qc']

  return adict


def get_output(adict,outfile,sidout):
  out = csv.writer(autofile(outfile, 'w'),dialect='excel-tab')
  out.writerow(['PID','PHENOGROUP','PHENOGROUP_DESC'])

  out2 = csv.writer(autofile(sidout, 'w'),dialect='excel-tab')
  pids = sorted(adict)
  for it in pids:
    sids = adict[it].pop('sid')
    if adict[it]:
      for sid in sids:
        out2.writerow([sid,it])
      if 'control' in ';'.join(adict[it].keys()):
        out.writerow([it,'0','CONTROL'])
      else:
        out.writerow([it,'1',';'.join(sorted(adict[it])).upper()])


def option_parser():
  import optparse

  usage  = 'usage: %prog inputfile outputfile'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option(      '--sidout', dest='sidout', metavar='FILE',
                    help='Specify an output file of sid include')
  #parser.add_option(     '--dupskip', dest='dupskip', metavar='N', type='int', default=0,
  #                  help='specifying number of skipped rows in the duplicates file (optional)')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  adict = get_pids(args[0])
  get_output(adict,args[1],options.sidout)


if __name__ == '__main__':
  main()
