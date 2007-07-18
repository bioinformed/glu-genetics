# -*- coding: utf-8 -*-
'''
File:          sampledef_prep.py

Authors:       Jun Lu (lujun@mail.nih.gov)

Created:

Abstract:      A utility for generating .map files where modifying sample/locus ids are necessary.
               It is necessary to modify mapper.py first before to run this scripts

Requires:      Python 2.5

Revision:
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   operator         import itemgetter

from   glu.lib.genodata import load_map
from   glu.lib.utils    import autofile


def create_new_sampledef(infile,outfile,reason_map,rankfile):
  report_map = {}
  rows = csv.reader(autofile(rankfile),dialect='excel-tab')

  report_map = dict((r[1],(r[0],r[2])) for r in rows if r[0])

  rows = csv.reader(autofile(infile),dialect='excel-tab')
  out  = csv.writer(autofile(outfile,'w'),dialect='excel-tab')
  out.writerow(rows.next())

  for r in rows:
    if r[8]in ('TRUE','QC'):
      out.writerow(r)
      continue
    reason = [reason_map.get(v,v) for v in r[8].split(';')]
    n_r = []
    for v in reason:
      if v in report_map:
        n_r.append(report_map[v])
      else:
        print r
        raise ValueError,'The value "%s" is not on the reason list' % v
    n_r.sort(key=lambda x: int(x[0]))
    r[8] = n_r[0][1] if n_r else ''
    out.writerow(r)

  return None


def option_parser():
  import optparse

  usage = 'usage: %prog [options] sampledef_in sampledef_out'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-m', '--reasonfile', dest='reasonfile', metavar='FILE',
                    help='file for mapping old exclude-reason list to a new list')

  parser.add_option('-r', '--rankfile', dest='rankfile', metavar='FILE',
                    help='file for ranking excldue-reason list')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  reason_map = {}
  if options.reasonfile:
    reason_map = load_map(options.reasonfile)

  create_new_sampledef(args[0],args[1],reason_map,options.rankfile)


if __name__ == '__main__':
  main()
