# -*- coding: utf-8 -*-
'''
File:          splitgenos.py

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


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-m', '--snpchr',   dest='snpchr',   metavar='FILE',   default= '-',
                    help='Input file to map between snp id and chr')
  parser.add_option('-M', '--sampop',   dest='sampop',   metavar='FILE',   default= '-',
                    help='Input file to map between sample id and pop id')
  parser.add_option('-g', '--genos',    dest='genos',    metavar='FILE',   default= '-',
                    help='Input genotype file')
  parser.add_option('-d', '--dir',      dest='dir',      type = 'str',
                    help='Output file directory')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  chrs = map(str,range(1,23))
  chrs.extend(['X','XY','Y','M'])

  pops = ['CONTROL','CASE_NONAGGRESSIVE','CASE_AGGRESSIVE','CEPH']

  snpchr = load_map(options.snpchr,0,1,1)
  sampop = load_map(options.sampop,12,8,1)
  for key,value in sampop.iteritems():
    if value == 'QC':
      sampop[key] = 'CEPH'

  #open the file handles and store them in a dictionary
  filedict = {}
  for chr in chrs:
    filename = '%s/genotypes_chr%s.txt.gz' % (options.dir,chr)
    w = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
    filedict[chr] = w

    for pop in pops:
      if pop not in filedict:
        filename = '%s/genotypes_%s.txt.gz' % (options.dir,pop)
        w = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
        filedict[pop] = w

    filename = '%s/genotypes_%s_chr%s.txt.gz' % (options.dir,pop,chr)
    w = csv.writer(autofile(filename, 'w'),dialect='excel-tab')
    filedict[(pop,chr)] = w

  genos = csv.reader(autofile(options.genos),dialect='excel-tab')

  header = genos.next()
  for w in filedict.itervalues():
    w.writerow(header)

  for row in genos:
    specid = row[2]
    snpid  = row[0]
    chr    = snpchr[snpid]
    pop    = sampop[specid]

    w = filedict[chr]
    w.writerow(row)

    w = filedict[pop]
    w.writerow(row)

    w = filedict[(pop,chr)]
    w.writerow(row)


if __name__=='__main__':
  main()

