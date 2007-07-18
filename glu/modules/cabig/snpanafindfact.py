# -*- coding: utf-8 -*-
'''
File:          snpanafindfact.py

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

from   operator          import itemgetter
from   itertools         import islice,count

from   glu.lib.fileutils import autofile,load_map


HEADER = ['ASSO_ANA_FACT_ID','ASSO_ANA_PVALUE','ASSO_ANA_RANK','ASSO_ANALYSIS_ID','SNPANNO_ID','OR_NONAGGRESSIVE_HETEROZYGOTE','OR_NONAGGRESSIVE_HOMOZYGOTE','OR_AGGRESSIVE_HETEROZYGOTE','OR_AGGRESSIVE_HOMOZYGOTE']



def write_rows(infile,out,skip,ids,snpmap,id1,id2):
  data = list(islice(csv.reader(autofile(infile),dialect='excel-tab'),skip,None))

  data.sort(key=itemgetter(8))
  # append the rank for the second p value to the original data list
  r = 0
  for d in data:
    if d[8] and d[8].find('nan') == -1:
      d[8] = '%12.6f' % float(d[8])
      r += 1
      d.append(r)
    else:
      d[8] = ''
      d.append('')

  data.sort(key=itemgetter(2))
  r = 0
  for d in data:
    if d[2] and d[2].find('nan') == -1:
      r += 1
      pvalue1 = '%12.6f' % float(d[2])
      rank1 = r
    else:
      pvalue1 = ''
      rank1 = ''

    pvalue2 = d[8]
    rank2   = d[14]

    rows= [[ids.next(),pvalue1,rank1,id1,snpmap[d[0]],d[4],d[3],d[6],d[5]],[ids.next(),pvalue2,rank2,id2,snpmap[d[0]],d[10],d[9],d[12],d[11]]]
    out.writerows(rows)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',   dest='outfile',   metavar='FILE',   default= '-',
                    help='Output file for table snp_analysis_finding_fact')
  parser.add_option('-i', '--ID',        dest='ID',        metavar='FILE',   default= '-',
                    help='Input file of association results from incidence density sampling')
  parser.add_option('-I', '--SS',        dest='SS',        metavar='FILE',   default= '-',
                    help='Input file of association results from selection sampling')
  parser.add_option('-m', '--snpmap',    dest='snpmap',    metavar='FILE',   default= '-',
                    help='Input file mapping dbsnp ids')
  parser.add_option('-c', '--inicount',  dest='inicount',  type = 'int')


  return parser



def main():

  parser = option_parser()
  options,args = parser.parse_args()

  snpmap = load_map(options.snpmap)

  #use a sequence number starting from the last one in the database for the asso_ana_fact_id
  ids = count(options.inicount)

  out = csv.writer(autofile(options.outfile,'w'))
  out.writerow(HEADER)

  write_rows(options.ID,out,1,ids,snpmap,5,6)
  write_rows(options.SS,out,0,ids,snpmap,7,8)


if __name__=='__main__':
  main()
