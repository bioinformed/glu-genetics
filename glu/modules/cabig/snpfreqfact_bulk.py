# -*- coding: utf-8 -*-
'''
File:          snpfreqfact_bulk.py

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

from associationfinding_bulk import load_snp_dim,load_gene_snp_asso
from glu.lib.utils           import autofile

HEADER = ['dbSNP ID', 'Chromosome', 'Physical Position (bp)', 'Associated Genes', 'Population', 'Completion Rate(N/M)',
          'Hardy Weinberg pValue', 'Allele','Allele Count(Frequency)', 'Genotype', 'Genotype Count(Frequency)']

POPS = ['CASE','CASE_NONAGGRESSIVE','CASE_AGGRESSIVE','CONTROL','CEPH']

def load_study_pop(filename):
  r = csv.reader(autofile(filename))
  r.next()
  map = {}
  for row in r:
    id,name = row[1],row[2]
    map[id] = name

  return map


def main():
  snpinfo = load_snp_dim(sys.argv[1])
  genemap = load_gene_snp_asso(sys.argv[2])
  popmap = load_study_pop(sys.argv[3])
  freqs = csv.reader(autofile(sys.argv[4]))
  freqs.next()
  out = csv.writer(autofile(sys.argv[5],'w'),dialect='excel-tab')
  out.writerow(HEADER)

  popouts = {}
  for pop in POPS:
    w = csv.writer(autofile('frequencies_%s.txt.gz' % pop,'w'),dialect='excel-tab')
    w.writerow(HEADER)
    popouts[pop] = w

  for row in freqs:
    snpannoid = row[1]
    dbSNPid,chr,loc = snpinfo[snpannoid]

    #FIXME: consolidate with the paragraph in "associationfinding.py"
    genes = genemap.get(snpannoid,None)
    if genes is not None and len(genes) == 1:
      genes = genes[0]
    elif genes is not None and len(genes) > 1:
      genes = '|'.join(genes)

    pop = popmap[row[2]]

    refallele,refacount,refhomcount,otherallele,otheracount,otherhomcount,hetcount = row[3:10]
    missinggcount,hwp,completion = row[12:15]
    totalgcount = int(refhomcount) + int(otherhomcount) + int(hetcount)
    totalattempted = int(missinggcount) + totalgcount
    rhomfreq = float(refhomcount)/totalgcount
    ohomfreq = float(otherhomcount)/totalgcount
    hetfreq = float(hetcount)/totalgcount

    if totalgcount == totalattempted:
      completion = '100%% (%d/%d)' % (totalgcount,totalattempted)
    else:
      completion = '%.4f%% (%d/%d)' % (float(completion) * 100, totalgcount,totalattempted)

    hwp = '%.6f' % float(hwp)
    allele = refallele + '|' + otherallele
    acount = '%s(%.3f)|%s(%.3f)' % (refacount,float(refacount)/(2*totalgcount),otheracount,\
              float(otheracount)/(2*totalgcount))
    genotype = '%s%s|%s%s|%s%s' % (refallele,refallele,refallele,otherallele,otherallele,otherallele)
    gcount = '%s(%.3f)|%s(%.3f)|%s(%.3f)' % (refhomcount,rhomfreq,hetcount,hetfreq,otherhomcount,ohomfreq)

    out.writerow([dbSNPid,chr,loc,genes,pop,completion,hwp,allele,acount,genotype,gcount])
    (popouts[pop]).writerow([dbSNPid,chr,loc,genes,pop,completion,hwp,allele,acount,genotype,gcount])


if __name__=='__main__':
  main()

