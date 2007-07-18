# -*- coding: utf-8 -*-
'''
File:          associationfinding_breast_bulk.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from itertools         import chain

from glu.lib.fileutils import autofile,load_map

def load_snp_dim(filename):
  r = csv.reader(autofile(filename))
  r.next()
  map = {}
  for row in r:
    snpannoid = row[7]
    chr = row[0][3:]
    dbsnpid = row[2]
    loc = row[4]
    map[snpannoid] = [dbsnpid,chr,loc]

  return map


def load_gene_snp_asso(filename):
  r = csv.reader(autofile(filename),dialect='excel-tab')
  r.next()
  map = {}
  for row in r:
    snpannoid = row[2]
    map.setdefault(snpannoid,[]).append(row[1])

  return map


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',        dest='outfile',        metavar='FILE',   default= '-',
                    help='Output file for bulk snp_analysis_finding_fact')
  parser.add_option('-s', '--snpmap',         dest='snpmap',         metavar='FILE',   default= '-',
                    help='the snp mapping file')
  parser.add_option('-g', '--genemap',        dest='genemap',        metavar='FILE',   default= '-',
                    help='the gene snp mapping file')
  parser.add_option('-a', '--analysis',       dest='analysis',       metavar='FILE',   default= '-',
                    help='the data for table snp_association_analysis')
  parser.add_option('-i', '--infile',         dest='infile',         metavar='FILE',   default= '-',
                    help='the data file for table snp_analysis_finding_fact')

  return parser


def main():

  parser = option_parser()
  options,args = parser.parse_args()


  snpmap  = load_snp_dim(options.snpmap)
  genemap = load_gene_snp_asso(options.genemap)
  anamap  = load_map(options.analysis)
  associations = csv.reader(autofile(options.infile))
  associations.next()
  out = csv.writer(autofile(options.outfile,'w'),dialect='excel-tab')
  out.writerow(['dbSNP ID', 'Chromosome', 'Physical Position (bp)', 'Associated Genes', 'Analysis Name', 'p-value',
                'Whole Genome Rank','OR Heterozygote','OR Homozygote'])
  for row in associations:
    snpannoid = row[4]
    snpinfo = snpmap.get(snpannoid)
    genes = genemap.get(snpannoid,None)
    name = anamap[row[3]]
    if genes is not None and len(genes) == 1:
      genes = genes[0]
    elif genes is not None and len(genes) > 1:
      genes = '|'.join(genes)
    out.writerow(list(chain(snpinfo,[genes,name,row[1],row[2]],row[5:7])))

if __name__ == '__main__':
  main()
