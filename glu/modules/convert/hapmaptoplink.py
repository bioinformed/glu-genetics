# -*- coding: utf-8 -*-
'''
File:          hapmaptoplink.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from itertools import dropwhile,islice,imap,chain

from glu.lib.fileutils import autofile,namefile,load_map


HAPMAP_HEADERS = ['rs# SNPalleles chrom pos strand genome_build center protLSID assayLSID panelLSID QC_code',
                  'rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode']


def load_hapmap_genotypes(file_or_name):
  '''
  Load the hampmap genotype data from file.
  @param file_or_name: file name or file object
  @type  file_or_name: str or file object
  @return: rows of genotypes with the first row being the sample names
  @rtype:  generator
  '''
  gfile = autofile(file_or_name)
  gfile = dropwhile(lambda s: s.startswith('#'), gfile)

  header = gfile.next()

  if not any(header.startswith(h) for h in HAPMAP_HEADERS):
    raise ValueError, "Input file '%s' does not appear to be in HapMap format." % namefile(file_or_name)

  header = list(islice(header.split(),11,None))
  n = len(header)

  yield header

  for line in gfile:
    fields = line.split()
    chr    = fields[2].strip('chr')
    rs     = fields[0]
    gd     = '%8.6f' % (float(fields[3])/1e8)
    pos    = fields[3]
    genos  = list(imap(format_geno, islice(fields,11,None)))
    assert len(genos) == n
    yield chain([chr,rs,gd,pos],genos)


def format_geno(genostr):
   return ' '.join(list(genostr.replace('N','0')))


def load_pedigree_file(pedfile,smapfile):
  '''
  @param smapfile: file providing a pedid_indid to sample(NA#) mapping
  '''
  smap  = load_map(smapfile,key_index=1,value_index=0)
  peds  = dict()
  pfile = autofile(pedfile)
  for line in pfile:
    fields = line.split()
    #FIXME: need to throw an exception if there is no sample map for an individual
    peds[smap.get(fields[0] + '_' + fields[1],None)] = fields
  return peds


def output_tped_file(rows,peds,tpedfile,tfamfile):
  tped = autofile(tpedfile,'w')
  tfam = autofile(tfamfile,'w')
  samples = rows.next()
  for sample in samples:
    tfam.write('%s\n' % ' '.join(peds[sample]))
  for row in rows:
    tped.write('%s\n' % ' '.join(row))


def option_parser():
  import optparse
  usage = 'Usage: %prog [options]  comparison...'

  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--genofile', dest='genofile',   metavar='FILE',
                     help='The file for genotype in HapMap format')
  parser.add_option('--pedfile',  dest='pedfile', metavar='FILE',
                     help='The pedigree file')
  parser.add_option('--smapfile', dest='smapfile',      metavar='FILE',
                     help='A map file providing map from the individual in pedfile to sample in genofile')
  parser.add_option('--tpedfile', dest='tpedfile',   metavar='FILE',
                     help='transposed ped file for plink')
  parser.add_option('--tfamfile', dest='tfamfile',   metavar='FILE',
                     help='transposed fam file for plink')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  genos = load_hapmap_genotypes(options.genofile)
  peds  = load_pedigree_file(options.pedfile,options.smapfile)

  output_tped_file(genos,peds,options.tpedfile,options.tfamfile)


if __name__ == '__main__':
  main()
