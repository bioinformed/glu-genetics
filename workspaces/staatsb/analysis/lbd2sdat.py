# -*- coding: utf-8 -*-
'''
File:          lbd2sdat.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-06-15

Abstract:      Parses an Illumina Locus by DNA report file and outputs a
               genotype matrix in sdat format.

Requires:      Python 2.5, biozilla

Revision:      $Id: lbd2sdat.py 335 2006-08-09 18:55:24Z staatsb $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv

from   itertools             import islice,chain,izip,imap,groupby
from   operator              import itemgetter

from   biozilla.utils        import autofile
from   biozilla.genodata     import save_genomatrix,filter_genomatrix_by_column,filter_genomatrix_by_row, \
                                    load_list
from   biozilla.genoarray    import snp_acgt


def load_lbd_file(filename, gcthreshold=None, locuslimit=None):
  data = csv.reader(autofile(filename), dialect='excel')

  def parseloci():
    row=islice(data.next(),3,None)
    if locuslimit:
      row = islice(row,locuslimit)
    return list(row)

  def skiprows(n):
    for i in xrange(n):
      data.next()

  def sample_generator():
    for key,rows in groupby(data,itemgetter(0,1,2,3,4)):
      rows = list(rows)
      assert len(rows) == 2
      genos,scores = rows
      assert  genos[6] == 'calls'
      assert scores[6] == 'Score_Call'

      sampleid = genos[0]

      if sampleid == 'blank':
        continue

      genos  = islice(genos,8,None)
      scores = imap(float,islice(scores,8,None))

      if locuslimit:
        genos = islice(genos, locuslimit) # limit genos for sample

      if gcthreshold:
        genos = [ (g if gc>gcthreshold else 'U') for g,gc in izip(genos,scores) ]

      yield sampleid,genos

  skiprows(13)
  strandhead  = parseloci()
  skiprows(1)
  locihead    = parseloci()
  skiprows(1)
  gencodehead = parseloci()
  skiprows(3)
  samples     = sample_generator()

  return izip(locihead,gencodehead,strandhead),samples


basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C',
                  'a': 't', 'c': 'g', 't': 't', 'g': 'c', '/': '/'}

def reversecomplement(seq):
  # TODO: generalize this to handle protein code, brackets,
  #       and move to a general genetics util package
  return ''.join(basecomplement[base] for base in reversed(seq))


def fix_strand(gencode, strand, targetstrand):
  '''Normalize genotypes to target strand'''
  if strand != targetstrand:
    gencode = reversecomplement(gencode)
  return gencode


def convert_ab_genos(samples, abmap, header):
  loci,lgenos = izip(*header)

  for sampleid,genos in samples:
    genos = (abmap[locus,geno] for locus,geno in izip(loci,genos))
    yield sampleid,snp_acgt.array(genos)


def build_abmap(header):
  abmap={}
  for h in header:
    locus,gencode = h
    a,b = gencode.split('/')
    abmap[locus,'A'] = a,a
    abmap[locus,'B'] = b,b
    abmap[locus,'H'] = a,b
    abmap[locus,'U'] = None
  return abmap


def option_parser():
  import optparse

  usage = 'usage: %prog [options] lbdfile...'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)
  parser.add_option('-o', '--outfile', dest='outfile', metavar='FILE', default='-',
                    help='Name of the outfile to be generated.')
  parser.add_option('-e', '--excludeloci', dest='excludeloci', metavar='FILE',
                    help='A file for excluding loci from the output')
  parser.add_option('-E', '--excludesamples', dest='excludesamples', metavar='FILE',
                    help='A file for excluding samples from the output')
  parser.add_option('-s', '--targetstrand', dest='targetstrand', metavar='N', type='string', default='T',
                    help='Strand the genotypes will be normalized to (T=Top, B=Bottom). Default=T')
  parser.add_option('-t', '--gcthreshold', dest='gcthreshold', type='float', metavar='N')
  parser.add_option('-l', '--samplelimit', dest='samplelimit', metavar='N', type='int', default=0,
                    help='Limit the number of samples considered to N for testing purposes (default=0 for unlimited)')
  parser.add_option('-L', '--locuslimit', dest='locuslimit', metavar='N', type='int', default=0,
                    help='Limit the number of loci considered to N for testing purposes (default=0 for unlimited)')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  args = iter(args)

  header,samples = load_lbd_file(args.next(),
                                 gcthreshold=options.gcthreshold,
                                 locuslimit=options.locuslimit)

  header = list(header)
  for arg in args:
    sub_header,more_samples = load_lbd_file(arg,
                                 gcthreshold=options.gcthreshold,
                                 locuslimit=options.locuslimit)

    if list(sub_header) != header:
      raise RuntimeError,'Genotype headers do not match'

    samples = chain(samples,more_samples)

  if options.samplelimit:
    samples = islice(samples,options.samplelimit)

  # normalize genotypes with respect to strand
  header = [ (locus,fix_strand(gencode,strand,options.targetstrand))
                for locus,gencode,strand in header ]

  # create mapping dictionary for locus
  abmap = build_abmap(header)

  samples = convert_ab_genos(samples, abmap, header)

  columns = [locus for locus,geno in header ]
  matrix = chain([columns],samples)

  if options.excludesamples:
    excludesamples = set(load_list(options.excludesamples))
    matrix = filter_genomatrix_by_row(matrix,excludesamples,exclude=True)

  if options.excludeloci:
    excludeloci = set(load_list(options.excludeloci))
    matrix = filter_genomatrix_by_column(matrix,excludeloci,exclude=True)

  # write sdat matrix file
  save_genomatrix(options.outfile,matrix)


if __name__ == '__main__':
  main()
