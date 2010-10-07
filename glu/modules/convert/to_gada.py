# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Convert an Illumina Locus by DNA report file into a GLU genotype file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys

from   itertools                 import izip

import h5py

from   glu.lib.illumina          import create_Illumina_abmap

from   glu.lib.fileutils         import table_writer

from   glu.lib.genolib.locus     import Genome
from   glu.lib.genolib.transform import GenoTransform


def split_fullname(filename):
  dirname  = os.path.dirname(filename)

  # Get filename
  filename = os.path.basename(filename)

  # Split filename into 1 or 2 parts up to the first '.'
  parts = filename.split('.',1)

  # Combine dirname and up to the first '.' of filename as prefix
  prefix = os.path.join(dirname,parts[0])

  # Suffix the remainder of filename after the first '.'
  suffix = '' if len(parts) == 1 else parts[1]

  return prefix,suffix


def option_parser():
  import optparse

  usage = 'usage: %prog [options] gdat'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('--includesamples', dest='includesamples', metavar='FILE', action='append',
                    help='List of samples to include, all others will be skipped')
  parser.add_option('--excludesamples', dest='excludesamples', metavar='FILE', action='append',
                    help='List of samples to exclude, only samples not present will be kept')
  parser.add_option('--renamesamples', dest='renamesamples', metavar='FILE',
                    help='Rename samples from a file containing rows of original name, tab, new name')
  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output genotype file name')
  parser.add_option('-M', '--manifest', dest='manifest', metavar='FILE',
                    help='Illumina manifest file (BPM or CSV)')
  parser.add_option('-w', '--warnings', action='store_true', dest='warnings',
                    help='Emit warnings and A/B calls for SNPs with invalid manifest data')
  return parser


def create_recode(loci,genome,abmap):
  empty      = {'DD':'AA', 'ID':'AB', 'DI':'AB', 'II':'BB','  ':'NC'}
  geno_cache = {}
  recode     = []
  locusmap   = genome.loci

  for lname in loci:
    loc     = locusmap[lname]
    alleles = abmap.get(lname)

    genomap = geno_cache.get(alleles)
    if not genomap:
      if alleles:
        a,b = alleles
        aa  = a+a
        ab  = a+b
        ba  = b+a
        bb  = b+b
        genomap = { aa:'AA', ab:'AB', ba:'AB', bb:'BB', '  ':'NC' }
        if alleles==('-','+') or alleles==('+','-'):
          genomap.update(empty)

      else:
        genomap = empty

      geno_cache[alleles] = genomap

    recode.append( (loc.name,loc.chromosome,loc.location,genomap) )

  return recode


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args)!=1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  if not options.manifest:
    sys.stderr.write('ERROR: An Illumina manifest file must be specified')
    sys.exit(2)

  if not options.output:
    sys.stderr.write('ERROR: An output filename template must be specified')
    sys.exit(2)

  prefix,suffix = split_fullname(options.output)

  transform = GenoTransform.from_object(options)

  genome = Genome()

  errorhandler = None
  if options.warnings:
    def errorhandler(msg):
      sys.stderr.write('WARNING: %s\n' % msg)

  sys.stderr.write('Processing Illumina manifest file...')
  abmap = create_Illumina_abmap(options.manifest,genome,targetstrand='forward',
                                                        errorhandler=errorhandler)
  sys.stderr.write('done.\n')

  indata = h5py.File(args[0])

  snps    = indata['SNPs'][:]
  samples = indata['Samples'][:]
  genos   = indata['Genotype']
  lrr     = indata['LRR']
  baf     = indata['BAF']

  recode  = create_recode(snps,genome,abmap)
  include = transform.samples.include
  exclude = transform.samples.exclude
  rename  = transform.samples.rename

  for i in xrange(len(samples)):
    sample = samples[i]

    if include is not None and sample not in include:
      continue

    if exclude is not None and sample in exclude:
      continue

    if rename is not None:
      sample = rename.get(sample,sample)

    out = table_writer('%s_%s.%s' % (prefix,sample,suffix), hyphen=sys.stdout)
    out.writerow( ('Name','Chr','Position','Log.R.Ratio','B.Allele.Freq','GType') )

    sample_data = izip(recode,lrr[i]/10000.,baf[i]/10000.,genos[i])
    for (lname,chrom,location,genomap),l,b,geno in sample_data:
      out.writerow( (lname,chrom,location,l,b,genomap[geno]) )


if __name__ == '__main__':
  main()
