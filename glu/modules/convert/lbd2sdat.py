# -*- coding: utf-8 -*-
'''
File:          lbd2sdat.py

Authors:       Brian Staats (staatsb@mail.nih.gov)
               Kevin Jacobs (jacobske@mail.nih.gov)

Created:       2006-06-15

Abstract:      Parses an Illumina Locus by DNA report file and outputs a
               genotype matrix or triple stream.

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   itertools            import islice,chain,izip,imap,groupby
from   operator             import itemgetter

from   glu.lib.utils        import autofile,hyphen
from   glu.lib.sections     import read_sections
from   glu.lib.sequence     import norm_snp_seq, complement_base

from   glu.lib.genolib.reprs     import get_genorepr
from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.io        import save_genostream
from   glu.lib.genolib.genoarray import model_from_alleles


def load_lbd_file(filename, gcthreshold=None):
  data = csv.reader(autofile(filename), dialect='excel')

  def parseloci():
    row = data.next()
    assert row[1] == 'locusNames'
    assert row[2] == ''
    return row[3:]

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

      if gcthreshold is not None:
        scores = imap(float,islice(scores,8,None))
        genos = [ (g if gc>gcthreshold else 'U') for g,gc in izip(genos,scores) ]

      yield sampleid,genos

  skiprows(15)
  loci    = parseloci()
  skiprows(5)
  samples = sample_generator()

  return loci,samples


def load_user_abmap(file_or_name,skip=0):
  '''
  Creates a dictionary representing a mapping from A and B probes for each locus to allele names.
  '''
  lfile = autofile(file_or_name)
  f = csv.reader(lfile, dialect='excel-tab')

  if skip:
    f = islice(f,skip,None)

  for row in f:
    if len(row) < 3 or '' in row[:3]:
      continue

    locus,a,b = row[:3]

    yield locus,(a,b)


def load_illumina_manifest(filename):
  ifile = csv.reader(autofile(filename),dialect='excel')
  sections = read_sections(ifile)

  heading,contents = sections.next()

  # OPA manifest
  if heading == 'data':
    heading,contents = sections.next()
    assert heading == 'Heading'

    headings = dict(c[:2] for c in islice(contents,10) if len(c) > 1)
    assert headings['Assay Format'] == 'Golden Gate'

  # Infinium manifest
  elif heading == 'Heading':
    contents = dict(c[:2] for c in contents if len(c) > 1)
    assert contents['Assay Format'] in ('Infinium','Infinium II')

    heading,contents = sections.next()
    assert heading == 'Assay'

  return contents


def find_index(header,headings):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass
  raise ValueError,'Cannot find heading index'


def extract_ab_from_manifest(manifest,targetstrand='customer'):
  targetstrand = targetstrand.lower()
  assert targetstrand in ('top','bottom','forward','reverse','customer','anticustomer','design','antidesign')

  manifest    = iter(manifest)
  header      = manifest.next()
  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  snp_idx     = find_index(header,['SNP'])
  cstrand_idx = find_index(header,['CustomerStrand'])
  dstrand_idx = find_index(header,['IlmnStrand','Ilmn Strand'])
  topseq_idx  = find_index(header,['TopGenomicSeq'])
  max_idx     = max(topseq_idx,dstrand_idx,cstrand_idx,snp_idx,name_idx,assayid_idx)

  for assay in manifest:
    assay  += ['']*(max_idx-len(assay)+1)
    locus   = assay[name_idx]
    snp     = assay[snp_idx]
    cstrand = assay[cstrand_idx].lower()
    dstrand = assay[dstrand_idx].lower()
    topseq  = assay[topseq_idx]

    assert cstrand in ('top','bot')
    assert dstrand in ('top','bot')
    assert (snp[0],snp[2],snp[4]) == ('[','/',']')

    # Alleles on the design strand
    aa,bb = snp[1],snp[3]

    try:
      tstrand,a,b = norm_snp_seq(topseq)
      tstrand     = tstrand.lower()
      assert tstrand == 'top'
    except ValueError:
      tstrand,a,b = dstrand,aa,bb

    if dstrand!=tstrand:
      aa = complement_base(aa)
      bb = complement_base(bb)

    if (a,b) != (aa,bb):
      raise ValueError('Sequence alleles do not match assay alleles')

    if targetstrand in ('forward','reverse'):
      # Get the strand orientation of the design sequence
      gstrand = assay[assayid_idx].split('_')[2]
      assert gstrand in 'FRU'
      if gstrand == 'U':
        raise ValueError("Unknown strand for assay '%s'" % locus)

      # Alleles are forward strand if the tstrand matches the design strand
      # and the design is on the forward strand or the converse of both
      # conditions is true.
      forward = (tstrand != dstrand) ^ (gstrand == 'F')

    flip =    ((targetstrand == 'customer'     and tstrand != cstrand)
           or  (targetstrand == 'anticustomer' and tstrand == cstrand)
           or  (targetstrand == 'design'       and tstrand != dstrand)
           or  (targetstrand == 'antidesign'   and tstrand == dstrand)
           or  (targetstrand == 'top'          and tstrand != 'top'  )
           or  (targetstrand == 'bottom'       and tstrand != 'bot'  )
           or  (targetstrand == 'forward'      and not forward       )
           or  (targetstrand == 'reverse'      and     forward       ))

    if flip:
      a,b = complement_base(a),complement_base(b)

    yield locus,(a,b)


def build_abmap(loci,usermap,modelmap):
  abmap={}
  for locus in loci:
    a,b = usermap.get(locus, ('A','B') )
    modelmap[locus]  = model_from_alleles([a,b])
    abmap[locus,'A'] = a,a
    abmap[locus,'B'] = b,b
    abmap[locus,'H'] = a,b
    abmap[locus,'U'] = None,None
  return abmap


def convert_ab_genos(loci, samples, abmap):
  for sampleid,genos in samples:
    genos = [abmap[locus,geno] for locus,geno in izip(loci,genos)]
    yield sampleid,genos


def option_parser():
  import optparse

  usage = 'usage: %prog [options] lbdfile...'
  parser = optparse.OptionParser(usage=usage)
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Name of the output file to be generated.')
  parser.add_option('-m', '--abmap', dest='abmap', metavar='FILE',
                    help='A file mapping A and B probes to alleles')
  parser.add_option('-M', '--manifest', dest='manifest', metavar='FILE',
                    help='Illumina manifest file')
  parser.add_option('-s', '--targetstrand', dest='targetstrand', metavar='T', default='customer',
                    help='Target strand based on Illumina manifest file: top, bottom, forward, reverse, customer (default), anticustomer, design, antidesign')
  parser.add_option('-F','--outformat',  dest='outformat', metavar='string',
                    help='The file output format for genotype data. Values=sdat, ldat, trip or genotriple')
  parser.add_option('-g', '--genorepr', dest='genorepr', metavar='REP', default='snp',
                    help='Input genotype representation.  Values=snp (default), hapmap, marker')
  parser.add_option('-t', '--gcthreshold', dest='gcthreshold', type='float', metavar='N')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  if options.output == '-' and not options.outformat:
    options.outformat = 'sdat'

  user_abmap = {}
  if options.manifest:
    print >> sys.stderr, 'Processing Illumina manifest file...',
    manifest = load_illumina_manifest(options.manifest)
    user_abmap.update(extract_ab_from_manifest(manifest,targetstrand=options.targetstrand))
    print >> sys.stderr, 'done.'

  if options.abmap:
    user_abmap.update(load_user_abmap(options.abmap))

  args = iter(args)
  loci,samples = load_lbd_file(args.next(),gcthreshold=options.gcthreshold)

  loci = list(loci)
  for arg in args:
    more_loci,more_samples = load_lbd_file(arg,gcthreshold=options.gcthreshold)

    if list(more_loci) != loci:
      raise RuntimeError,'Genotype headers do not match'

    samples = chain(samples,more_samples)

  modelmap = {}
  abmap    = build_abmap(loci,user_abmap,modelmap)
  samples  = convert_ab_genos(loci, samples, abmap)
  genos    = GenomatrixStream.from_tuples(samples, 'sdat', loci=loci, modelmap=modelmap)

  genorepr = get_genorepr(options.genorepr)
  save_genostream(hyphen(options.output,sys.stdout),genos,genorepr=genorepr)


if __name__ == '__main__':
  main()
