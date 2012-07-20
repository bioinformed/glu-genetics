# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Variant Call Format parser'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import sys

from   collections                  import defaultdict, OrderedDict

import pysam

from   glu.lib.fileutils            import table_reader, autofile, hyphen
from   glu.lib.recordtype           import recordtype


VCFRecord = recordtype('VCFRecord',    'chrom start end names ref var qual filter info format genos')


strcache    = {}
def _intern(x,strcache=strcache.setdefault): return strcache(x,x)


def _encode_vcf_records(rows,normalize_indels=True):
  for row in rows:
    chrom      = _intern(row[0])
    start      = int(row[1])-1
    names      = row[2].split(';') if row[2]!='.' else []
    ref        = _intern(row[3])
    end        = start+len(ref)
    var        = [ _intern(v) for v in row[4].split(',') ]
    qual       = row[5]
    filter     = [ _intern(f) for f in row[6].split(';') ] if row[6]!='.' else []
    info       = row[7].split(';') if row[7]!='.' else []

    n          = len(row)
    if n>8:
      format     = _intern(row[8])

      if n>9:
        if ':' in format:
          genos      = [ g.split(':') for g in row[9:] ] if len(row)>9 else None
        else:
          genos      = [ [g]          for g in row[9:] ] if len(row)>9 else None
      else:
        genos = None
    else:
      format = genos = None

    # VCF codes indels with an extra left reference base, which we strip
    if normalize_indels:
      r          = ref[0]
      if all(a.startswith(r) for a in var):
        start   += 1
        ref      = _intern(ref[1:])
        var      = [ _intern(v[1:]) for v in var ]

    yield VCFRecord(chrom,start,end,names,ref,var,qual,filter,info,format,genos)


class VCFReader(object):
  def __init__(self, filename, hyphen=None, normalize_indels=True, field_size_limit=1024*1024):
    if csv.field_size_limit()<field_size_limit:
      csv.field_size_limit(field_size_limit)

    self.filename         = filename
    self.normalize_indels = normalize_indels
    self.tabixfile        = None
    self.data = data      = table_reader(filename,hyphen=sys.stdin)

    self.metadata         = metadata       = OrderedDict()
    self.header           = None
    self.samples          = None

    for row in data:
      if not row:
        continue
      elif row[0].startswith('##'):
        if len(row)!=1:
          raise ValueError('Invalid VCF header line')

        meta,value = row[0].split('=',1)
        meta = meta[2:]

        if meta not in metadata:
          metadata[meta] = []

        metadata[meta].append(row[0])

      elif row[0].startswith('#'):
        self.header = header = list(row)
        header[0]   = header[0][1:]

        self.samples = [ s.split('.')[0] for s in header[9:] ]
        break
      else:
        raise ValueError('Invalid VCF file detected')


  def __iter__(self):
    return _encode_vcf_records(self.data,self.normalize_indels)

  def fetch(self, chromosome, start, stop):
    tabixfile = self.tabixfile

    if tabixfile is None:
      self.tabixfile = tabixfile = pysam.Tabixfile(self.filename,cache_size=128*1024*1024)

    records = [ r.split('\t') for r in tabixfile.fetch(chromosome, start, stop) ]

    return _encode_vcf_records(records,self.normalize_indels)


class VCFWriter(object):
  def __init__(self, filename, metadata, names, reference=None):
    self.reference = pysam.Fastafile(reference) if reference is not None else None
    self.out = out = autofile(hyphen(filename,sys.stdout),'w')

    for meta in metadata:
      for m in metadata[meta]:
        out.write(m)
        out.write('\n')

    out.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+names))
    out.write('\n')

  def write_locus(self, vcfvar):
    # FORMAT: chrom start end names ref var filter info format genos
    pos = vcfvar.start+1
    ref = vcfvar.ref
    var = vcfvar.var

    # VCF codes indels with an extra left reference base
    if not vcfvar.ref or '' in vcfvar.var:
      if self.reference is None:
        raise RuntimeError('Reference sequence required to write VCF indel loci')

      pos -= 1
      r    = self.reference.fetch(vcfvar.chrom,pos-1,pos).upper()
      ref  = r+ref
      var  = [ r+a for a in var ]

    row = [ vcfvar.chrom,
            str(pos),
            ';'.join(vcfvar.names) or '.',
            ref,
            ','.join(var),
            vcfvar.qual,
            ';'.join(sorted(vcfvar.filter)) or '.',
            ';'.join(vcfvar.info),
            vcfvar.format ] + [ ':'.join(g) for g in vcfvar.genos ]

    out  = self.out

    text = '\t'.join(map(str,row))

    if ' ' in text:
      text = text.replace(' ','_')

    out.write(text)
    out.write('\n')
