# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Variant Call Format parser'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   collections                  import OrderedDict

import pysam

from   glu.lib                      import fileutils
from   glu.lib.recordtype           import recordtype


VCFRecordBase = recordtype('VCFRecord',    'chrom start end names ref var qual filter info format genostr genolist')


class VCFRecord(VCFRecordBase):
  __slots__ = ()

  @property
  def genos(self):
    genos = self.genolist

    if genos is not None:
      return genos

    format  = self.format
    genostr = self.genostr

    fields  = genostr.split('\t') if genostr else []

    if not fields:
      genos      = fields
    elif ':' in format:
      genos      = [ g.split(':') for g in fields ]
    else:
      genos      = [ [g]          for g in fields ]

    self.genolist = genos

    return genos


strcache    = {}
def _intern(x,strcache=strcache.setdefault): return strcache(x,x)


def _encode_vcf_records(data,normalize_indels=True):
  for line in data:
    line       = line.rstrip()
    fields     = line.split('\t',9)
    n          = len(fields)

    chrom      = _intern(fields[0])
    start      = int(fields[1])-1
    names      = fields[2].split(';') if fields[2]!='.' else []
    ref        = _intern(fields[3])
    end        = start+len(ref)
    var        = [ _intern(v) for v in fields[4].split(',') ]
    qual       = float(fields[5]) if fields[5]!='.' else None
    filter     = [ _intern(f) for f in fields[6].split(';') ] if fields[6]!='.' else []
    info       = fields[7].split(';') if fields[7]!='.' else []
    format     = _intern(fields[8]) if n>8 else None
    genostr    = fields[9] if n>9 else None

    # VCF codes indels with an extra left reference base, which we strip
    if normalize_indels:
      r          = ref[0]
      if all(a.startswith(r) for a in var):
        start   += 1
        ref      = _intern(ref[1:])
        var      = [ _intern(v[1:]) for v in var ]

    yield VCFRecord(chrom,start,end,names,ref,var,qual,filter,info,format,genostr,None)


class VCFReader(object):
  def __init__(self, filename, hyphen=sys.stdin, normalize_indels=True, field_size_limit=1024*1024):
    self.filename         = filename
    self.normalize_indels = normalize_indels
    self.tabixfile        = None
    self.data = data      = fileutils.autofile(fileutils.hyphen(filename,hyphen))

    self.metadata         = metadata       = OrderedDict()
    self.header           = None
    self.samples          = None

    for line in data:
      line = line.rstrip()

      if line.startswith('##'):
        meta,value = line.split('=',1)
        meta       = meta[2:]

        if meta not in metadata:
          metadata[meta] = []

        metadata[meta].append(line)

      elif line.startswith('#'):
        self.header  = header = line.split('\t')
        header[0]    = header[0][1:]
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

    records = tabixfile.fetch(chromosome, start, stop)

    return _encode_vcf_records(records,self.normalize_indels)


class VCFWriter(object):
  def __init__(self, filename, metadata, names, reference=None,hyphen=sys.stdout):
    self.reference = pysam.Fastafile(reference) if reference is not None else None
    self.out = out = fileutils.autofile(fileutils.hyphen(filename,hyphen),'w')

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
