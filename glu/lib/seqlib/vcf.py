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


class VCFReader(object):
  def __init__(self, filename, hyphen=None, field_size_limit=1024*1024):
    if csv.field_size_limit()<field_size_limit:
      csv.field_size_limit(field_size_limit)

    self.data = data    = table_reader(filename,hyphen=sys.stdin)

    self.metadata       = metadata       = OrderedDict()
    self.header         = None
    self.samples        = None

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
    strcache    = {}
    def intern(x,strcache=strcache.setdefault): return strcache(x,x)

    for row in self.data:
      chrom      = intern(row[0])
      start      = int(row[1])-1
      names      = row[2].split(';') if row[2]!='.' else []
      ref        = intern(row[3])
      end        = start+len(ref)
      var        = [ intern(v) for v in row[4].split(',') ]
      qual       = row[5]
      filter     = [ intern(f) for f in row[6].split(';') ] if row[6]!='.' else []
      info       = row[7].split(';')
      format     = intern(row[8])

      if ':' in format:
        genos      = [ g.split(':') for g in row[9:] ] if len(row)>9 else None
      else:
        genos      = [ [g]          for g in row[9:] ] if len(row)>9 else None

      # VCF codes indels with an extra left reference base, which we strip
      r          = ref[0]
      if all(a.startswith(r) for a in var):
        start   += 1
        ref      = intern(ref[1:])
        var      = [ intern(v[1:]) for v in var ]

      yield VCFRecord(chrom,start,end,names,ref,var,qual,filter,info,format,genos)


class VCFWriter(object):
  def __init__(self, filename, metadata, names, reference):
    self.reference = pysam.Fastafile(reference)
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

    out = self.out
    out.write('\t'.join(map(str,row)))
    out.write('\n')
