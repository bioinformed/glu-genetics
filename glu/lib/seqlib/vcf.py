# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Variant Call Format parser'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys

from   collections                  import defaultdict

from   glu.lib.fileutils            import table_reader
from   glu.lib.recordtype           import recordtype


VCFRecord = recordtype('VCFRecord',    'chrom start end names ref var qual filter info format genos')


class VCFReader(object):
  def __init__(self, filename, hyphen=None):
    self.data = data    = table_reader(filename,hyphen=sys.stdin)

    self.metadata_order = metadata_order = []
    self.metadata       = metadata       = defaultdict(list)
    self.header         = None
    self.samples        = None

    for row in data:
      if not row:
        continue
      elif row[0].startswith('##'):
        meta,value = row[0].split('=',1)
        meta = meta[2:]

        if meta not in metadata:
          metadata_order.append(meta)

        metadata[meta].append(row)

      elif row[0].startswith('#'):
        self.header = header = list(row)
        header[0]   = header[0][1:]

        self.samples = [ s.split('.')[0] for i,s in enumerate(header[9:]) ]
        break
      else:
        raise ValueError('Invalid VCF file detected')


  def __iter__(self):
    strcache    = {}
    def intern(x,strcache=strcache.setdefault): return strcache(x,x)

    for row in self.data:
      chrom      = intern(row[0])
      end        = int(row[1])
      start      = end-1
      names      = row[2].split(',') if row[2]!='.' else []
      ref        = intern(row[3])
      var        = [ intern(v) for v in row[4].split(',') ]
      qual       = row[5]
      filter     = [ intern(f) for f in row[6].split(';') ] if row[6]!='.' else []
      info       = row[7].split(';')
      format     = intern(row[8])
      genos      = [ g.split(':') for g in row[9:] ] if len(row)>9 else None

      yield VCFRecord(chrom,start,end,names,ref,var,qual,filter,info,format,genos)


