# -*- coding: utf-8 -*-

from __future__ import division


import pysam

from   glu.lib.fileutils    import list_reader


class ReferenceVariants(object):
  def __init__(self, refvariant, ingroupfile=None):
    self.vars      = pysam.Tabixfile(refvariant,cache_size=128*1024*1024)
    self.ingroup   = set(s.strip() for s in list_reader(ingroupfile)) if ingroupfile else ()

  def get(self, chromosome, start, stop, var):
    var     = set(var)
    records = self.vars.fetch(chromosome, start, stop+1)

    i = set()
    o = set()

    for row in records:
      vchrom,vstart,vstop,vgeno,vsubjects = row.strip().split('\t')
      if vchrom!=chromosome or int(vstart)!=start or int(vstop)!=stop:
        continue

      vvars = set(vgeno.split('/')[1:])

      if not vvars&var:
        continue

      ingroup = self.ingroup
      vsubjects = [ s.strip() for s in vsubjects.split(',') ]
      for s in vsubjects:
        if s in ingroup:
          i.add(s)
        else:
          o.add(s)

    i = sorted(i) if i else []
    o = sorted(o) if o else []

    return i,o
