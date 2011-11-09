# -*- coding: utf-8 -*-

from __future__ import division


from   collections          import defaultdict, namedtuple

import pysam

from   glu.lib.recordtype   import recordtype


VarInfo    = namedtuple('VarInfo',   'exact_vars inexact_vars common_score function_score')
VarRecord  = recordtype('VarRecord', 'chromosome start stop allele common_score function_score source')
VariantKey = namedtuple('VariantKey', 'chromosome start stop source')


class CGFVariants(object):
  def __init__(self, cgfvariant, reference_fasta):
    self.vars      = pysam.Tabixfile(cgfvariant,cache_size=128*1024*1024)
    self.reference = pysam.Fastafile(reference_fasta)

  def query_variants(self, chromosome, start, stop):
    if chromosome.startswith('chr'):
      chromosome = chromosome[3:]

    chrmap = {'X':23,'Y':24,'MT':25,'M':25}
    score  = self.vars.fetch(chrmap.get(chromosome,chromosome), start, stop)

    for s in score:
      chrom,vstart,vstop,allele,common_score,function_score,source = s.split('\t')
      source = [ s.strip() for s in source.split(',') ]

      yield VarRecord(chrom, int(vstart), int(vstop), allele,
                             float(common_score), int(function_score),
                             source)

  def build_variants_lookup(self, chromosome, start, stop):
    vdict = defaultdict(list)
    vdata = self.query_variants(chromosome, start, stop)

    for d in vdata:
      for source in d.source:
        vdict[VariantKey(d.chromosome,d.start,d.stop,source)].append(d)

    return vdict

  def get_refseq(self, chromosome, start, stop):
    if not chromosome.startswith('chr'):
      chromosome='chr'+chromosome
    return self.reference.fetch(chromosome, start, stop).upper()

  def score_and_classify(self, chromosome, start, stop, geno):
    geno           = [ g.strip().upper() for g in geno ]
    qref           = self.get_refseq(chromosome, start, stop)
    vdata          = self.build_variants_lookup(chromosome, start, stop)
    qvar_alleles   = set(a for a in geno if a!=qref)

    exact_vars     = []
    inexact_vars   = []
    common_score   = defaultdict(float)
    function_score = 0

    for v,vrecs in vdata.iteritems():
      if v.start!=start or v.stop!=stop:
        if v.source:
          inexact_vars.append(v.source)
        continue

      valleles = set(va.allele for va in vrecs)
      for vrec in vrecs:
        a = vrec.allele
        if a in geno:
          common_score[a] = max(common_score[a],vrec.common_score)
          function_score |= vrec.function_score

        if v.source:
          if qvar_alleles-valleles:
            inexact_vars.append(v.source)
          else:
            exact_vars.append(v.source)

    exact_vars   = sorted(set(exact_vars))
    inexact_vars = sorted(set(inexact_vars))
    common_score = min(common_score.values()) if common_score else 0

    return VarInfo(exact_vars,inexact_vars,common_score,function_score)
