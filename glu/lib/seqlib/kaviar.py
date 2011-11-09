# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   glu.lib.fileutils            import table_reader


def kaviar_reader(filename):
  for row in table_reader(filename):
    chrom = row[0]
    if not chrom.startswith('chr'):
      chrom = 'chr'+chrom

    loc     = int(row[1])
    mallele = row[3]
    maf     = float(row[4] or 0)

    for allelestuff in row[5:]:
      parts = allelestuff.split(':',1)
      if len(parts)!=2:
        continue

      allele,stuff = parts

      if allele=='rsids':
        continue

      stuff = stuff.strip().replace(', ',',')

      if not stuff.startswith('reference'):
        yield chrom,loc,mallele,maf,allele,stuff
