# -*- coding: utf-8 -*-

__abstract__  = 'WTCCC Raw genotype parser'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import csv

from   itertools                 import islice

from   glu.lib.utils             import is_str
from   glu.lib.fileutils         import autofile,namefile,parse_augmented_filename,get_arg

from   glu.lib.genolib.streams   import GenomatrixStream
from   glu.lib.genolib.reprs     import get_genorepr
from   glu.lib.genolib.locus     import Genome


__all__ = ['load_wtccc_raw']


def load_wtccc_raw(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a WTCCC Raw genotype data file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  genorepr = get_arg(args, ['genorepr']) or 'snp'
  unique   = get_arg(args, ['unique'], True)
  gc       = float(get_arg(args, ['gc'], 0))

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if is_str(genorepr):
    genorepr = get_genorepr(genorepr)

  gfile = autofile(filename)
  rows = csv.reader(gfile,dialect='tsv')

  try:
    columns = rows.next()
  except StopIteration:
    raise ValueError('Input file "%s" is empty' % namefile(filename))

  loci = [ intern(h) for h in islice(columns,1,None) ]

  if genome is None:
    genome = Genome()

  def _load_wtccc_raw(rows):
    n = len(columns)

    # Micro-optimization
    local_intern = intern
    local_strip  = str.strip

    gmap = {'-':'  '}

    for row in rows:
      if len(row) != n:
        raise ValueError('Invalid WTCCC raw row on line %d of %s' % (rows.line_num+1,namefile(filename)))

      sample = local_intern(local_strip(row[0]))
      data   = [ d.split(';') for d in islice(row,1,None) ]

      if gc:
        genos = [ d[0].replace('-',' ') if float(d[1])>=gc else '  ' for d in data ]
      else:
        genos = [ d[0].replace('-',' ') for d in data ]

      yield sample,genos

  genos = GenomatrixStream.from_strings(_load_wtccc_raw(rows),'sdat',genorepr=genorepr,loci=loci,
                                                                genome=genome,phenome=phenome,
                                                                unique=unique)

  if unique:
    genos = genos.unique_checked()

  return genos
