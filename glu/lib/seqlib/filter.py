# -*- coding: utf-8 -*-

from __future__ import with_statement, division

__abstract__  = 'functions for alignment filtering'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


__all__ = ['filter_alignments','alignment_filter_options']

#             NAME            POSITIVE      NEGATIVE
FLAGS = { 'paired_end'      : (0x001,       0x000),
          'single_end'      : (0x000,       0x001),
          'proper_pair'     : (0x002|0x001, 0x000),
          'improper_pair'   : (      0x001, 0x002),
          'unmapped'        : (0x004,       0x000),
          'mapped'          : (0x000,       0x004),
          'mate_unmapped'   : (0x008|0x001, 0x000),
          'mate_mapped'     : (      0x001, 0x008),
          'reverse'         : (0x010,       0x000),
          'forward'         : (0x000,       0x010),
          'mate_reverse'    : (0x020|0x001, 0x000),
          'mate_forward'    : (      0x001, 0x020),
          'first_pair'      : (0x040|0x001, 0x000),
          'second_pair'     : (0x080|0x001, 0x000),
          'secondary_align' : (0x100,       0x000),
          'primary_align'   : (0x000,       0x100),
          'qc-'             : (0x200,       0x000),
          'qc+'             : (0x000,       0x200),
          'nondup'          : (0x000,       0x400),
          'dup'             : (0x400,       0x000),
        }


SUPPORTED_OPTIONS = set(FLAGS)
ALL_OPTIONS       = SUPPORTED_OPTIONS|set(['read_group','sample','library','center',
                                           'qname','has_seq','has_qual','query_len','seq_len'])


def normalize_options(opts):
  opts = opts or []

  if isinstance(opts,basestring):
    opts = [opts]

  opts = (o.strip().lower() for opt in opts for o in opt.split(','))
  return set(o for o in opts if o)


# Internal helper generator function for filter_alignments
def _filter_alignments(alignments, positive_flags, negative_flags, min_mapq):
  for alignment in alignments:
    if min_mapq and alignment.mapq < min_mapq:
      continue

    flags = alignment.flag
    pos_mismatch = (flags&positive_flags)!=positive_flags
    neg_mismatch = (flags&negative_flags)!=0

    if pos_mismatch or neg_mismatch:
      continue
      
    yield alignment


def filter_alignments(alignments, include, exclude, min_mapq=None):
  '''
  Filter alignments based on a set of specified options

  >>> class A(object):
  ...   def __init__(self,flag): self.flag = flag
  ...   def __repr__(self): return hex(self.flag)

  >>> aligns = [ A(0x000|0x004), A(0x001|0x200) ]
  >>> print list(filter_alignments(aligns, 'single_end',None))
  [0x4]
  >>> print list(filter_alignments(aligns, [],['single_end']))
  [0x201]
  >>> print list(filter_alignments(aligns, ['paired_end'],[]))
  [0x201]
  >>> print list(filter_alignments(aligns, [],['paired_end']))
  [0x4]
  >>> print list(filter_alignments(aligns, ['mapped'],[]))
  [0x201]
  >>> print list(filter_alignments(aligns, [],['mapped']))
  [0x4]
  >>> print list(filter_alignments(aligns, ['unmapped'],[]))
  [0x4]
  >>> print list(filter_alignments(aligns, [],['unmapped']))
  [0x201]
  '''
  include = normalize_options(include)
  exclude = normalize_options(exclude)

  unknown_options = (include|exclude)-ALL_OPTIONS

  if unknown_options:
    unknown_options = ', '.join(sorted(unknown_options))
    raise ValueError('Alignment filter does not recognize the following options: %s' % unknown_options)

  unsupported_options = (include|exclude)-SUPPORTED_OPTIONS

  if unsupported_options:
    unsupported_options = ', '.join(sorted(unsupported_options))
    raise ValueError('Alignment filter does not currently support the following options: %s' % unsupported_options)

  positive_flags = 0
  negative_flags = 0

  for flag in include:
    pos,neg = FLAGS[flag]
    positive_flags |= pos
    negative_flags |= neg

  for flag in exclude:
    neg,pos = FLAGS[flag]
    positive_flags |= pos
    negative_flags |= neg

  if not positive_flags and not negative_flags:
    return alignments

  return _filter_alignments(alignments,positive_flags,negative_flags)


def alignment_filter_options(group):
  opts = ', '.join(sorted(SUPPORTED_OPTIONS))

  group.add_argument('--includealign', action='append', metavar='OPTS',
                   help='Include alignments that meet the specified comma separated criteria. '
                        'These include: ' + opts)
  group.add_argument('--excludealign', action='append', metavar='OPTS',
                   help='Exclude alignments that meet the specified comma separated criteria. '
                        ' See --includealign for supported criteria.')
  group.add_argument('--minmapq', type=int, default=None, metavar='N', 
                  help='Exclude alignments with MAPQ < N')


def _test():
  import doctest
  doctest.testmod()


if __name__ == '__main__':
  _test()
