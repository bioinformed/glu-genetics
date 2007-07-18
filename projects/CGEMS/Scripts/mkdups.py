__version__ = '0.1'

import sys
import csv

from   itertools  import islice, izip
from   textwrap   import fill

from   biozilla.utils      import any, all, autofile, percent


def load_smap(filename):
  samples = csv.reader(autofile(filename),dialect='excel-tab')
  return dict( (s,i) for s,i in samples )

def load_samples(filename):
  return (row[0] for row in csv.reader(autofile(filename),dialect='excel-tab'))


def option_parser():
  import optparse

  usage = 'usage: %prog [options] samples smap'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 2:
    parser.print_help()
    return

  samples = load_samples(args[0])
  smap    = load_smap(args[1])

  dups = {}

  for sample in samples:
    sb = sample[-8:]
    exemplar = smap.get(sb,sb)
    dups.setdefault(exemplar,set()).add(sample)

  for dupset in dups.itervalues():
    if len(dupset) > 1:
      print '\t'.join(dupset)


if __name__ == '__main__':
  main()
