import csv
import sys

from operator              import itemgetter
from itertools             import islice, chain

from biozilla.fileutils    import autofile, hyphen
from biozilla.xtab         import xtab_list


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',   dest='outfile',   metavar='FILE',   default= '-',
                    help='Output file for formatted data')
  parser.add_option('-m', '--measure',  dest='measure', default='r2',
                    help="Measure of LD: r2 (default) or D'")

  return parser


def merge(i,j,x):
  if x:
    return x[0]
  else:
    return ''


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    parser.print_help()
    return

  out = csv.writer(autofile(hyphen(options.outfile,sys.stdout),'w'),dialect='excel-tab')

  datain = [ islice(csv.reader(autofile(hyphen(arg,sys.stdin)),dialect='excel-tab'),1,None) for arg in args ]
  datain = chain(*datain)

  if options.measure.lower() == 'r2':
    col = 2
  elif options.measure.lower() == "d'":
    col = 3
  else:
    raise ValueError, 'Unknown or unsupported LD measure specified: %s' % options.measure

  rows = xtab_list(datain,itemgetter(1),itemgetter(0),itemgetter(col),merge)
  out.writerow(['']+rows.next())
  for label,row in rows:
    out.writerow([label]+row)

if __name__=='__main__':
  main()
