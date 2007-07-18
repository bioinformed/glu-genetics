
import csv
import sys
from   biozilla.utils    import autofile
from   biozilla.genodata import load_map


CCOVERALL = {'1':0,'2':1,'3':2,'4':2}

def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output',      dest='output',      metavar='FILE',   default= '-',
                    help='Output file for table stdpt_analysis_grp_asso')
  parser.add_option('-p', '--pheno',       dest='pheno',       metavar='FILE',   default= '-',
                    help='Input pheno data file')
  parser.add_option('-m', '--partidmap',   dest='partidmap',   metavar='FILE',   default= '-',
                    help='Input file mapping study participant ids')

  return parser


def main():

  parser = option_parser()
  options,args = parser.parse_args()

  if not options.pheno:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify phenotype file'
    return

  if not options.partidmap:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify the participant id mapping file'
    return

  if not options.output:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify output file'
    return



  infile = csv.reader(autofile(options.pheno))
  infile.next()

  mapfile = csv.reader(autofile(options.partidmap),dialect='excel-tab')
  map = dict((row[0],row[1]) for row in mapfile)

  outfile = csv.writer(autofile(options.output,'w'))
  rows = set()
  for row in infile:
    if row[0] not in map:
      continue
    pid = map[row[0]]
    #stat column is 0 for control,1 for early, 2 for advanced
    anagrp = stat = int(row[22])
    #create a set of lists and then write out to eliminate repeated rows
    rows.add( (anagrp+13,pid) )
    rows.add( (anagrp+16,pid) )

    anagrp = ccoverall = CCOVERALL[row[5]]
    rows.add( (anagrp+19,pid) )
    rows.add( (anagrp+22,pid) )

  outfile.writerows(rows)


if __name__ == '__main__':
  main()
