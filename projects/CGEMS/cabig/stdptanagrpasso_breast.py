
import csv
import sys
from   biozilla.fileutils    import autofile,load_list

HEADER    = ['ANA_GRP_ID','PARTICIPANT_DID']
CCOVERALL = {'1':0,'2':1,'3':2,'4':2}

def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output',      dest='output',      metavar='FILE',   default= '-',
                    help='Output file for table stdpt_analysis_grp_asso')
  parser.add_option('-p', '--pheno',       dest='pheno',       metavar='FILE',   default= '-',
                    help='Input pheno data file')
  parser.add_option('-m', '--partid',   dest='partid',   metavar='FILE',   default= '-',
                    help='Input file containing study participant ids')

  return parser


def main():

  parser = option_parser()
  options,args = parser.parse_args()

  if not options.pheno:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify phenotype file'
    return

  if not options.partid:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify the participant id file'
    return

  if not options.output:
    parser.print_help()
    print >> sys.stderr, 'Error: Must specify output file'
    return

  infile = csv.reader(autofile(options.pheno),dialect='excel-tab')
  infile.next()

  participants = load_list(options.partid)

  outfile = csv.writer(autofile(options.output,'w'))
  outfile.writerow(HEADER)
  rows = set()

  for row in infile:
    if row[0] not in participants:
      continue
    #stat column is 0 for control,1 for case
    anagrp = stat = int(row[1])
    #create a set of lists and then write out to eliminate repeated rows
    rows.add( (anagrp+25,row[0]) )
    rows.add( (anagrp+27,row[0]) )


  outfile.writerows(rows)


if __name__ == '__main__':
  main()
