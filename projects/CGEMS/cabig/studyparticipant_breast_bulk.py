import csv
import sys

from biozilla.fileutils import autofile,load_map

AGE    = {'0':'>55','1':'55-59','2':'60-64','3':'65-69','4':'70-74','5':'>74','':''}
HEADER = ['Study Participant de-identifier ID','Gender','Age','Affection Status','Family History','Population']


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',        dest='outfile',        metavar='FILE',   default= '-',
                    help='Output file for bulk study_participant')
  parser.add_option('-i', '--infile',         dest='infile',         metavar='FILE',   default= '-',
                    help='the table data for study_participant')
  parser.add_option('-p', '--popmap',         dest='popmap',         metavar='FILE',   default= '-',
                    help='the table data for study_population')
  return parser


def main():

  parser = option_parser()
  options,args = parser.parse_args()

  popmap = load_map(options.popmap)
  #read the phenotyic file
  r = csv.reader(autofile(options.infile))
  r.next()
  out = csv.writer(autofile(options.outfile,'w'),dialect='excel-tab')
  out.writerow(HEADER)

  for row in r:
    pid     = row[1]
    age     = AGE[row[2]]
    gender  = row[4]
    status  = row[8]
    famhist = row[7]
    pop     = popmap[row[5]]
    out.writerow([pid,gender,age,status,famhist,pop])

if __name__=='__main__':
  main()
