
import csv
import sys

from operator              import itemgetter
from itertools             import islice,count
from biozilla.fileutils    import autofile,load_map


HEADER_OUT_P = ['ASSO_ANA_FACT_ID','ASSO_ANA_PVALUE','ASSO_ANA_RANK','ASSO_ANALYSIS_ID','SNPANNO_ID','OR_NONAGGRESSIVE_HETEROZYGOTE','OR_NONAGGRESSIVE_HOMOZYGOTE','OR_AGGRESSIVE_HETEROZYGOTE','OR_AGGRESSIVE_HOMOZYGOTE']

HEADER_OUT_B = ['ASSO_ANA_FACT_ID','ASSO_ANA_PVALUE','ASSO_ANA_RANK','ASSO_ANALYSIS_ID','SNPANNO_ID','OR_HETEROZYGOTE','OR_HOMOZYGOTE']

HEADER_IN_B = ['Rank','Locus','Unadjusted','p-value','HetOR','HomOR','Adjusted','p-value','HetOR','HomOR','CHROMOSOME','LOCATION','GENE NEIGHBORHOOD']

HEADER_IN_P = ['Locus','Contingency Table','p-value','HetOR1','HomOR1','HetOR2','HomOR2','Score with covariates','p-value','HetOR1','HomOR1','HetOR2','HomOR2','DF']

def write_rows(infile,out,skip,ids,snpmap,id1,id2):
  data = list(islice(csv.reader(autofile(infile),dialect='excel-tab'),skip,None))

  data.sort(key=itemgetter(9))
  # append the rank for the second p value to the original data list
  r = 0
  for d in data:
    if d[9] and d[9].find('nan') == -1:
      d[9] = '%12.6f' % float(d[9])
      r += 1
      d.append(r)
    else:
      d[9] = ''
      d.append('')

  data.sort(key=itemgetter(4))
  r = 0
  for d in data:
    if d[4] and d[4].find('nan') == -1:
      r += 1
      pvalue1 = '%12.6f' % float(d[4])
      rank1 = r
    else:
      pvalue1 = ''
      rank1 = ''

    pvalue2 = d[9]
    rank2   = d[15]

    rows= [[ids.next(),pvalue1,rank1,id1,snpmap[d[1]],d[5],d[6]],[ids.next(),pvalue2,rank2,id2,snpmap[d[1]],d[10],d[11]]]
    out.writerows(rows)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',   dest='outfile',   metavar='FILE',   default= '-',
                    help='Output file for table snp_analysis_finding_fact')
  parser.add_option('-i', '--ID',        dest='ID',        metavar='FILE',   default= '-',
                    help='Input file of association results from incidence density sampling')
  parser.add_option('-I', '--SS',        dest='SS',        metavar='FILE',   default= '-',
                    help='Input file of association results from selection sampling')
  parser.add_option('-m', '--snpmap',    dest='snpmap',    metavar='FILE',   default= '-',
                    help='Input file mapping dbsnp ids')
  parser.add_option('-c', '--inicount',  dest='inicount',  type = 'int')


  return parser



def main():

  parser = option_parser()
  options,args = parser.parse_args()

  snpmap = load_map(options.snpmap)

  #use a sequence number starting from the last one in the database for the asso_ana_fact_id
  ids = count(options.inicount)

  out = csv.writer(autofile(options.outfile,'w'))
  out.writerow(HEADER_OUT_B)

  write_rows(options.ID,out,1,ids,snpmap,9,10)
  #write_rows(options.SS,out,0,ids,snpmap,7,8)


if __name__=='__main__':
  main()
