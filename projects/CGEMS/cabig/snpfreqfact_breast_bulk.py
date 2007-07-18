import csv
import sys
from biozilla.fileutils      import autofile,load_map

HEADER = ['dbSNP ID', 'Chromosome', 'Physical Position (bp)', 'Associated Genes', 'Population', 'Completion Rate(N/M)',
          'Hardy Weinberg pValue', 'Allele','Allele Count(Frequency)', 'Genotype', 'Genotype Count(Frequency)']

POPS   = ['CASE','CONTROL']


def load_snp_dim(filename):
  r = csv.reader(autofile(filename))
  r.next()
  map = {}
  for row in r:
    snpannoid = row[7]
    chr       = row[0][3:]
    dbsnpid   = row[2]
    loc       = row[4]
    map[snpannoid] = [dbsnpid,chr,loc]

  return map


def load_gene_snp_asso(filename):
  r = csv.reader(autofile(filename),dialect='excel-tab')
  r.next()
  map = {}
  for row in r:
    snpannoid = row[2]
    map.setdefault(snpannoid,[]).append(row[1])

  return map


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--outfile',        dest='outfile',        metavar='FILE',   default= '-',
                    help='Output file for bulk snp_frequency_fact')
  parser.add_option('-s', '--snpmap',         dest='snpmap',         metavar='FILE',   default= '-',
                    help='the snp mapping file')
  parser.add_option('-g', '--genemap',        dest='genemap',        metavar='FILE',   default= '-',
                    help='the gene snp mapping file')
  parser.add_option('-p', '--popmap',         dest='popmap',         metavar='FILE',   default= '-',
                    help='the population mapping file')
  parser.add_option('-i', '--infile',         dest='infile',         metavar='FILE',   default= '-',
                    help='the data file for table snp_frequency_fact')

  return parser



def main():

  parser = option_parser()
  options,args = parser.parse_args()

  snpinfo = load_snp_dim(options.snpmap)
  genemap = load_gene_snp_asso(options.genemap)
  popmap  = load_map(options.popmap)
  freqs   = csv.reader(autofile(options.infile))
  freqs.next()
  out = csv.writer(autofile(options.outfile,'w'),dialect='excel-tab')
  out.writerow(HEADER)

  popouts = {}
  for pop in POPS:
    w = csv.writer(autofile('frequencies_%s.txt.gz' % pop,'w'),dialect='excel-tab')
    w.writerow(HEADER)
    popouts[pop] = w

  for row in freqs:
    snpannoid = row[1]
    dbSNPid,chr,loc = snpinfo[snpannoid]

    #FIXME: consolidate with the paragraph in "associationfinding.py"
    genes = genemap.get(snpannoid,None)
    if genes is not None and len(genes) == 1:
      genes = genes[0]
    elif genes is not None and len(genes) > 1:
      genes = '|'.join(genes)

    pop = popmap[row[2]]

    refallele,refacount,refhomcount,otherallele,otheracount,otherhomcount,hetcount = row[3:10]
    missinggcount,hwp,completion = row[12:15]
    totalgcount    = int(refhomcount) + int(otherhomcount) + int(hetcount)
    totalattempted = int(missinggcount) + totalgcount
    rhomfreq = float(refhomcount)/totalgcount
    ohomfreq = float(otherhomcount)/totalgcount
    hetfreq  = float(hetcount)/totalgcount

    if totalgcount == totalattempted:
      completion = '100%% (%d/%d)' % (totalgcount,totalattempted)
    else:
      completion = '%.4f%% (%d/%d)' % (float(completion) * 100, totalgcount,totalattempted)

    hwp = '%.6f' % float(hwp)
    allele = refallele + '|' + otherallele
    acount = '%s(%.3f)|%s(%.3f)' % (refacount,float(refacount)/(2*totalgcount),otheracount,\
              float(otheracount)/(2*totalgcount))
    genotype = '%s%s|%s%s|%s%s' % (refallele,refallele,refallele,otherallele,otherallele,otherallele)
    gcount = '%s(%.3f)|%s(%.3f)|%s(%.3f)' % (refhomcount,rhomfreq,hetcount,hetfreq,otherhomcount,ohomfreq)

    out.writerow([dbSNPid,chr,loc,genes,pop,completion,hwp,allele,acount,genotype,gcount])
    (popouts[pop]).writerow([dbSNPid,chr,loc,genes,pop,completion,hwp,allele,acount,genotype,gcount])


if __name__=='__main__':
  main()

