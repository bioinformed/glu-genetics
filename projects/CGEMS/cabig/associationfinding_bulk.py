import sys
import csv
from itertools      import chain

from biozilla.utils import autofile
from utils          import load_map

def load_snp_dim(filename):
  r = csv.reader(autofile(filename))
  r.next()
  map = {}
  for row in r:
    snpannoid = row[7]
    chr = row[0][3:]
    dbsnpid = row[2]
    loc = row[4]
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


def main():
  snpmap = load_snp_dim(sys.argv[1])
  genemap = load_gene_snp_asso(sys.argv[2])
  anamap = load_map(sys.argv[3],0,1,1,dialect='excel')
  associations = csv.reader(autofile(sys.argv[4]))
  associations.next()
  out = csv.writer(autofile(sys.argv[5],'w'),dialect='excel-tab')
  out.writerow(['dbSNP ID', 'Chromosome', 'Physical Position (bp)', 'Associated Genes', 'Analysis Name', 'p-value',
                'Whole Genome Rank','OR Nonaggressive Heterozygote','OR Nonaggressive Homozygote',
                'OR Aggressive Heterozygote','OR Aggressive Homozygote'])
  for row in associations:
    snpannoid = row[4]
    snpinfo = snpmap.get(snpannoid)
    genes = genemap.get(snpannoid,None)
    name = anamap[row[3]]
    if genes is not None and len(genes) == 1:
      genes = genes[0]
    elif genes is not None and len(genes) > 1:
      genes = '|'.join(genes)
    out.writerow(list(chain(snpinfo,[genes,name,row[1],row[2]],row[5:9])))

if __name__ == '__main__':
  main()
