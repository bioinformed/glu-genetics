'''
File:          phase2marga.py

Authors:       Zhaoming Wang  (wangzha@mail.nih.gov)

Created:       2007-09-21

Abstract:      Sample a pair of haploytypes for each DNA according to the frequency of each haplotype reconstructed by PHASE
               Usage: python phase2marga.py arg1 arg2 arg3 arg4 arg5
               arg1:  phase output file
               arg2:  phenotype map file
               arg3:  physical or genetic distances map file
               arg4:  allele map file
               arg5:  output file in Margarita input format

Requires:      Python 2.5, glu

Revision:      $Id: $
'''

import sys
import re
import csv
from itertools         import chain,izip
from random            import Random
from glu.lib.fileutils import load_map,load_list

spaces = re.compile('[\t ,]+')

def pickone(id,choices,rand):
  alist = []
  for i,choice in enumerate(choices):
    alist.extend([i] * int(100 * float(choice[2])))
  rand.shuffle(alist)
  return id,choices[alist[0]][0:2]


def haplosampler(filename):
  rand = Random()
  infile = file(filename)
  state = 'start'
  choices = []
  id = None
  for line in infile:
    if line.startswith('IND') and state == 'start':
      id = line.strip().split()[1]
      state = 'id'
    elif line.startswith('IND') and state == 'haplo':
      one = pickone(id,choices,rand)
      choices = []
      id = line.strip().split()[1]
      state = 'id'
      yield one
    elif not line.startswith('IND') and state == 'haplo':
      fields = spaces.split(line.strip())
      choices.append(fields)
    elif not line.startswith('IND') and state == 'id':
      fields = spaces.split(line.strip())
      choices.append(fields)
      state = 'haplo'
    else:
      print >> sys.stderr, 'Wrong format of the input file'


def recode(pairs,alleles):
  val1 = []
  val2 = []
  hap1,hap2 = pairs
  for h1,h2,allele in izip(hap1,hap2,alleles):
    a1,a2 = allele.split(',')
    if (h1,h2) == (a1,a2):
      val1.append('0')
      val2.append('1')
    elif (h1,h2) == (a1,a1):
      val1.append('0')
      val2.append('0')
    elif (h1,h2) == (a2,a2):
      val1.append('1')
      val2.append('1')
    elif (h1,h2) == (a2,a1):
      val1.append('1')
      val2.append('0')
    else:
      print >> sys.stderr, 'Encounter %s and %s, expecting %s and %s' % (h1,h2,a1,a2)
  return [''.join(val1),''.join(val2)]


def load_alleles(filename):
  rows = csv.reader(file(filename),dialect='excel-tab')
  return [row[1] for row in rows]


def main():
  haplos = haplosampler(sys.argv[1])  # haploytpe file
  #todo: add a list of individual id to filter  them out
  phenos = load_map(sys.argv[2],key_index=0,value_index=8) # subject.def column 0 and 8, need to be generalized!
  distances = load_list(sys.argv[3],index=1,skip=1) # try to use physical distance
  alleles = load_alleles(sys.argv[4]) # note the order has to be the same as SNPs in haplotype

  cases=[]
  controls=[]
  for id, pairs in haplos:
    if phenos[id] == 'CASE':
      cases.extend(recode(pairs,alleles))
    elif phenos[id] == 'CONTROL':
      controls.extend(recode(pairs,alleles))

  nummarkers = len(distances)
  numcases = len(cases)
  numcontrols = len(controls)

  outfile = file(sys.argv[5],'w')
  outfile.write('%d %d %d\n' % (numcases,numcontrols,nummarkers))
  outfile.write('%s\n' % '\n'.join(distances))
  outfile.write('%s\n' % '\n'.join(cases))
  outfile.write('%s\n' % '\n'.join(controls))


if __name__ == '__main__':
  main()
