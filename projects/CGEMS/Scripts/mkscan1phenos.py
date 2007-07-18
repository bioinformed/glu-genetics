import csv
import sys
import random

from   numpy          import array,matrix,exp,vsplit,zeros,where,hstack,vstack,log,exp,maximum,nan
from   itertools      import izip,imap,repeat,islice

from   biozilla.utils import pick,tally,autofile


def normalize(row,n):
  m = len(row)

  if n == m:
    return row
  elif n < m:
    return row[:n]
  else:
    return row + ['']*(n-m)


def load_phenos(filename):
  phenos = csv.reader(autofile(filename),dialect='excel-tab')
  header = phenos.next()
  return header,imap(normalize, phenos, repeat(len(header)))


def make_index(header):
  return dict( (f,i) for i,f in enumerate(header) )


def filter_pheno_columns(header, phenos, keep):
  indices = [ header.index(k) for k in keep ]

  def _filter_phenos():
    for pheno in phenos:
      yield [ pheno[i] for i in indices ]

  return keep,_filter_phenos()


def expand_phenos(header, phenos):
  new_header = header + ['YR','STAT','PREVALENT']
  idx = make_index(new_header)

  ccmap = {'1':'0','2':'1','3':'2','4':'2'}
  drops = set(r[0] for r in csv.reader(file('drop_samples'),dialect='excel-tab'))

  def _expand_phenos():
    for pheno in phenos:
      cc    = ccmap[pheno[idx['CC_OVERALL']]]
      yrs   = pheno[idx['YR1']],pheno[idx['YR2']],pheno[idx['YR3']]
      stats = pheno[idx['STAT1']],pheno[idx['STAT2']],pheno[idx['STAT3']]

      pid = pheno[idx['PID']]

      if pid in drops:
        continue

      selected = False
      for yr,stat in zip(yrs,stats):
        prev = '1' if yr=='0' else '0'
        if stat == '1':
          yield pheno + [yr,'0',prev]
        elif stat == '2':
          assert cc in ('1','2')
          selected = True
          yield pheno + [yr,cc,prev]
          break

      if cc in ('1','2') and not selected:
        print >> sys.stderr, pheno[0],cc

  return new_header,_expand_phenos()


def expand_phenos_simple(header, phenos):
  new_header = header + ['YR','STAT','PREVALENT']
  idx = make_index(new_header)

  ccmap = {'1':'0','2':'1','3':'2','4':'2'}

  drops = set(r[0] for r in csv.reader(file('drop_samples'),dialect='excel-tab'))

  def _expand_phenos():
    for pheno in phenos:
      cc    = ccmap[pheno[idx['CC_OVERALL']]]
      yrs   = pheno[idx['YR1']],pheno[idx['YR2']],pheno[idx['YR3']]
      stats = pheno[idx['STAT1']],pheno[idx['STAT2']],pheno[idx['STAT3']]

      prev = 0
      for yr,stat in zip(yrs,stats):
        if yr == 0:
          prev = '1'
          break

      pid = pheno[idx['PID']]
      if pid not in drops:
        yield pheno + ['0',cc,prev]

  return new_header,_expand_phenos()


def filter_pheno_records(header,phenos):
  def _filter_pheno_records():
    idx = make_index(header)

    for pheno in phenos:
      consent = pheno[idx['HAS_CONSENT']]
      samples = pheno[idx['GENOTYPED_SAMPLE_IDS']]
      center  = pheno[idx['CENTER']]

      if consent == '1' and samples and center != '3':
        yield pheno

  return header,_filter_pheno_records()


def expand_nominal(header,phenos,field,valuemap):
  values = sorted(set(valuemap.values()))
  header = header + values
  pheno_idx = make_index(header)
  value_idx = make_index(values)

  def _expand_nominal():
    for pheno in phenos:
      vals = ['0']*len(values)
      value_field = valuemap[pheno[pheno_idx[field]]]
      vals[ value_idx[value_field] ] = '1'
      yield pheno + vals

  return header,_expand_nominal()


def compute_maf(counts):
  hom1 = hom2 = het = 0

  for g,n in counts.iteritems():
    if g[0] != g[1]:
      het = n
    elif hom1:
      hom2 = n
    else:
      hom1 = n

  hom1,hom2 = min(hom1,hom2),max(hom1,hom2)

  maf = float(2*hom1 + het)/(hom1+het+hom2)/2
  return maf


def main():
  header,phenos = load_phenos(sys.argv[1])
  header = map(str.upper,header)
  header,phenos = filter_pheno_records(header,phenos)

  if 1:
    header,phenos = expand_phenos(header,phenos)
  else:
    header,phenos = expand_phenos_simple(header,phenos)

  agemap = dict( (str(i),'AGELEVEL%s' % i) for i in range(4) )
  header,phenos = expand_nominal(header,phenos,'AGELEVEL',agemap)

  centermap = dict( (str(i),'CENTER%02d' % i) for i in range(12) )
  header,phenos = expand_nominal(header,phenos,'CENTER',centermap)

  yrmap = dict( (str(i),'YRIND%d' % i) for i in range(9) )
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  yrmap = dict( (str(i),'INCIDENT') for i in range(9) )
  yrmap['0'] = 'PREV'
  header,phenos = expand_nominal(header,phenos,'YR',yrmap)

  keep = ['PID','STAT',
          'CENTER01','CENTER02','CENTER04','CENTER05','CENTER06',
          'CENTER08','CENTER09','CENTER10', #'PREV',
          'AGELEVEL1','AGELEVEL2','AGELEVEL3']

  #keep = ['PID','STAT']
  keep = ['PID','STAT','CENTER','AGELEVEL']

  header,phenos = filter_pheno_columns(header,phenos,keep)
  phenos = list(phenos)

  pids = set(csv.reader(autofile(sys.argv[2]),dialect='excel-tab').next()[1:])

  f = csv.writer(file('phenos.csv','w'),dialect='excel-tab')
  f.writerow(header)
  f.writerows(pheno for pheno in phenos if pheno[0] in pids)


if __name__ == '__main__':
  main()
