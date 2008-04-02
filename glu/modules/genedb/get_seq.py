# -*- coding: utf-8 -*-
'''
File:          get_seq.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       Wed May 17 08:51:57 EDT 2006

Abstract:

Compatibility: Python 2.5 and above

Requires:      No external dependencies, yet...

Revision:      $Id$
'''

__program__   = 'get_seq'
__authors__   = ['Kevin Jacobs (jacobs@bioinformed.com)']
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sqlite3
import time
import csv
import os

import seq_comp

from tagzilla import *

HAPMAP_PATH= '/home/jacobske/projects/CGEMS/hapmap'
GENOTYPES  = os.path.join(HAPMAP_PATH,'build20/non-redundant')
PEDIGREES  = os.path.join(HAPMAP_PATH,'pedigrees', 'peds')
POPS       = ['CEU','YRI','JPT','CHB','JPT+CHB']

def escape(s):
  return "'%s'" % s.replace("'","''")


def get_sequences(con,snps):
  cur = con.cursor()

  try:
    cur.execute('CREATE TABLE snpsequence (lname TEXT PRIMARY KEY, sequence TEXT);')
  except:
    pass

  sql = '''
  SELECT  lname,sequence
  FROM    snpsequence
  WHERE   lname IN (%s);
  '''
  sqlsnps = ','.join(escape(snp) for snp in snps if snp)
  cur.execute(sql % sqlsnps)
  results  = cur.fetchall()

  missing  = snps - set(lname for lname,seq in results)
  results2 = [ s for s in get_sequences_from_genewindow(missing) if len(s) == 2 ]

  missing -= set(lname for lname,seq in results2)
  results3 = ((lname,'') for lname in missing)

  sql = 'INSERT INTO snpsequence VALUES (?,?);'
  cur.executemany(sql, chain(results2,results3))

  con.commit()

  return ( (lname,seq) for lname,seq in chain(results,results2) if seq)


def get_sequences_from_genewindow(snps):
  snps = list(snps)
  if not snps:
    return []

  command = 'java -Xmx1000M -classpath sequence/gw_tools.jar:sequence/classes12.jar:sequence nci/cgf/annotator/tools/export/BatchExport'
  snpname = 'tmpFoo%d' % time.time()
  sequencename = 'tmpBar%d' % time.time()
  try:
    snpfile = file(snpname,'w')
    for snp in snps:
      snpfile.write('%s\n' % snp)
    snpfile.close()
    args = '%s sequence/sublist.txt sequence/oblig.txt %s' % (snpname,sequencename)
    os.popen('%s %s 2>/dev/null' % (command,args)).read()
    return csv.reader(open(sequencename),dialect='excel-tab')
  finally:
    os.unlink(snpname)
    try:
      os.unlink(sequencename)
    except:
      pass



def extend(s,n):
  m = len(s)
  if m < n:
    s = list(s)
    s.extend(['']*(n-m))
  return s


def option_parser():
  import optparse

  usage = 'usage: %prog [options] genofile...'
  parser = optparse.OptionParser(usage=usage, add_help_option=False)

  parser.add_option('-h', '--help', dest='help', action='store_true',
                          help='show this help message and exit')
  parser.add_option('--license', dest='license', action='store_true',
                          help="show program's copyright and license terms and exit")
  return parser


def get_seq(options,args):
  con = sqlite3.connect('genome35.db')

  snps = set()
  for arg in args:
    snps.update(row[0] for row in csv.reader(autofile(arg),dialect='excel-tab') if row)

  degenerate = set('wymkwsbdhvnx')
  dcount = 0
  n = 16
  for rs,seq in get_sequences(con, snps):
    i = seq.index('[')
    left = seq[i-n:i].lower()
    i = seq.index(']')
    right = seq[i+1:i+n+1].lower()
    bases = set(left+right)
    if bases & degenerate:
      dcount += 1
    print '>%s\n%s' % (rs,seq)

  print 'Degenerate:',dcount

def main():
  launcher(get_seq, option_parser, **globals())


if __name__ == '__main__':
  main()
