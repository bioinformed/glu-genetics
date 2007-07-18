# -*- coding: utf-8 -*-
'''
File:          load_db.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import csv
import glob
import sqlite3

from   itertools         import islice,chain

from   glu.lib.utils     import autofile
from   glu.lib.sections  import read_sections


HAPMAP='/usr/local/share/hapmap/build22/rs_strand/non-redundant/geno*'
MANIFEST='/usr/local/share/manifests/HumanHap550_11218540_C.csv'


def clean_alias(a):
  a=a.strip()
  if len(a)>2 and a[0] in ('"',"'") and a[0] == a[-1]:
    a = a[1:-1]
  if a == '-':
    a = ''
  return a


def get_aliases(genes, filename):
  geneids = set()
  for gene in genes:
    geneid=gene[0]
    symbol=gene[1]
    geneids.add(geneid)
    yield symbol,geneid,symbol

  aliasfile = csv.reader(autofile(filename), dialect='excel-tab')
  aliasfile.next()

  for record in aliasfile:
    symbol = record[2].upper()
    geneid  = int(record[1])

    if record[0] != '9606' or record[-1] == '-' or record[2] != symbol or geneid not in geneids:
      continue

    yield symbol,geneid,symbol

  aliasfile = csv.reader(autofile(filename), dialect='excel-tab')
  aliasfile.next()

  for record in aliasfile:
    symbol = record[2].upper()
    geneid  = int(record[1])

    if record[0] != '9606' or record[-1] == '-' or record[2] != symbol or geneid not in geneids:
      continue

    aliases = record[4].upper()

    aliases = [ clean_alias(a) for a in aliases.split('|') ]

    for alias in aliases:
      if alias and alias != symbol and alias == alias.upper():
        yield alias,geneid,symbol


def filter_aliases(aliases):
  aliases = list(aliases)

  seen = {}
  for a,gid,s in aliases:
    seen.setdefault(a,gid)

  used = set()
  for a,gid,s in aliases:
    if a not in used and seen[a] == gid:
      used.add(a)
      yield a,gid


def get_genes(filename):
  genes = csv.reader(autofile(filename),dialect='excel-tab')
  header = genes.next()

  for taxid,chromosome,chrStart,chrEnd,orientation,contig,cnt_start,     \
      cnt_stop,cnt_orient,featureName,featureId,featureType,groupLabel,  \
      transcript,evidence in genes:

    if groupLabel != 'reference' or featureType != 'GENE':
      continue

    if '|' in chromosome:
      chromosome = chromosome.split('|')[0]
      chrStart = chrEnd = None

    assert featureId.startswith('GeneID:')
    geneid     = int(featureId[7:])
    if chrStart:
      chrStart = int(chrStart)
    if chrEnd:
      chrEnd   = int(chrEnd)
    cnt_start  = int(cnt_start)
    cnt_stop   = int(cnt_stop)

    yield geneid,featureName,chromosome,chrStart,chrEnd,orientation,contig,cnt_start,\
          cnt_stop,cnt_orient,featureType,groupLabel,transcript,evidence


def load_genes(con,genes):
  cur = con.cursor()

  try:
    cur.execute('DROP TABLE GENE;')
  except:
    pass

  sql = '''
  create table GENE (geneid        INTEGER,
                     featureName   TEXT,
                     chromosome    TEXT,
                     chrStart      INTEGER,
                     chrEnd        INTEGER,
                     orientation   TEXT,
                     contig        TEXT,
                     cnt_start     INTEGER,
                     cnt_end       INTEGER,
                     cnt_orient    TEXT,
                     featureType   TEXT,
                     submitGroup   TEXT,
                     transcript    TEXT,
                     evidence      TEXT);'''


  cur.execute(sql)
  sql = 'INSERT INTO GENE VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
  cur.executemany(sql, genes)

  sql = 'CREATE INDEX idx_gene_id ON GENE (geneid,featureType)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_name ON GENE (featureName,featureType)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_loc ON GENE (chromosome,chrStart,chrEnd)';
  cur.execute(sql)

  sql = 'SELECT COUNT(*) FROM GENE;'
  cur.execute(sql)
  print cur.fetchall()

  con.commit()



def load_aliases(con, aliases):
  cur = con.cursor()

  try:
    cur.execute('DROP TABLE ALIAS;')
  except:
    pass

  cur.execute('CREATE TABLE ALIAS (alias TEXT PRIMARY KEY, geneid INTEGER);')

  sql = 'INSERT INTO ALIAS VALUES (?,?);'
  cur.executemany(sql, aliases)

  sql = 'SELECT COUNT(*) FROM ALIAS;'
  cur.execute(sql)
  print cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_aliases ON ALIAS (alias);'
  cur.execute(sql)

  con.commit()


def fix_aliases(con):
  sql = '''
  SELECT   s.featureName, s.geneID
  FROM     gene s
  WHERE    NOT EXISTS (SELECT * FROM alias a WHERE a.Alias == s.featureName)
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference';
  '''

  cur = con.cursor()
  cur.execute(sql)
  unmapped = cur.fetchall()

  counts = {}
  for symbol,geneid in unmapped:
    counts[symbol] = counts.get(symbol,0) + 1

  unmapped = [ (s,gid) for s,gid in unmapped if counts[s] == 1 ]

  print 'UNMAPPED LEN: %d' % len(unmapped)

  cur = con.cursor()
  sql = "INSERT INTO ALIAS VALUES (?,?);"
  cur.executemany(sql, unmapped)

  sql = 'SELECT COUNT(*) FROM ALIAS;'
  cur.execute(sql)
  print cur.fetchall()[0][0]

  con.commit()


def get_snps(filename):
  filenames = glob.glob(filename)
  for filename in filenames:
    print filename
    gfile = autofile(filename)
    header = gfile.next()
    for line in gfile:
      fields     = line.split(' ')
      rs         = fields[0]
      chromosome = fields[2]
      strand     = fields[4]

      if chromosome.startswith('chr'):
        chromosome = chromosome[3:]
      position = int(fields[3])
      yield rs,chromosome,position,strand


def squash_dups(snps):
  seen = set()
  for s in snps:
    if s[0] not in seen:
      seen.add(s[0])
      yield s


def load_snps(con, snps):
  cur = con.cursor()

  try:
    cur.execute('DROP TABLE SNP;')
  except:
    pass

  cur.execute('CREATE TABLE SNP (lname TEXT PRIMARY KEY, chromosome TEXT, location INTEGER, strand TEXT);')

  sql = 'INSERT INTO SNP VALUES (?,?,?,?);'
  cur.executemany(sql, snps)

  sql = 'SELECT COUNT(*) FROM SNP;'
  cur.execute(sql)
  print cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_snp ON SNP (chromosome,location);'
  cur.execute(sql)

  con.commit()


def load_illumina_manifest(filename):
  ifile = csv.reader(autofile(filename),dialect='excel')
  sections = read_sections(ifile)

  heading,contents = sections.next()

  # OPA manifest
  if heading == 'data':
    heading,contents = sections.next()
    assert heading == 'Heading'

    headings = dict(islice(contents,10))
    assert headings['Assay Format'] == 'Golden Gate'

  elif heading == 'Heading':
    contents = dict(contents)
    assert contents['Assay Format'] == 'Infinium II'

    heading,contents = sections.next()
    assert heading == 'Assay'

  return contents


def find_index(header,headings):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass
  raise ValueError,'Cannot find heading index'


def extract_illumina_snps(manifest):
  manifest    = iter(manifest)
  header      = manifest.next()

  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  chr_idx     = find_index(header,['Chr'])
  loc_idx     = find_index(header,['MapInfo'])

  strandmap   = {'F':'+','R':'-','U':'+'}

  for assay in manifest:
    rs         = assay[name_idx]
    chromosome = assay[chr_idx]
    position   = int(assay[loc_idx])
    strand     = strandmap[assay[assayid_idx].split('_')[2]]
    yield rs,chromosome,position,strand


def main():
  con = sqlite3.connect('genome36.db')

  if 1:
    genes = list(get_genes('data/seq_gene.md.b36.2.gz'))
    load_genes(con,genes)

  if 1:
    aliases = get_aliases(genes,'data/gene_info.gz')
    aliases = filter_aliases(aliases)
    load_aliases(con,aliases)

    fix_aliases(con)

  if 1:
    hapmap_snps   = get_snps(HAPMAP)
    #illumina_snps = extract_illumina_snps(load_illumina_manifest(MANIFEST))
    illumina_snps = []

    snps = chain(hapmap_snps,illumina_snps)
    snps = squash_dups(snps)
    load_snps(con,snps)


if __name__ == '__main__':
  main()
