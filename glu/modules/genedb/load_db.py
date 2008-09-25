# -*- coding: utf-8 -*-

__abstract__  = 'Build a genedb database'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import glob
import sqlite3

from   itertools         import chain,islice

from   glu.lib.fileutils import autofile, load_table, tryint

from   glu.modules.convert.from_lbd import load_illumina_manifest


DBSNP     = ['data/snp129.txt.gz', 'data/snp128.txt.gz', 'data/snp126.txt.gz']
ARRAYS    = glob.glob('data/snpArray*')
HAPMAP    = glob.glob('/usr/local/share/hapmap/build23/rs_strand/non-redundant/geno*')
MANIFESTS = ['/usr/local/share/manifests/Human1M-Duov3_B.csv',
             '/usr/local/share/manifests/Human1Mv1_C.csv',
             '/usr/local/share/manifests/Human610-Quadv1_B.csv']


def clean_alias(a):
  a=a.strip()
  if len(a)>2 and a[0] in ('"',"'") and a[0] == a[-1]:
    a = a[1:-1]
  if a == '-':
    a = ''
  return a


def get_aliases(con,filename):
  sql = 'SELECT geneid,featureName FROM GENE;'
  cur = con.cursor()
  cur.execute(sql)
  geneids = dict( (tryint(geneid),symbol) for geneid,symbol in cur.fetchall() )

  aliasfile = load_table(filename,want_header=True)
  aliasfile.next()

  updates = []
  for record in aliasfile:
    symbol = record[2]
    geneid  = int(record[1])

    if record[0] != '9606' or record[-1] == '-' or geneid not in geneids:
      continue

    if symbol != geneids[geneid]:
      geneids[geneid] = symbol
      updates.append( (symbol,geneid) )

  print 'Updating map with %d new official feature names...' % len(updates),
  sql = 'UPDATE GENE SET featureName=? WHERE geneid=?;'
  cur.executemany(sql,updates)
  con.commit()
  print 'Done.'

  for geneid,symbol in geneids.iteritems():
    yield symbol,geneid,symbol
    if symbol != symbol.upper():
      yield symbol.upper(),geneid,symbol.upper()
    yield str(geneid),geneid,symbol

  aliasfile = load_table(filename,want_header=True)
  aliasfile.next()

  for record in aliasfile:
    symbol = record[2]
    geneid  = int(record[1])

    if record[0] != '9606' or record[-1] == '-' or geneid not in geneids:
      continue

    aliases = record[4]

    aliases = [ clean_alias(a) for a in aliases.split('|') ]

    for alias in aliases:
      if alias and alias != symbol:
        yield alias,geneid,symbol
        if alias != alias.upper():
          yield alias.upper(),geneid,symbol


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
  genes  = load_table(filename,want_header=True)
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
  print 'GENES:',cur.fetchall()[0][0]

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
  print 'ALIASES:',cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_aliases ON ALIAS (alias);'
  cur.execute(sql)

  con.commit()


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

  cur.execute('''
    CREATE TABLE SNP (lname      TEXT PRIMARY KEY,
                      chromosome TEXT,
                      location   INTEGER,
                      strand     TEXT,
                      alleles    TEXT);''')

  sql = 'INSERT INTO SNP VALUES (?,?,?,?,?);'

  while 1:
    rows = list(islice(snps,50000))
    if not rows:
      break
    cur.executemany(sql, rows)
    con.commit()

  sql = 'SELECT COUNT(*) FROM SNP;'
  cur.execute(sql)
  print 'SNPS:',cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_snp ON SNP (chromosome,location);'
  cur.execute(sql)

  con.commit()


def get_hapmap_snps(filename):
  print 'LOADING HAPMAP:',filename
  gfile  = autofile(filename)
  header = gfile.next()
  for line in gfile:
    fields     = line.split(' ')
    rs         = fields[0]
    chromosome = fields[2]
    strand     = fields[4]
    alleles    = fields[1]
    position   = int(fields[3])

    if chromosome.startswith('chr'):
      chromosome = chromosome[3:]

    yield rs,chromosome,position,strand,alleles


def find_index(header,headings):
  for h in headings:
    try:
      return header.index(h)
    except ValueError:
      pass
  raise ValueError('Cannot find heading index')


def extract_illumina_snps(manifest):
  manifest    = iter(manifest)
  header      = manifest.next()

  assayid_idx = find_index(header,['IlmnID','Ilmn ID'])
  name_idx    = find_index(header,['Name'])
  chr_idx     = find_index(header,['Chr'])
  loc_idx     = find_index(header,['MapInfo'])
  snp_idx     = find_index(header,['SNP'])

  strandmap   = {'F':'+','R':'-','U':'+'}

  for assay in manifest:
    rs         = assay[name_idx]
    chromosome = assay[chr_idx]
    position   = int(assay[loc_idx])
    strand     = strandmap[assay[assayid_idx].split('_')[-2]]
    alleles    = assay[snp_idx][1:-1]

    if chromosome=='Mt':
      chromosome='M'

    yield rs,chromosome,position,strand,alleles


def get_goldenpath_dbsnp(snpfile):
  print 'LOADING DBSNP:',snpfile
  for row in load_table(snpfile):
    rs         = row[4]
    chromosome = row[1]
    position   = tryint(row[2])
    strand     = row[6]
    alleles    = row[9]

    # Adjust for 0-based indexing
    if position is not None:
      position += 1

    if chromosome.startswith('chr'):
      chromosome = chromosome[3:]

    yield rs,chromosome,position,strand,alleles


def get_goldenpath_arrays(arrayfile):
  print 'LOADING ARRAY:',arrayfile
  for row in load_table(arrayfile):
    name       = row[4]
    chromosome = row[1]
    position   = tryint(row[2])
    strand     = row[6]
    alleles    = row[7]

    # Adjust for 0-based indexing
    if position is not None:
      position += 1

    if chromosome.startswith('chr'):
      chromosome = chromosome[3:]

    # If rsId is available, output only that
    if len(row)>=9:
      rs = row[8]
      yield rs,chromosome,position,strand,alleles
    else:
      yield name,chromosome,position,strand,alleles


def get_cytobands(cytofile):
  print 'LOADING CYTOBANDS:',cytofile
  seen = set()
  for row in load_table(cytofile):
    band       = row[3]
    chromosome = row[0]
    start      = tryint(row[1])
    stop       = tryint(row[2])
    color      = row[4]

    # Adjust for 0-based indexing
    if start is not None:
      start += 1

    if chromosome.startswith('chr'):
      chromosome = chromosome[3:]

    name = chromosome + band

    if name in seen:
      print name,row

    seen.add(name)

    yield name,chromosome,start,stop,color


def load_cytobands(con,bands):
  cur = con.cursor()

  try:
    cur.execute('DROP TABLE CYTOBAND;')
  except:
    pass

  cur.execute('''
    CREATE TABLE CYTOBAND (band       TEXT PRIMARY KEY,
                           chromosome TEXT,
                           start      INTEGER,
                           stop       INTEGER,
                           color      TEXT);''')

  sql = 'INSERT INTO CYTOBAND VALUES (?,?,?,?,?);'
  cur.executemany(sql, bands)

  sql = 'SELECT COUNT(*) FROM CYTOBAND;'
  cur.execute(sql)
  print 'BANDS:',cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_cytoband ON CYTOBAND (chromosome,start,stop);'
  cur.execute(sql)

  con.commit()


def get_mirbase(mfile):
  print 'LOADING MIRBASE'
  for row in load_table(mfile):
    if not row or row[0].startswith('#'):
      continue

    names       = [ tuple(n.strip().replace('"','').split('=')) for n in row[8].strip().split(';') ]
    names       = dict(names[:-1])
    geneid      = names['ACC']
    featureName = names['ID']
    chromosome  = row[0]
    chrStart    = int(row[3])-1
    chrEnd      = int(row[4])-1
    orientation = row[6]
    contig      = None
    cnt_start   = None
    cnt_stop    = None
    cnt_orient  = None
    featureType = row[2]
    groupLabel  = 'reference'
    transcript  = None
    evidence    = None

    yield geneid,featureName,chromosome,chrStart,chrEnd,orientation,contig,cnt_start,\
          cnt_stop,cnt_orient,featureType,groupLabel,transcript,evidence


def main():
  con = sqlite3.connect('db/b2/genedb_ncbi36.3_huge.db')

  con.execute('PRAGMA synchronous=OFF;')
  con.execute('PRAGMA journal_mode=OFF;')
  con.execute('PRAGMA cache_size=20000;')

  if 1:
    genes = [ get_genes('data/seq_gene.md.b36.3.gz'),
              get_mirbase('data/hsa-v12.gff') ]
    genes = list(chain(*genes))
    load_genes(con,genes)

  if 1:
    aliases = get_aliases(con,'data/gene_info.gz')
    aliases = list(filter_aliases(aliases))
    load_aliases(con,aliases)

  if 1:
    bands = get_cytobands('data/cytoBand.txt.gz')
    load_cytobands(con,bands)

  if 1:
    streams  = []
    streams += [ get_goldenpath_dbsnp(s)  for s in DBSNP  ]
    streams += [ get_hapmap_snps(h)       for h in HAPMAP ]
    streams += [ get_goldenpath_arrays(a) for a in ARRAYS ]
    streams += [ extract_illumina_snps(load_illumina_manifest(m))
                   for m in MANIFESTS ]

    snps = chain(*streams)
    snps = squash_dups(snps)
    load_snps(con,snps)


if __name__ == '__main__':
  main()
