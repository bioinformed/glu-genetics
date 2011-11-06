# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Build a genedb database'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id: load_db.py 960 2009-09-22 13:18:52Z bioinformed $'


import glob
import time
import sqlite3

from   itertools            import chain, islice, imap
from   collections          import defaultdict, namedtuple

from   glu.lib.utils        import gcdisabled
from   glu.lib.fileutils    import autofile, table_reader, map_reader, list_reader, tryint
from   glu.lib.glu_launcher import format_elapsed_time
from   glu.lib.illumina     import IlluminaManifest


DBSNP     = ['data/snp130.txt.gz', 'data/snp129.txt.gz', 'data/snp128.txt.gz', 'data/snp126.txt.gz']
ARRAYS    = glob.glob('data/snpArray*')
HAPMAP    = glob.glob('/usr/local/share/hapmap/build23/rs_strand/non-redundant/geno*')
MANIFESTS = ['/usr/local/share/manifests/Human1M-Duov3_B.csv',
             '/usr/local/share/manifests/Human1Mv1_C.csv',
             '/usr/local/share/manifests/Human610-Quadv1_B.csv']


COMP = {'A':'T','T':'A','C':'G','G':'C','N':'N','-':'-'}


def rev_comp(s):
  comp =  [ COMP[c] for c in s ]
  comp.reverse()
  return ''.join(comp)


def sqlite_magic(con):
  con.execute('PRAGMA synchronous=OFF;')
  con.execute('PRAGMA journal_mode=OFF;')
  con.execute('PRAGMA count_changes=OFF;')
  con.execute('PRAGMA cache_size=200000;')
  con.execute('PRAGMA default_cache_size=200000;')


def batch_write(con,sql,values,chunk=100000):
  cur = con.cursor()
  sqlite_magic(con)
  total_rows = 0
  very_start = time.time()

  while 1:
    start = time.time()
    rows = list(islice(values,chunk))

    if not rows:
      break

    cur.executemany(sql, rows)
    t1 = time.time()
    con.commit()
    t2 = time.time()

    total_rows   += len(rows)
    total_time    = format_elapsed_time(t2-very_start)
    total_rate    = total_rows/(t2-very_start)
    execute_time  = format_elapsed_time(t1-start)
    commit_time   = format_elapsed_time(t2-t1)
    batch_rate    = len(rows)/(t2-start)

    print '[total] time: %s rows/s %.1f; [batch] exectute: %s, commit: %s, rows/s: %.1f' % \
              (total_time,total_rate,execute_time,commit_time,batch_rate)


def get_goldenpath_aliases(con,filename,seen):
  for alias,name in table_reader(filename, columns=[1,0]):
    seen.add(alias)
    yield alias,name


def clean_alias(a):
  a=a.strip()
  if len(a)>2 and a[0] in ('"',"'") and a[0] == a[-1]:
    a = a[1:-1]
  if a == '-':
    a = ''
  return a


def get_entrez_aliases(con,filename,seen):
  sql = 'SELECT geneid,name FROM GENE;'
  cur = con.cursor()
  cur.execute(sql)

  geneids = defaultdict(list)
  for geneid,name in cur.fetchall():
    geneids[tryint(geneid)].append(name)

  aliasfile = table_reader(filename,want_header=True)
  aliasfile.next()

  for record in aliasfile:
    symbol = record[2]
    geneid = int(record[1])

    if record[0] != '9606' or record[-1] == '-' or geneid not in geneids:
      continue

    names = geneids[geneid]

    aliases = [symbol]+[ clean_alias(a) for a in record[4].split('|') ]

    for alias in aliases:
      if not alias:
        continue
      for a in [alias,alias.upper()]:
        if a not in seen:
          seen.add(a)
          for name in names:
            yield a,name


def get_goldenpath_genes(dirname):
  known_header    = ['name','chrom','strand','txStart','txEnd','cdsStart','cdsEnd',
                     'exonCount','exonStarts','exonEnds','proteinID','alignID']

  xref_header     = ['kgID','mRNA','spID','spDisplayID','geneSymbol','refseq','protAcc',
                     'description']

  genes         = table_reader(dirname+'/knownGene.txt.gz',header=known_header)
  geneid_link   = map_reader(dirname+'/knownToLocusLink.txt.gz')
  symbol_link   = map_reader(dirname+'/kgXref.txt.gz', columns=[0,4])
  mrna_link     = map_reader(dirname+'/kgXref.txt.gz', columns=[0,1])
  sp_link       = map_reader(dirname+'/kgXref.txt.gz', columns=[0,3])
  ncbi_link     = map_reader(dirname+'/kgXref.txt.gz', columns=[0,6])
  canonical_set = set(list_reader(dirname+'/knownCanonical.txt.gz', columns=[4]))

  for name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd, \
      exonCount,exonStarts,exonEnds,proteinID,alignID in genes:

    txStart     = int(txStart)   if txStart   else None
    txEnd       = int(txEnd)     if txEnd     else None
    cdsStart    = int(cdsStart)  if cdsStart  else None
    cdsEnd      = int(cdsEnd)    if cdsEnd    else None
    exonCount   = int(exonCount) if exonCount else None

    geneID      = geneid_link.get(name)
    symbol      = symbol_link.get(name)
    mrnaAcc     = mrna_link.get(name)
    protAcc     = proteinID or sp_link.get(name) or ncbi_link.get(name)
    random      = chrom.endswith('_random')
    canonical   = (name in canonical_set) | (not random)<<1 | (not random and '_' not in chrom)<<2

    if not symbol:
      print 'No gene symbol for:',name
    if not mrnaAcc:
      print 'No mrnaAcc for:',name

    yield name,symbol,geneID,mrnaAcc,proteinID,canonical,chrom,strand,txStart,txEnd, \
          cdsStart,cdsEnd,exonCount,exonStarts,exonEnds


def load_genes(con,genes):
  cur = con.cursor()

  try:
    cur.execute('DROP INDEX idx_gene_name;')
  except:
    pass

  try:
    cur.execute('DROP INDEX idx_gene_id;')
  except:
    pass

  try:
    cur.execute('DROP INDEX idx_gene_symbol;')
  except:
    pass

  try:
    cur.execute('DROP INDEX idx_gene_loc;')
  except:
    pass

  try:
    cur.execute('DROP TABLE GENE;')
  except:
    pass

  sql = '''
  CREATE TABLE GENE (id            INTEGER PRIMARY KEY,
                     name          TEXT,
                     symbol        TEXT,
                     geneID        INTEGER,
                     mRNA          TEXT,
                     protein       TEXT,
                     canonical     INTEGER,
                     chrom         TEXT,
                     strand        TEXT,
                     txStart       INTEGER,
                     txEnd         INTEGER,
                     cdsStart      INTEGER,
                     cdsEnd        INTEGER,
                     exonCount     INTEGER,
                     exonStarts    TEXT,
                     exonEnds      TEXT);'''

  cur.execute(sql)
  sql = 'INSERT INTO GENE VALUES (NULL,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'
  batch_write(con,sql,genes)

  sql = 'CREATE UNIQUE INDEX idx_gene_id ON GENE (id)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_name ON GENE (name)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_geneid ON GENE (geneid)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_symbol ON GENE (symbol)';
  cur.execute(sql)

  sql = 'CREATE INDEX idx_gene_loc ON GENE (chrom,txStart,txEnd)';
  cur.execute(sql)

  sql = 'SELECT COUNT(*) FROM GENE;'
  cur.execute(sql)
  print 'GENES:',cur.fetchall()[0][0]

  con.commit()

  try:
    cur.execute('DROP TABLE gene_index;')
  except:
    pass

  print 'Creating GENE rtree index:'
  sql = '''CREATE VIRTUAL TABLE gene_index USING rtree_i32(id,txStart,txEnd);'''
  cur.execute(sql)

  sql = '''
  INSERT INTO gene_index(id,txStart,txEnd)
  SELECT id,txStart,txEnd
  FROM   gene;
  '''
  cur.execute(sql)

  con.commit()


def load_gene_aliases(con, aliases):
  cur = con.cursor()

  try:
    cur.execute('DROP INDEX idx_aliases;')
  except:
    pass

  try:
    cur.execute('DROP TABLE ALIAS;')
  except:
    pass

  cur.execute('CREATE TABLE ALIAS (alias TEXT, name TEXT);')

  sql = 'INSERT INTO ALIAS VALUES (?,?);'
  batch_write(con,sql,aliases)

  sql = 'SELECT COUNT(*) FROM ALIAS;'
  cur.execute(sql)
  print 'ALIASES:',cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_aliases ON ALIAS (alias);'
  cur.execute(sql)

  con.commit()


def squash_dups(snps):
  seen = set()
  for s in snps:
    name = s[0]
    if name not in seen:
      seen.add(name)
      yield s


def squash_dups2(aliases):
  seen = set()
  for a,n in aliases:
    key = a,n
    if key not in seen:
      seen.add(key)
      yield key


def squash_dups4(rows):
  seen = set()
  for row in rows:
    key = tuple(row[:4])
    if key not in seen:
      seen.add(key)
      yield row


def load_snps(con, snps):
  cur = con.cursor()

  try:
    cur.execute('DROP INDEX idx_snp_name;')
  except:
    pass

  try:
    cur.execute('DROP INDEX idx_snp_loc;')
  except:
    pass

  try:
    cur.execute('DROP TABLE SNP;')
  except:
    pass

  con.commit()

  cur.execute('''
    CREATE TABLE snp (id         INTEGER PRIMARY KEY,
                      name       TEXT,
                      chrom      TEXT,
                      start      INTEGER,
                      end        INTEGER,
                      strand     TEXT,
                      refAllele  TEXT,
                      alleles    TEXT,
                      vclass     TEXT,
                      func       TEXT,
                      weight     INTEGER);''')

  sql = 'INSERT INTO snp VALUES (NULL,?,?,?,?,?,?,?,?,?,?);'
  batch_write(con,sql,snps)

  sql = 'SELECT COUNT(*) FROM snp;'
  cur.execute(sql)
  print 'SNPS:',cur.fetchall()[0][0]

  con.commit()

  sql = 'CREATE INDEX idx_snp_name ON snp (name);'
  print sql
  cur.execute(sql)

  sql = 'CREATE INDEX idx_snp_loc ON snp (chrom,start,end);'
  print sql
  cur.execute(sql)

  con.commit()

  try:
    cur.execute('DROP TABLE snp_index;')
  except:
    pass

  print 'Creating SNP rtree index:'
  sql = '''CREATE VIRTUAL TABLE snp_index USING rtree_i32(id,start,end);'''
  cur.execute(sql)

  sql = '''
  INSERT INTO snp_index(id,start,end)
  SELECT id,start,end
  FROM   snp;
  '''
  cur.execute(sql)

  con.commit()


def get_goldenpath_dbsnp(filename):
  print 'LOADING UCSC DBSNP:',filename

  snps = table_reader(filename)

  SNPRecord = namedtuple('SNPRecord','bin chrom start end name score strand refNCBI refUCSC observed molType '
                                     'vclass valid avHet avHetSE func locType weight')

  for row in imap(SNPRecord._make,snps):
    name      = row.name
    chrom     = intern(row.chrom)
    start     = int(row.start)
    end       = int(row.end)
    strand    = row.strand
    refAllele = row.refUCSC.replace('-','')
    alleles   = row.observed.replace('-','')
    vclass    = intern(row.vclass)
    func      = intern(row.func if row.func!='unknown' else '')
    weight    = intern(row.weight)

    if len(refAllele)<20:
      refAllele = intern(refAllele)
    if len(alleles)<20:
      alleles = intern(alleles)

    yield name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight


def get_goldenpath_arrays(filename):
  import pysam

  reference = pysam.Fastafile('hg18/hg18.fa')

  print 'LOADING UCSC SNP ARRAY:',filename

  vclass    = ''
  func      = ''
  weight    = 1

  for row in table_reader(filename):
    name      = row[4]
    chrom     = row[1]
    start     = int(row[2])
    end       = int(row[3])
    strand    = row[6]

    try:
      refAllele = reference.fetch(chrom, start, end).upper()
      #assert len(refAllele)==(end-start)
      # refAllele is always relative to the + strand
      #if strand=='-':
      #  refAllele = rev_comp(refAllele)

    except IndexError:
      print 'Invalid region: %s:%d-%d' % (chrom, start+1, end)
      refAllele = None

    alleles = row[7].replace('-','')

    if len(row)==9:
      rsID = row[8]
      #print rsID,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
      yield rsID,chrom,start,end,strand,refAllele,alleles,vclass,func,weight

    #print name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
    yield name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight


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
    chrom = assay[chr_idx]
    position   = int(assay[loc_idx])-1
    strand     = strandmap[assay[assayid_idx].split('_')[-2]]
    alleles    = assay[snp_idx][1:-1]

    yield rs,chrom,position,strand,alleles


def get_cytobands(cytofile):
  print 'LOADING CYTOBANDS:',cytofile
  seen = set()
  for row in table_reader(cytofile):
    band  = row[3]
    chrom = row[0]
    start = tryint(row[1])
    stop  = tryint(row[2])
    color = row[4]

    if chrom.startswith('chr'):
      chrom = chrom[3:]
    if chrom.upper()=='MT':
      chrom = 'M'

    name  = chrom + band

    if name in seen:
      print name,row

    seen.add(name)

    yield name,chrom,start,stop,color


def load_cytobands(con,bands):
  cur = con.cursor()

  try:
    cur.execute('DROP TABLE CYTOBAND;')
  except:
    pass

  cur.execute('''
    CREATE TABLE CYTOBAND (band   TEXT,
                           chrom  TEXT,
                           start  INTEGER,
                           stop   INTEGER,
                           color  TEXT);''')

  sql = 'INSERT INTO CYTOBAND VALUES (?,?,?,?,?);'
  batch_write(con,sql,bands)

  sql = 'SELECT COUNT(*) FROM CYTOBAND;'
  cur.execute(sql)
  print 'BANDS:',cur.fetchall()[0][0]

  sql = 'CREATE INDEX idx_cytoband ON CYTOBAND (chrom,start,stop);'
  cur.execute(sql)

  con.commit()


def build_hg18():
  con = sqlite3.connect('out/genedb_hg18_snp130_rtree.db')

  sqlite_magic(con)

  if 1:
    genes = get_goldenpath_genes('hg18')
    load_genes(con,genes)

  if 1:
    bands = get_cytobands('hg18/cytoBand.txt.gz')
    load_cytobands(con,bands)

  if 1:
    seen = set()
    aliases = [ get_goldenpath_aliases(con,'hg18/kgAlias.txt.gz',seen),
                    get_entrez_aliases(con,'entrez/gene_info.gz',seen) ]
    aliases = chain.from_iterable(aliases)
    aliases = squash_dups2(aliases)
    load_gene_aliases(con,aliases)

  if 1:
    streams   = []
    streams  += [ get_goldenpath_dbsnp('hg18/snp130.txt.gz') ]
    streams  += [ get_goldenpath_arrays(a) for a in glob.glob('hg18/snpArray*.txt.gz') ]
    #streams += [ extract_illumina_snps(IlluminaManifest(m))
    #               for m in MANIFESTS ]

    with gcdisabled():
      snps = chain.from_iterable(streams)
      snps = squash_dups4(snps)
      load_snps(con,snps)


def build_hg19():
  con = sqlite3.connect('out/genedb_hg19_snp131_rtree.db')

  sqlite_magic(con)

  if 1:
    genes = get_goldenpath_genes('hg19')
    load_genes(con,genes)

  if 1:
    bands = get_cytobands('hg19/cytoBand.txt.gz')
    load_cytobands(con,bands)

  if 1:
    seen = set()
    aliases = [ get_goldenpath_aliases(con,'hg19/kgAlias.txt.gz',seen),
                    get_entrez_aliases(con,'entrez/gene_info.gz',seen) ]
    aliases = chain.from_iterable(aliases)
    aliases = squash_dups2(aliases)
    load_gene_aliases(con,aliases)

  if 1:
    streams   = []
    streams  += [ get_goldenpath_dbsnp('hg19/snp131.txt.gz') ]
    #streams  += [ get_goldenpath_arrays(a) for a in glob.glob('hg19/snpArray*.txt.gz') ]
    #streams += [ extract_illumina_snps(IlluminaManifest(m))
    #               for m in MANIFESTS ]

    with gcdisabled():
      snps = chain.from_iterable(streams)
      snps = squash_dups4(snps)
      load_snps(con,snps)


if __name__ == '__main__':
   main2()
