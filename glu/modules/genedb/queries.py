# -*- coding: utf-8 -*-

__abstract__  = 'Query genedb database'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from glu.lib.utils import namedtuple


CANONICAL_TRANSCRIPT = 1
CANONICAL_CHROMOSOME = 1<<1
CANONICAL_ASSEMBLY   = 1<<2
CANONICAL_ALL        = CANONICAL_TRANSCRIPT|CANONICAL_CHROMOSOME|CANONICAL_ASSEMBLY


def query_genes_by_name(con, gene, canonical=None, mapped=None):
  sql = '''
  SELECT   a.Alias,g.symbol,g.chrom,MIN(g.txStart) as start,MAX(g.txEnd) as end,g.strand,"GENE"
  FROM     alias a, gene g
  WHERE    %s
  GROUP BY g.symbol
  ORDER BY g.chrom,start,end,g.symbol,g.canonical DESC
  '''

  conditions = ['g.name = a.name']

  if canonical is not None:
    conditions.append('g.canonical & %d' % canonical)

  if mapped:
    conditions += [ 'g.txStart<>""', 'g.txEnd<>""' ]
  elif mapped is not None:
    conditions += [ 'g.txStart=""', 'g.txEnd=""' ]

  if '%' in gene:
    conditions.append('a.alias LIKE ?')
  else:
    conditions.append('a.alias = ?')

  sql = sql % '\n       AND '.join(conditions)

  cur = con.cursor()
  cur.execute(sql, (gene,))
  return cur.fetchall()


def query_gene_by_name(con,gene):
  genes = query_genes_by_name(con,gene)
  if not genes:
    raise KeyError('Cannot find gene "%s"' % gene)
  elif len(genes) > 1:
    raise KeyError('Gene not unique "%s"' % gene)
  return genes[0]


def query_gene_neighborhood(con,chrom,start,end,up,dn):
  sql = '''
  SELECT   symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand
  FROM     gene
  WHERE    chrom  = ?
    AND  ((strand = '+' AND txStart<? AND txEnd>?)
       OR (strand = '-' AND txStart<? AND txEnd>?))
  GROUP BY symbol
  ORDER BY chrom,start,end,symbol,canonical DESC
  '''
  cur = con.cursor()
  cur.execute(sql, (chrom,end+dn,start-up,end+up,start-dn))
  return cur.fetchall()


def query_genes_by_location(con,chr,start,end):
  sql = '''
  SELECT   symbol as alias,symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand,"GENE"
  SELECT   symbol,symbol,chrom,MIN(txStart) as start, MAX(txEnd) as end,strand
  FROM     gene
  WHERE    chrom = ?
    AND    txStart<? AND txEnd>?
  GROUP BY symbol
  ORDER BY chrom,start,end,symbol,canonical DESC
  '''
  cur = con.cursor()
  cur.execute(sql, (chr, loc))
  return cur.fetchall()


def query_cytoband_by_location(con,chrom,loc):
  sql = '''
  SELECT   band,start,stop,color
  FROM     cytoband
  WHERE    chrom = ?
    AND    ? BETWEEN start AND stop
  ORDER BY chrom,MIN(start,stop);
  '''
  if chrom.startswith('chr'):
    chrom = chrom[3:]
  cur = con.cursor()
  cur.execute(sql, (chrom, loc))
  return cur.fetchall()


def query_cytoband_by_name(con,name):
  sql1 = '''
  SELECT   chrom,start,stop,color
  FROM     cytoband
  WHERE    band = ?;
  '''

  cur = con.cursor()
  cur.execute(sql1, (name,))
  results = cur.fetchall()

  if len(results) == 1:
    return results[0]
  elif len(results) > 1:
    raise KeyError('Ambiguous cytoband "%s"' % name)

  sql2 = '''
  SELECT   chrom,MIN(start),MAX(stop),""
  FROM     cytoband
  WHERE    band LIKE (? || "%")
  GROUP BY chrom;
  '''

  if not name or ('p' not in name and 'q' not in name):
    return None

  cname = name
  if name[-1] in '0123456789':
    cname += '.'

  cur = con.cursor()
  cur.execute(sql2, (cname,))
  results = cur.fetchall()

  if not results:
    return None
  if len(results) == 1:
    return results[0]
  else:
    raise KeyError('Ambiguous cytoband "%s"' % name)


def query_contig_by_location(con,chrom,loc):
  sql = '''
  SELECT   name,start,stop
  FROM     contig
  WHERE    chrom = ?
    AND    ? BETWEEN start AND stop
  ORDER BY chrom,MIN(start,stop);
  '''
  if chrom.startswith('chr'):
    chrom = chrom[3:]
  cur = con.cursor()
  cur.execute(sql, (chrom, loc))
  return cur.fetchall()


def query_contig_by_name(con,name):
  sql1 = '''
  SELECT   chrom,start,stop
  FROM     contig
  WHERE    name = ?;
  '''

  cur = con.cursor()
  cur.execute(sql1, (name,))
  results = cur.fetchall()

  if not results:
    return None
  elif len(results) == 1:
    return results[0]
  elif len(results) > 1:
    raise KeyError('Ambiguous contig "%s"' % name)


def query_snps_by_name(con,name):
  sql = '''
  SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
  FROM     snp
  WHERE    name %s
  '''
  if '%' in name:
    sql = sql % 'LIKE ?'
  else:
    sql = sql % '= ?'

  cur = con.cursor()
  cur.execute(sql,(name,))
  return cur.fetchall()


def query_snp_by_name(con,name):
  snps = query_snps_by_name(con,name)
  if not snps:
    raise KeyError('Cannot find snp "%s"' % name)
  elif len(snps) > 1:
    raise KeyError('SNP not unique "%s"' % name)
  return snps[0]


def query_snps_by_location(con,chrom,start,end):
  sql = '''
  SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
  FROM     snp
  WHERE    chrom  = ?
    AND    start <= ?
    AND    end   >= ?
  ORDER BY start;
  '''
  cur = con.cursor()
  cur.execute(sql,(chrom,end,start))
  return cur.fetchall()


def query_snps_by_location_rtree(con,chrom,start,end):

  import time

  t0 = time.time()
  cur = con.cursor()
  cur.execute('SELECT id FROM snp_index WHERE start<=? AND end>=?', (end,start))
  ids = cur.fetchall()
  print '!!! Got %d ids in %.2f seconds' % (len(ids),time.time()-t0)

  if not ids:
    return []

  ids = ','.join([ str(i[0]) for i in ids ])

  sql = '''
    SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
    FROM     snp
    WHERE    chrom = '%s'
      AND    id in (%s)
    ORDER BY start;'''

  t0 = time.time()

  cur.execute(sql % (chrom,ids))
  rows = cur.fetchall()

  print '!!! Got %d SNPs in %.2f seconds' % (len(rows),time.time()-t0)

  return rows
