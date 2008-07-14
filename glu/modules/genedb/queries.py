# -*- coding: utf-8 -*-

__abstract__  = 'Query genedb database'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


def query_genes_by_name(con,gene):
  sql = '''
  SELECT   a.Alias,g.featureName,g.chromosome,g.chrStart,g.chrEnd,g.orientation,g.featureType
  FROM     alias a, gene g
  WHERE    g.geneID = a.geneID
    AND    a.Alias %s
  ORDER BY g.chromosome,MIN(g.chrStart,g.chrEnd);
  '''
  if '%' in gene:
    sql = sql % 'LIKE ?'
  else:
    sql = sql % '= ?'

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


def query_gene_neighborhood(con,chromosome,location,up,dn):
  sql = '''
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '+'
    AND    chromosome = ?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  UNION
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '-'
    AND    chromosome = ?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  ORDER BY chromosome,chrStart;
  '''
  cur = con.cursor()
  cur.execute(sql, (chromosome,location,up,dn,chromosome,location,dn,up))
  return cur.fetchall()


def query_genes_by_location(con,chr,loc):
  sql = '''
  SELECT   featureName,featureName,chromosome,chrStart,chrEnd,orientation,featureType
  FROM     gene
  WHERE    chromosome = %s
    AND    %d BETWEEN chrStart AND chrEnd
  ORDER BY chromosome,MIN(chrStart,chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, chr, loc)
  return cur.fetchall()


def query_snps_by_name(con,name):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     snp
  WHERE    lname %s
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


def query_snps_by_location(con,chr,start,end):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     snp
  WHERE    chromosome = ?
    AND    location BETWEEN ? AND ?
  ORDER BY location;
  '''
  cur = con.cursor()
  cur.execute(sql,(chr,start,end))
  return cur.fetchall()
