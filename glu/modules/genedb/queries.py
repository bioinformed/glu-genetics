# -*- coding: utf-8 -*-

__abstract__  = 'Query genedb database'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


def query_genes_by_name(con,gene):
  sql = '''
  SELECT   a.Alias,g.featureName,g.chromosome,g.chrStart,g.chrEnd,g.orientation,g.featureType
  FROM     alias a, gene g
  WHERE    g.geneID = a.geneID
    AND    g.chrStart<>""
    AND    g.chrEnd<>""
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
  WHERE    chromosome  = ?
    AND  ((orientation = '+' AND ? BETWEEN (chrStart - ?) AND (chrEnd + ?))
       OR (orientation = '-' AND ? BETWEEN (chrStart - ?) AND (chrEnd + ?)))
  ORDER BY chromosome,chrStart;
  '''
  cur = con.cursor()
  cur.execute(sql, (chromosome,location,up,dn,location,dn,up))
  return cur.fetchall()


def query_genes_by_location(con,chr,loc):
  sql = '''
  SELECT   featureName,featureName,chromosome,chrStart,chrEnd,orientation,featureType
  FROM     gene
  WHERE    chromosome = ?
    AND    ? BETWEEN chrStart AND chrEnd
  ORDER BY chromosome,MIN(chrStart,chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, (chr, loc))
  return cur.fetchall()


def query_cytoband_by_location(con,chr,loc):
  sql = '''
  SELECT   band,start,stop,color
  FROM     cytoband
  WHERE    chromosome = ?
    AND    ? BETWEEN start AND stop
  ORDER BY chromosome,MIN(start,stop);
  '''
  cur = con.cursor()
  cur.execute(sql, (chr, loc))
  return cur.fetchall()


def query_cytoband_by_name(con,name):
  sql1 = '''
  SELECT   chromosome,start,stop,color
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
  SELECT   chromosome,MIN(start),MAX(stop),""
  FROM     cytoband
  WHERE    band LIKE (? || "%")
  GROUP BY chromosome;
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
