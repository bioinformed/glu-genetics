# -*- coding: utf-8 -*-

__abstract__  = 'Query genedb database'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


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


def query_gene(con,name):
  sql = '''
  SELECT   a.Alias,s.featureName,s.chromosome,s.orientation,s.chrStart,s.chrEnd,s.featureType
  FROM     alias a, gene s
  WHERE    s.geneID = a.geneID
    AND    a.Alias = ?
  ORDER BY s.chromosome,MIN(s.chrStart,s.chrEnd);
  '''
  cur = con.cursor()
  cur.execute(sql, (name,))
  return cur.fetchall()


def query_snp(con,lname):
  sql = '''
  SELECT   lname,chromosome,location,strand
  FROM     snp
  WHERE    lname = ?
  '''
  cur = con.cursor()
  cur.execute(sql,(lname,))
  return cur.fetchall()


def query_gene_neighborhood(con,chromosome,location,up,dn):
  sql = '''
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '+'
    AND    chromosome=?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  UNION
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '-'
    AND    chromosome=?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  ORDER BY chromosome,chrStart;
  '''
  cur = con.cursor()
  cur.execute(sql, (chromosome,location,up,dn,chromosome,location,dn,up))
  return cur.fetchall()
