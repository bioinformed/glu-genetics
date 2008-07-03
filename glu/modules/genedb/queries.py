# -*- coding: utf-8 -*-
'''
File:          queries.py

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       Thr Aug  3 14:45:03 EDT 2006

Abstract:      Query genedb database

Compatibility: Python 2.5 and above

Requires:      No external dependencies, yet...

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


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
  SELECT   a.Alias,s.featureName,s.chromosome,s.orientation,s.chrStart,s.chrEnd
  FROM     alias a, gene s
  WHERE    s.geneID = a.geneID
    AND    a.Alias = ?
    AND    s.featureType='GENE'
    AND    s.submitGroup='reference'
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
    AND    featureType='GENE'
    AND    submitGroup='reference'
    AND    chromosome=?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  UNION
  SELECT   featureName,chromosome,orientation,chrStart,chrEnd
  FROM     gene
  WHERE    orientation = '-'
    AND    featureType='GENE'
    AND    submitGroup='reference'
    AND    chromosome=?
    AND    ? BETWEEN (chrStart - ?) AND (chrEnd + ?)
  ORDER BY chromosome,chrStart;
  '''
  cur = con.cursor()
  cur.execute(sql, (chromosome,location,up,dn,chromosome,location,dn,up))
  return cur.fetchall()
