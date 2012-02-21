# -*- coding: utf-8 -*-

__abstract__  = 'Query genedb database'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import re

from   glu.lib.utils import namedtuple


CANONICAL_TRANSCRIPT = 1
CANONICAL_CHROMOSOME = 1<<1
CANONICAL_ASSEMBLY   = 1<<2
CANONICAL_ALL        = CANONICAL_TRANSCRIPT|CANONICAL_CHROMOSOME|CANONICAL_ASSEMBLY

KGENOME_CHRMAP       = dict( (i,i) for i in range(1,23) )
KGENOME_CHRMAP.update({23:'X',24:'Y',25:'M'})

cyto_re = re.compile('(\d+|X|Y)(?:([p|q])(?:(\d+)(.\d+)?)?)?$')


def query_genes_by_name(con, gene, canonical_contig=True, canonical_transcript=None, mapped=None):
  sql = '''
  SELECT   a.Alias,g.symbol,g.chrom,MIN(g.txStart) as start,MAX(g.txEnd) as end,g.strand,"GENE"
  FROM     alias a, gene g
  WHERE    %s
  GROUP BY g.symbol
  ORDER BY g.chrom,start,end,g.symbol,g.canonical DESC
  '''

  conditions = ['g.name = a.name']

  if canonical_transcript is not None:
    conditions.append('g.canonical & %d' % canonical)

  if canonical_contig:
    conditions.append("g.chrom NOT LIKE '%!_%' ESCAPE '!'")

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
  genes = cur.fetchall()

  if len(genes)>1:
    new_genes = [ g for g in genes if g[1]==gene ]
    if new_genes:
      genes = new_genes

  if len(genes)>1:
    new_genes = [ g for g in genes if g[1].upper()==gene.upper() ]
    if new_genes:
      genes = new_genes

  return genes


def query_gene_by_name(con,gene,canonical_contig=True,canonical_transcript=None,mapped=None):
  genes = query_genes_by_name(con,gene,canonical_contig=canonical_contig,
                                       canonical_transcript=canonical_transcript,
                                       mapped=mapped)

  if not genes:
    raise KeyError('Cannot find gene "%s"' % gene)
  elif len(genes) > 1:
    raise KeyError('Gene not unique "%s"' % gene)

  return genes[0]


def query_gene_neighborhood(con,chrom,start,end,up,dn):
  cur = con.cursor()

  if 'rtree' not in con.filename:
    sql = '''
    SELECT   symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand
    FROM     gene
    WHERE    chrom  = ?
      AND  ((strand = '+' AND txStart<? AND txEnd>?)
         OR (strand = '-' AND txStart<? AND txEnd>?))
    GROUP BY symbol
    ORDER BY chrom,start,end,symbol,canonical DESC
    '''
    cur.execute(sql, (chrom,end+dn,start-up,end+up,start-dn))
  else:
    sql = '''
    SELECT   symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand
    FROM     gene
    WHERE    chrom  = ?
      AND    id in (SELECT id
                    FROM   gene_index
                    WHERE  txStart < ?
                      AND  txEnd   > ?)
      AND  ((strand = '+' AND txStart<? AND txEnd>?)
         OR (strand = '-' AND txStart<? AND txEnd>?))
    GROUP BY symbol
    ORDER BY chrom,start,end,symbol,canonical DESC
    '''
    d = max(dn,up)
    cur.execute(sql, (chrom,end+d,start-d,end+dn,start-up,end+up,start-dn))

  return cur.fetchall()


def query_genes_by_location(con,chrom,start,end):
  if 'rtree' not in con.filename:
    sql = '''
    SELECT   symbol as alias,symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand,"GENE"
    SELECT   symbol,symbol,chrom,MIN(txStart) as start, MAX(txEnd) as end,strand
    FROM     gene
    WHERE    chrom = ?
      AND    txStart<? AND txEnd>?
    GROUP BY symbol
    ORDER BY chrom,start,end,symbol,canonical DESC
    '''
  else:
    sql = '''
    SELECT   symbol as alias,symbol,chrom,MIN(txStart) as start,MAX(txEnd) as end,strand,"GENE"
    SELECT   symbol,symbol,chrom,MIN(txStart) as start, MAX(txEnd) as end,strand
    FROM     gene
    WHERE    chrom = ?
      AND    id in (SELECT id
                    FROM   gene_index
                    WHERE  txStart < ?
                      AND  txEnd   > ?)
    GROUP BY symbol
    ORDER BY chrom,start,end,symbol,canonical DESC
    '''

  cur = con.cursor()
  cur.execute(sql, (chrom, end, start))
  return cur.fetchall()


def split_cytoband(band):
  m = cyto_re.match(str(band))
  return m.groups() if m is not None else None


def cytoband_name(bands):
  if len(bands)==1:
    return bands[0][0]

  bands = [ split_cytoband(b[0]) for b in bands if b ]

  if not bands:
    return ''

  ref   = bands[0]
  match = 0

  for i in range(4):
    if not all(b[i] is not None and b[i]==ref[i] for b in bands):
      break
    match = i+1

  if not match:
    return ''

  return ''.join(ref[:match])


def query_cytoband_by_location(con,chrom,loc):
  sql = '''
  SELECT   band,start,stop,color
  FROM     cytoband
  WHERE    chrom = ?
    AND    ? BETWEEN start AND stop
  ORDER BY chrom,MIN(start,stop);
  '''
  if chrom is None:
    return '',[]
  if chrom.startswith('chr'):
    chrom = chrom[3:]
  cur = con.cursor()
  cur.execute(sql, (chrom, loc))

  bands = cur.fetchall()

  return cytoband_name(bands),bands


def query_cytobands_by_location(con,chrom,start,end):
  sql = '''
  SELECT   band,start,stop,color
  FROM     cytoband
  WHERE    chrom = ?
    AND    start<? AND stop>?
  ORDER BY chrom,start,stop;
  '''
  if chrom is None:
    return '',[]
  if chrom.startswith('chr'):
    chrom = chrom[3:]
  cur = con.cursor()
  cur.execute(sql, (chrom,end,start))

  bands = cur.fetchall()

  return cytoband_name(bands),bands


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


def query_snps_by_name(con,name,canonical=True,resolve_1kgenome=True):

  if resolve_1kgenome and name.startswith('SNP') and '-' in name:
    try:
      chrom,end = name[3:].split('-')
      chrom = 'chr%s' % KGENOME_CHRMAP[int(chrom)]
      end   = int(end)
      return [ (name,chrom,end-1,end,'+','?','?','?','?',None) ]
    except (ValueError,IndexError):
      pass

  sql = '''
  SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
  FROM     snp
  WHERE    name %s'''

  if '%' in name:
    sql = sql % 'LIKE ?'
  else:
    sql = sql % '= ?'

  if canonical:
    sql += "\n    AND    chrom NOT LIKE '%!_%' ESCAPE '!'"

  cur = con.cursor()
  cur.execute(sql,(name,))
  return cur.fetchall()


def query_snp_by_name(con,name,canonical=True):
  snps = query_snps_by_name(con,name,canonical)
  if not snps:
    raise KeyError('Cannot find snp "%s"' % name)
  elif len(snps) > 1:
    raise KeyError('SNP not unique "%s"' % name)
  return snps[0]


def query_snps_by_location(con,chrom,start,end):
  cur = con.cursor()

  if 'rtree' not in con.filename:
    sql = '''
    SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
    FROM     snp
    WHERE    chrom = ?
      AND    start < ?
      AND    end   > ?
    ORDER BY start,name;
    '''
  else:
    sql = '''
    SELECT   name,chrom,start,end,strand,refAllele,alleles,vclass,func,weight
    FROM     snp
    WHERE    chrom = ?
      AND    id in (SELECT id
                    FROM   snp_index
                    WHERE  start < ?
                      AND  end   > ?)
    ORDER BY start,name;
    '''

  cur.execute(sql,(chrom,end,start))
  return cur.fetchall()
