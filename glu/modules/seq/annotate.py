# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import time
import sys

from   operator                     import itemgetter
from   itertools                    import imap, count, izip, groupby
from   collections                  import defaultdict, namedtuple, OrderedDict

from   pysam                        import Fastafile
from   Bio.Seq                      import Seq
from   Bio.Data.CodonTable          import TranslationError

from   glu.lib.utils                import gcdisabled
from   glu.lib.fileutils            import table_writer, table_reader
from   glu.lib.progressbar          import progress_loop

from   glu.lib.seqlib.edits         import reduce_match
from   glu.lib.seqlib.intervaltree  import IntervalTree

from   glu.lib.genedb               import open_genedb
from   glu.lib.genedb.queries       import query_snps_by_location, query_cytoband_by_location


Gene = namedtuple('Gene', 'id name symbol geneid mRNA protein canonical chrom strand txStart txEnd '
                          'cdsStart cdsEnd exonCount exonStarts exonEnds')

SNP  = namedtuple('SNP',  'name chrom start end strand refAllele alleles vclass func weight')

GeneFeature = namedtuple('GeneFeature', 'chrom start end cds_index exon_num type')


def decode_gene(gene):
  # Decode CDS boundaries and exon start and end coordinates
  cds_start = int(gene.cdsStart)
  cds_end   = int(gene.cdsEnd)
  starts    = [ int(s) for s in gene.exonStarts.split(',') if s ]
  ends      = [ int(e) for e in gene.exonEnds.split(',')   if e ]

  if gene.strand=='+':
    # Genes transcribed on the forward genomic strand are numbered in
    # increasing order and is bordered by the 5' UTR to the "left" and 3'
    # UTR to the "right" in genome orientation.
    exon_nums     = count(1)
    left,right    = "5' UTR","3' UTR"
    intron_offset = -1
  else:
    # Otherwise, genes transcribed on the reverse genomic strand are
    # numbered in increasing order and is bordered by the 3' UTR to the
    # "left" and 5' UTR to the "right" in genome orientation.
    exon_nums     = xrange(len(starts),0,-1)
    left,right    = "3' UTR","5' UTR"
    intron_offset = 0

  last = gene.txStart
  cds_index = 0

  for exon_num,exon_start,exon_end in izip(exon_nums,starts,ends):
    # Computer the intersection between the current exon and the CDS. They
    # intersect on the half-open interval [c_start,c_end)if and only if
    # c_start<c_end
    c_start = max(exon_start,cds_start)
    c_end   = min(exon_end,  cds_end)

    # Output left UTR if the exon starts before the CDS
    if exon_start<cds_start:
      if last!=exon_start:
        yield GeneFeature(gene.chrom,last,exon_start,None,exon_num+intron_offset,'intron')

      last = min(exon_end,cds_start)
      yield GeneFeature(gene.chrom,exon_start,last,None,exon_num,left)

    # Output CDS if the exon intersects the CDS
    if c_start<c_end:
      if last!=c_start:
        yield GeneFeature(gene.chrom,last,c_start,None,exon_num+intron_offset,'intron')

      last = c_end
      yield GeneFeature(gene.chrom,c_start,c_end,cds_index,exon_num,'CDS')
      cds_index += 1

    # Output right UTR if the exon ends after the CDS
    if cds_end<exon_end:
      start = max(exon_start,cds_end)
      if last!=start:
        yield GeneFeature(gene.chrom,last,start,None,exon_num+intron_offset,'intron')

      last = exon_end
      yield GeneFeature(gene.chrom,start,exon_end,None,exon_num,right)

  if last!=gene.txEnd:
      yield GeneFeature(gene.chrom,last,gene.txEnd,None,exon_num,'intron')


def get_transcripts(con):
  cur = con.cursor()
  cur.execute('SELECT * FROM GENE;')
  return imap(Gene._make, iter(cur))


def get_snps(con):
  cur = con.cursor()
  #cur.execute('SELECT * FROM SNP WHERE func<>"";')
  cur.execute('SELECT * FROM SNP;')

  make  = SNP._make
  strcache    = {}
  allelecache = {}
  funccache   = {}

  def intern(x,strcache=strcache.setdefault): return strcache(x,x)

  empty = ()

  for row in cur:
    # name chrom start end strand refAllele alleles vclass func weight
    row    = list(row)
    row[1] = intern(row[1])
    row[4] = intern(row[4])
    row[5] = intern(row[5]) if len(row[5])<4 else row[5]

    alleles = allelecache.get(row[6])
    if alleles is not None:
      row[6] = alleles
    else:
      alleles = [a.replace('-','') for a in row[6].split('/')]
      if max(len(a) for a in alleles)>3:
        alleles = tuple(intern(a) if len(a)<3 else a for a in alleles)
        row[6]  = alleles
      else:
        alleles = tuple(intern(a) for a in alleles)
        row[6]  = allelecache[row[6]] = alleles

    row[7]  = intern(row[7])

    func = funccache.get(row[8] or '')
    if func is not None:
      row[8] = func
    else:
      func = tuple(intern(f) for f in row[8].split(',')) if row[8] else empty
      row[8] = funccache[row[8]] = func

    row[9]  = int(row[9])

    yield make(row)


def get_snps_interval(con, chrom, ref_start, ref_end):
  snps = query_snps_by_location(con, chrom, ref_start, ref_end)

  make = SNP._make
  for row in snps:
    # name chrom start end strand refAllele alleles vclass func weight
    row    = list(row)
    row[4] = row[4].replace('-','')
    row[5] = set(a.replace('-','') for a in row[5].split('/'))
    row[7] = set(row[7].split(',')) if row[7] else set()
    yield make(row)


def group_evidence(orig):
  key_func = itemgetter(0,1,3,4,5,6,7,8)
  orig.sort(key=key_func)
  new = []
  for key,values in groupby(orig, key_func):
    values = list(values)
    if len(values)>1:
      transcripts = ','.join(sorted(set(v[2] for v in values)))
      values[0][2] = transcripts
    new.append(values[0])

  return new


class VariantAnnotator(object):
  def __init__(self, gene_db, reference_fasta):
    self.reference  = Fastafile(reference_fasta)
    self.con        = open_genedb(gene_db)
    self.gene_cache = OrderedDict()

    trans = get_transcripts(self.con)
    trans = progress_loop(trans, label='Loading transcripts: ', units='transcripts')

    self.feature_map = feature_map = defaultdict(IntervalTree)
    for gene in trans:
      feature_map[gene.chrom].insert(gene.txStart,gene.txEnd,gene)

    sys.stderr.write('Loading complete.\n')


  def decode_gene(self,gene):
    gene_cache = self.gene_cache
    key = gene.id

    try:
      result = gene_cache.pop(key)
    except KeyError:
      gene_parts = list(decode_gene(gene))

      features = IntervalTree()
      for feature in gene_parts:
        features.insert(feature.start,feature.end,feature)

      result = gene_parts,features

      if len(gene_cache)>=200:
        gene_cache.popitem(0)

    # Add result to end of LRU
    gene_cache[key] = result

    return result


  def classify(self, chrom, ref_start, ref_end, variant, nsonly=False):
    #print chrom,ref_start,ref_end,variant

    variant = variant.replace('-','')

    ref_nuc = self.reference.fetch(chrom,ref_start,ref_end).upper()
    var_nuc = variant.upper()

    evidence = []
    for feature in self.feature_map[chrom].find(ref_start, ref_end):
      evidence.extend( self.classify_feature(feature.value, ref_start, ref_end, ref_nuc, var_nuc) )

    ns = any('NON-SYNONYMOUS' in e[3] for e in evidence)
    if nsonly and not ns:
      return []

    # If not in a gene, check to see if there are any genes nearby
    if not evidence:
      five_prime  = set()
      three_prime = set()

      for feature in self.feature_map[chrom].find(ref_start-2000, ref_end+2000):
        gene = feature.value
        if (0<ref_end-gene.txStart<=2000) ^ (gene.strand=='-'):
          five_prime.add(gene)
        else:
          three_prime.add(gene)

      for gene in five_prime:
        evidence.append( ["5' of gene",gene.symbol,gene.mRNA,'','',
                          ref_nuc,var_nuc,'',''] )

      for gene in three_prime:
        evidence.append( ["3' of gene",gene.symbol,gene.mRNA,'','',
                          ref_nuc,var_nuc,'',''] )

    if not evidence:
      evidence.append( ['intergenic','','','','',ref_nuc,var_nuc,'',''] )

    evidence = group_evidence(evidence)
    dbsnp    = self.get_dbsnp(chrom, ref_start, ref_end, ref_nuc, var_nuc)
    evidence = [ e+dbsnp for e in evidence ]

    cytoband = query_cytoband_by_location(self.con, chrom, ref_start)
    cytoband = ','.join(c[0] for c in cytoband)

    context  = [ chrom,cytoband,ref_start,ref_end ]
    evidence = [ context+e for e in evidence ]

    return evidence


  def classify_feature(self, gene, ref_start, ref_end, ref_nuc, var_nuc):
    #print gene.name,gene.symbol,gene.strand,gene.txStart,gene.txEnd,gene.exonCount,gene.accession

    gene_parts,features = self.decode_gene(gene)

    intersect = defaultdict(list)
    for inter in features.find(ref_start, ref_end):
      intersect[inter.value.type].append(inter.value)

    evidence = []
    accession = gene.protein or gene.mRNA

    parts    = set(intersect)
    mut_type = set()

    if gene.strand=='+':
      left,right = "5'","3'"
    else:
      left,right = "3'","5'"

    for splice in features.find(ref_start-10,ref_end+10):
      splice = splice.value
      if splice.type=='CDS' or 'UTR' in splice.type:
        if 0<splice.start-ref_end<=10:
          mut_type.add('POSSIBLE %s INTRONIC SPLICE VARIANT' % left)
        if 0<ref_start-splice.end<=10:
          mut_type.add('POSSIBLE %s INTRONIC SPLICE VARIANT' % right)

    parts    = ','.join(sorted(parts))
    mut_type = ','.join(sorted(mut_type))

    if len(intersect)==1 and len(intersect['CDS'])==1:
      e = self.classify_exonic_variant(gene, gene_parts, intersect['CDS'][0],
                                       ref_start, ref_end, ref_nuc, var_nuc)
      evidence.append(e)
    elif len(intersect['CDS'])>=1:
      evidence.append([parts,gene.symbol,accession,
                       'NON-SYNONYMOUS',mut_type,ref_nuc,var_nuc,'',''])
    else:
      evidence.append([parts,gene.symbol,accession,
                       '',mut_type,ref_nuc,var_nuc,'',''])

    return evidence


  def classify_exonic_variant(self, gene, gene_parts, cds, ref_start, ref_end, ref_nuc, var_nuc):
    accession = gene.mRNA or gene.protein
    result = ['CDS',gene.symbol,'%s:exon%d:strand=%s' % (accession,cds.exon_num,gene.strand)]
    exon_start = ref_start - cds.start
    exon_end   = ref_end   - cds.start

    # FIXME: Report ref and var nuc relative to gene strand

    var_nuc = var_nuc.upper()

    #print gene.chrom,ref_start,ref_end,ref_nuc,var_nuc
    #assert len(ref_nuc)==(ref_end-ref_start)

    if ref_nuc==var_nuc:
      result += ['SYNONYMOUS','FALSE POSITIVE',ref_nuc,var_nuc,'','']
      return result

    ref_frame = len(ref_nuc)%3
    var_frame = len(var_nuc)%3

    if 0:
      print '  REF_FRAME: %d' % ref_frame
      print '  VAR_FRAME: %d' % var_frame

    mut_type = []

    if exon_start<5:
      mut_type.append("POSSIBLE 5' EXONIC SPLICE VARIANT")
    if cds.end-exon_end<5:
      mut_type.append("POSSIBLE 3' EXONIC SPLICE VARIANT")

    if ref_frame!=var_frame:
      mut_type.append('FRAMESHIFT')
      mut_type = ','.join(sorted(mut_type))
      result += ['NON-SYNONYMOUS',mut_type,ref_nuc,var_nuc,'','']
      return result

    ref_var_start = 0
    ref_cds_seq   = []
    for chrom,start,end,cds_index,exon_num,label in gene_parts:
      if label=='CDS':
        seq = Seq(self.reference.fetch(chrom,start,end))
        #assert len(seq)==(end-start)
        ref_cds_seq.append(seq)
        if cds_index<cds.cds_index:
          ref_var_start += len(seq)
        elif cds_index==cds.cds_index:
          ref_var_start += exon_start

    #assert ref_nuc==str(ref_cds_seq[cds.cds_index][exon_start:exon_end]).upper()

    if 0:
      print '  CDS  : %d-%d' % (cds.start,cds.end)
      print '  VAR  : %d-%d' % (ref_start,ref_end)
      print '  LOCAL: %d-%d (size=%d)' % (exon_start,exon_end,len(ref_cds_seq[cds.cds_index]))

    var_cds_seq = ref_cds_seq[:]

    v = list(var_cds_seq[cds.cds_index])
    v[exon_start:exon_end] = list(var_nuc)
    var_cds_seq[cds.cds_index] = ''.join(v)

    ref_cds = Seq(''.join(str(s) for s in ref_cds_seq))
    var_cds = Seq(''.join(str(s) for s in var_cds_seq))

    if gene.strand=='-':
      ref_var_start = len(ref_cds)-ref_var_start-1
      ref_cds       = ref_cds.reverse_complement()
      var_cds       = var_cds.reverse_complement()
      ref_cds_nuc   = str(Seq(ref_nuc).reverse_complement())
      var_cds_nuc   = str(Seq(var_nuc).reverse_complement())
    else:
      ref_cds_nuc   = ref_nuc
      var_cds_nuc   = var_nuc

    try:
      ref_cds_aa = ref_cds.translate()
      var_cds_aa = var_cds.translate()
    except TranslationError:
      mut_type.append('INVALID TRANSLATION')
      mut_type = ','.join(sorted(mut_type))
      result += ['PRESUMED NON-SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,'','']
      return result

    ref_aa,var_aa,aa_position = reduce_match(str(ref_cds_aa),str(var_cds_aa))

    if not ref_aa and not var_aa:
      mut_type = ','.join(sorted(mut_type))

      codon_start  = ref_var_start-ref_var_start%3
      codon_end    = ref_var_start+len(ref_nuc)
      if codon_end%3:
        codon_end += 3-codon_end%3

      aa_position = codon_start//3
      ref_frame   = ref_cds[codon_start:codon_end]
      ref_aa      = ref_frame.translate()

      assert len(ref_aa)

      result[-1] += ':aa=%d' % (aa_position+1)
      result += ['SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,str(ref_aa),str(ref_aa)]
      return result

    # Classify non-synonymous change by comparing AA sequences

    ref_stop = ref_aa.find('*')
    var_stop = var_aa.find('*')

    if ref_stop==-1:
      ref_stop = len(ref_aa)
    if var_stop==-1:
      var_stop = len(var_aa)

    if var_stop<ref_stop:
      mut_type.append('PREMATURE STOP')
    elif var_stop>ref_stop:
      mut_type.append('LOSS OF STOP')

    if len(ref_aa)==len(var_aa):
      mut_type.append('SUBSTITUTION')
    elif len(ref_aa)>len(var_aa):
      mut_type.append('DELETION')
    else:
      mut_type.append('INSERTION')

    if 0:
      print '  REF_NUC:',ref_cds_nuc
      print '  VAR_NUC:',var_cds_nuc
      print '   REF_AA:',ref_aa
      print '   VAR_AA:',var_aa
      #print '  NUC_DIFF:',levenshtein_sequence(str(ref_cds),str(var_cds))
      #print '  AA_DIFF: ',levenshtein_sequence(str(ref_aa), str(var_aa) )

      ref_size = ref_end-ref_start
      cds_size = len(ref_cds)
      print '  CDS_SIZE=%d (%.1f codons)' % (cds_size,cds_size/3.0)
      print '  CDS SEQ=%s' % ref_cds

      assert not ref_cds or str(ref_cds[:3])=='ATG'

    mut_type = ','.join(sorted(mut_type))
    result[-1] += ':aa=%d' % (aa_position+1)
    result += ['NON-SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,str(ref_aa),str(var_aa)]

    return result


  def get_dbsnp(self, chrom, ref_start, ref_end, ref_nuc, var_nuc):
    snps     = list(get_snps_interval(self.con, chrom, ref_start, ref_end))

    var_comp = str(Seq(var_nuc).complement())
    var      = set([var_nuc,var_comp])

    exact    = set()
    inexact  = set()

    for snp in snps:
      alleles = set(str(snp.alleles).split('/'))
      match   = (snp.start==ref_start and snp.end==ref_end and var&alleles)
      #print snp.name,snp.start,ref_start,snp.end,ref_end,snp.refAllele,ref_nuc,var,alleles,var&alleles,match
      #print snp.name,snp.start==ref_start,snp.end==ref_end,snp.refAllele,ref_nuc,

      if match:
        exact.add(snp.name)
      else:
        inexact.add(snp.name)

    return [ ','.join(sorted(exact)), ','.join(sorted(inexact)) ]


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()

  vs = VariantAnnotator(options.genedb, options.reference)

  if 0:
    print vs.classify('chr1',1110293,1110294,'A')
    print vs.classify('chr1',960530,960532,'')
    print vs.classify('chr1',968624,968625,'A')
    print vs.classify('chr1',1839388,1839389,'G')
    print vs.classify('chr1',1840513,1840522,'')
    print vs.classify('chr1',11651480,11651481,'A')
    print vs.classify('chr17',38449841,38449842,'')
    print vs.classify('chr10',51219559,51219560,'')
    print vs.classify('chr10',123347961,123347962,'')
    print vs.classify('chr1',12830291,12830295,'CTTGG')
    print vs.classify('chr1',12842477,12842478,'T')
    print vs.classify('chr1',12842478,12842478,'TAA')
    print vs.classify('chr1',12862238,12862239,'CA')
    print vs.classify('chr4',106375635,106375636,'T')
    print vs.classify('chr4',106376408,106376410,'GGA')
    print vs.classify('chr8', 31616810,31616820,'')
    print vs.classify('chr8', 31615810,31615820,'')
    print vs.classify('chr8', 32743315,32743317,'')

    return

  if 0:
    sys.stderr.write('Loading SNPs...\n')
    snps = get_snps('out/genedb_hg18_snp130.db')

    out = table_writer(sys.stdout)
    for snp in snps:
      alleles = set(snp.alleles)
      alleles.discard(snp.refAllele)
      if len(alleles)==1:
        var_nuc  = next(iter(alleles))
        func     = ','.join(sorted(snp.func))
        evidence = vs.classify(snp.chrom,snp.start,snp.end,var_nuc)
        if not evidence:
          evidence.append(['UNKNOWN'])

        for row in evidence:
          out.writerow([snp.name,func]+row)

  if 0:
    from merge_diffs import load_diffs, VARIANT_HEADER
    filename = sys.argv[1]
    variants = load_diffs(filename)
    header   = list(VARIANT_HEADER)
    extra    = ['CHROM','REF_START','REF_END','INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                'REF_NUC_NEW','VAR_NUC_NEW','REF_AA_NEW','VAR_AA_NEW', 'dbSNP_exact','dbSNP_inexact']

    keep     = itemgetter(3,4,5,6,7,10,11,12,13)

    out      = table_writer(sys.stdout)
    out.writerow(header + list(keep(extra)))

    stats = defaultdict(int)

    for v in variants:
      evidence = list(vs.classify(v.chrom, v.start, v.end, v.var_nuc))

      newbler_ns = bool(v.ref_aa and (v.ref_aa!=v.var_aa or len(v.ref_nuc)!=len(v.var_nuc)))
      glu_ns     = 'NON-SYNONYMOUS' in [ e[6] for e in evidence ]
      if glu_ns:
        evidence = [ e for e in evidence if e[6]=='NON-SYNONYMOUS' ]

      newbler_dbsnp = bool(v.known)
      glu_dbsnp     = any(1 for e in evidence if e[12])

      #stats[newbler_dbsnp,glu_dbsnp] += 1
      #if newbler_dbsnp!=glu_dbsnp:
      for e in evidence:
        out.writerow(list(v)+list(keep(e)))

    #print >> sys.stderr,'Newbler_NS,GLU_NS,count'
    #for key in sorted(stats):
    #  print >> sys.stderr,key,stats[key]

  if 0:
    filename = options.variants
    variants = table_reader(filename)
    header   = next(variants)
    extra    = ['CHROM','REF_START','REF_END','INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                'REF_NUC_NEW','VAR_NUC_NEW','REF_AA_NEW','VAR_AA_NEW', 'dbSNP_exact','dbSNP_inexact']
    out      = table_writer(options.output)
    out.writerow(header + extra[3:])

    stats = defaultdict(int)

    for v in variants:
      evidence = list(vs.classify(v[0], int(v[1]), int(v[2]), v[7]))
      evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[6] ] or [[]]

      for e in evidence:
        out.writerow(v+e[3:])
  if 0:
    filename = options.variants
    variants = table_reader(filename)
    header   = next(variants)
    extra    = ['CHROM','REF_START','REF_END','INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS','FUNC_TYPE',
                'REF_NUC_NEW','VAR_NUC_NEW','REF_AA_NEW','VAR_AA_NEW', 'dbSNP_exact','dbSNP_inexact']
    out      = table_writer(options.output)
    out.writerow(header + extra[3:])

    stats = defaultdict(int)

    for v in variants:
      evidence = list(vs.classify(v[0], int(v[1]), int(v[2]), v[7]))
      evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[6] ] or [[]]

      for e in evidence:
        out.writerow(v+e[3:])

  if 0:
    filename = options.variants
    variants = table_reader(filename,hyphen=sys.stdin)
    header   = next(variants)
    out      = table_writer(options.output,hyphen=sys.stdout)


    header    = header+['CHROM','CYTOBAND','REF_START','REF_END','INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS',
                 'FUNC_TYPE','REF_NUC','VAR_NUC','REF_AA','VAR_AA', 'dbSNP_exact','dbSNP_inexact']
    out.writerow(header)

    for row in variants:
      if not row or row[0].startswith('#'):
        continue

      chrom = row[0]
      start = int(row[1])
      end   = int(row[2])
      ref   = row[3]
      var   = row[4].split(',')[0]

      evidence = list(vs.classify(chrom, start, end, var, nsonly=True))
      evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[7] ]

      for e in evidence:
        out.writerow(row+e)

  if 1:  # *** VCF, SNPs only ***
    filename = options.variants
    variants = table_reader(filename)
    out      = table_writer(options.output,hyphen=sys.stdout)


    header    = ['CHROM','CYTOBAND','REF_START','REF_END','INTERSECT','SYMBOL','ACCESSION','FUNC_CLASS',
                 'FUNC_TYPE','REF_NUC','VAR_NUC','REF_AA','VAR_AA', 'dbSNP_exact','dbSNP_inexact']
    out.writerow(header)

    for row in variants:
      if not row or row[0].startswith('#'):
        continue

      chrom = row[0]
      end   = int(row[1])
      start = end-1
      ref   = row[3]
      var   = row[4].split(',')[0]

      evidence = list(vs.classify(chrom, start, end, var, nsonly=True))
      evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[7] ]

      for e in evidence:
        out.writerow(e)


if __name__=='__main__':
  if 1:
    main()
  else:
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      prof.runcall(main)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
