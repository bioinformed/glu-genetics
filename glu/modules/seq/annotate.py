# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                     import itemgetter
from   itertools                    import imap, count, izip, groupby
from   collections                  import defaultdict, namedtuple, OrderedDict

from   pysam                        import Fastafile
from   Bio.Seq                      import Seq
from   Bio.Data.CodonTable          import TranslationError

from   glu.lib.utils                import gcdisabled
from   glu.lib.fileutils            import table_writer, table_reader, autofile, hyphen
from   glu.lib.progressbar          import progress_loop
from   glu.lib.recordtype           import recordtype

from   glu.lib.seqlib.cgfvariants   import CGFVariants
from   glu.lib.seqlib.edits         import reduce_match
from   glu.lib.seqlib.intervaltree  import IntervalTree

from   glu.lib.genedb               import open_genedb
from   glu.lib.genedb.queries       import query_snps_by_location, query_cytoband_by_location


GeneRecord   = namedtuple('GeneRecord',   'id name symbol geneid mRNA protein canonical chrom strand '
                                          'txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds')
SNPRecord    = namedtuple('SNPRecord',    'name chrom start end strand refAllele alleles vclass func weight')
GeneFeature  = namedtuple('GeneFeature',  'chrom start end cds_index exon_num type')
GeneEvidence = recordtype('GeneEvidence', 'chrom cytoband ref_start ref_end intersect gene details '
                                          'func_class func_type ref_nuc var_nuc ref_aa var_aa')
VCFRecord    = recordtype('VCFRecord',    'chrom start end names ref var qual filter info format genos')


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
  return imap(GeneRecord._make, iter(cur))


def get_snps(con):
  cur = con.cursor()
  cur.execute('SELECT * FROM SNP;')

  make  = SNPRecord._make
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

  make = SNPRecord._make
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
      parts = gene_cache.pop(key)
    except KeyError:
      partlist = list(decode_gene(gene))

      parts = IntervalTree()
      for part in partlist:
        parts.insert(part.start,part.end,part)

      if len(gene_cache)>=300:
        gene_cache.popitem(0)

    # Add result to end of LRU
    gene_cache[key] = parts

    return parts


  def classify(self, chrom, ref_start, ref_end, variant, nsonly=False):
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
        evidence.append( ["5' of gene",gene,'','','',ref_nuc,var_nuc,'',''] )

      for gene in three_prime:
        evidence.append( ["3' of gene",gene,'','','',ref_nuc,var_nuc,'',''] )

    if not evidence:
      evidence.append( ['intergenic','','','','',ref_nuc,var_nuc,'',''] )

    evidence = group_evidence(evidence)
    cytoband = query_cytoband_by_location(self.con, chrom, ref_start)
    cytoband = ','.join(c[0] for c in cytoband)
    context  = [ chrom,cytoband,ref_start,ref_end ]

    if 0: # evidence:
      values = context+evidence[0]
      for f,v in zip(GeneEvidence.__slots__,values):
        print '%15s = %s' % (f,v)
      print

    evidence = [ GeneEvidence._make(context+e) for e in evidence ]

    return evidence


  def classify_feature(self, gene, ref_start, ref_end, ref_nuc, var_nuc):
    #print gene.name,gene.symbol,gene.strand,gene.txStart,gene.txEnd,gene.exonCount,gene.accession

    gene_parts = self.decode_gene(gene)

    intersect = defaultdict(list)
    for part in gene_parts.find_values(ref_start, ref_end):
      intersect[part.type].append(part)

    evidence = []

    parts    = set(intersect)
    mut_type = set()

    if gene.strand=='+':
      left,right = "5'","3'"
    else:
      left,right = "3'","5'"

    for splice in gene_parts.find_values(ref_start-10,ref_end+10):
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
      evidence.append([parts,gene,'','NON-SYNONYMOUS',mut_type,ref_nuc,var_nuc,'',''])
    else:
      evidence.append([parts,gene,'','',mut_type,ref_nuc,var_nuc,'',''])

    return evidence


  def classify_exonic_variant(self, gene, gene_parts, cds, ref_start, ref_end, ref_nuc, var_nuc):
    result = ['CDS',gene,'mRNA=%s:protein=%s:exon=%d:strand=%s' % \
                         (gene.mRNA,gene.protein,cds.exon_num,gene.strand)]

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
    for part in gene_parts:
      if part.type=='CDS':
        seq = Seq(self.reference.fetch(part.chrom,part.start,part.end))
        #assert len(seq)==(end-start)
        ref_cds_seq.append(seq)
        if part.cds_index<cds.cds_index:
          ref_var_start += len(seq)
        elif part.cds_index==cds.cds_index:
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

    return sorted(exact), sorted(inexact)


class VCFReader(object):
  def __init__(self, filename, hyphen=None):
    self.data = data    = table_reader(filename,hyphen=sys.stdin)

    self.metadata_order = metadata_order = []
    self.metadata       = metadata       = defaultdict(list)
    self.header         = None
    self.samples        = None

    for row in data:
      if not row:
        continue
      elif row[0].startswith('##'):
        meta,value = row[0].split('=',1)
        meta = meta[2:]

        if meta not in metadata:
          metadata_order.append(meta)

        metadata[meta].append(row)

      elif row[0].startswith('#'):
        self.header = header = list(row)
        header[0]   = header[0][1:]

        self.samples = [ (i,s.split('.')[0]) for i,s in enumerate(header[9:]) ]
        break
      else:
        raise ValueError('Invalid VCF file detected')


  def __iter__(self):
    strcache    = {}
    def intern(x,strcache=strcache.setdefault): return strcache(x,x)

    for row in self.data:
      chrom      = intern(row[0])
      end        = int(row[1])
      start      = end-1
      names      = row[2].split(',') if row[2]!='.' else []
      ref        = intern(row[3])
      var        = [ intern(v) for v in row[4].split(',') ]
      qual       = row[5]
      filter     = [ intern(f) for f in row[6].split(';') ] if row[6]!='.' else []
      info       = row[7].split(';')
      format     = intern(row[8])
      genos      = [ g.split(':') for g in row[9:] ] if len(row)>9 else None

      yield VCFRecord(chrom,start,end,names,ref,var,qual,filter,info,format,genos)


def load_kaviar(filename):
  kaviar = {}
  for row in table_reader(filename):
    chrom = row[0]
    if not chrom.startswith('chr'):
      chrom = 'chr'+chrom
    loc = int(row[1])

    for allelestuff in row[2:]:
      parts = allelestuff.split(':',1)
      if len(parts)!=2:
        continue
      allele,stuff = parts
      if allele=='rsids':
        continue
      stuff = stuff.strip().replace(', ',',')
      if not stuff.startswith('reference'):
        kaviar[chrom,loc,allele] = stuff

  return kaviar


def update_vcf_annotation(v, vs, cv, kaviar, options):
  if vs:
    # FIXME: Order genes and evidence consistently
    evidence   = list(vs.classify(v.chrom, v.start, v.end, v.var[0], nsonly=False))
    #v.names   = sorted(set(str(v) for e in evidence for v in e.varid_exact)|set(v.names))
    genes      = sorted(set(e.gene.symbol for e in evidence if e.gene and e.gene.symbol))
    geneids    = sorted(set(e.gene.geneid for e in evidence if e.gene and e.gene.geneid))
    location   = sorted(set(e.intersect   for e in evidence if e.intersect))
    function   = sorted(set(e.func_class or e.func_type  for e in evidence if e.func_class or e.func_type))
    nsevidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e.func_class ]
    nsinfo     = [ '%s:%s:%s->%s' % (e.gene.symbol,e.details,e.ref_aa,e.var_aa)
                   for e in nsevidence ]
    if not genes:
      v.filter.append('NonGenic')
    elif not nsevidence:
      v.filter.append('Synonymous')

    new_info = ['GENE_NAME=%s'             % (','.join(genes   )),
                'GENE_ID=%s'               % (','.join(str(g) for g in geneids)),
                'GENE_LOCATION=%s'         % (','.join(location)),
                'GENE_FUNCTION=%s'         % (','.join(function)),
                'GENE_FUNCTION_DETAILS=%s' % (','.join(nsinfo  )) ]

  if cv:
    cvinfo  = cv.score_and_classify(v.chrom,v.start,v.end,[v.ref,v.var[0]])
    if cvinfo.exact_vars:
      v.names = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)|set(v.names))
    inexact = ','.join(sorted(set(cvinfo.inexact_vars))) if cvinfo.inexact_vars else ''

    new_info.append('COMMON_SCORE=%.2f' % cvinfo.common_score)
    new_info.append('FUNCTION_SCORE=%d' % cvinfo.function_score)
    new_info.append('INEXACT_VARIANTS=%s' % inexact)

    if cvinfo.common_score>options.commonscore:
      v.filter.append('Common')

  if options.kaviar:
    kaviarkey = (v.chrom,v.start,v.var[0])
    if kaviarkey in kaviar:
      new_info.append('KAVIAR=%s' % kaviar[kaviarkey].replace(';',','))
      v.filter.append('Kaviar')
    else:
      new_info.append('KAVIAR=')

  if 'tgp' in v.names:
    v.names.remove('tgp')
    v.filter.append('1000G')

  if './.' in v.genos:
    filter.append('PartiallyCalled')

  # Remove any old fields that have been replaced by a new field
  new_info_fields = set(f.split('=',1)[0] for f in new_info)
  v.info          = new_info+[ f for f in v.info if f.split('=',1)[0] not in new_info_fields ]

  if len(v.filter)>1 and 'PASS' in v.filter:
    v.filter.remove('PASS')

  return v


def annotate_vcf(options):
  vs       = VariantAnnotator(options.genedb, options.reference)
  vcf      = VCFReader(options.variants,sys.stdin)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  out      = autofile(hyphen(options.output,sys.stdout),'w')
  kaviar   = load_kaviar(options.kaviar) if options.kaviar else None

  metadata = vcf.metadata
  metadata['FILTER'].append(['##FILTER=<PartiallyCalled,Description="Variant is not called for one or more samples">'])
  metadata['FILTER'].append(['##FILTER=<NonGenic,Description="Variant not in or near a gene">'])
  metadata['FILTER'].append(['##FILTER=<Synonymous,Description="Variant does not alter an amino-acid">'])

  if cv and options.commonscore:
    metadata['FILTER'].append(['##FILTER=<Common,Description="Variant is likely common with common score>%f">' % options.commonscore])
    metadata['FILTER'].append(['##FILTER=<1000G,Description="Variant was reported by 1000 Genomes project">'])

  if options.kaviar:
    metadata['FILTER'].append(['##FILTER=<Kaviar,Description="Variant appears in the Kaviar database">'])

  #metadata['FILTER'].append(['##FILTER=<NotDominant,Description="Variant does not fit dominant heritibility model">'])
  #metadata['FILTER'].append(['##FILTER=<NotRecessive,Description="Variant does not fit recessive heritibility model">'])

  metadata['INFO'].append(['##INFO=<ID=GENE_NAME,Number=.,Type=String,Description="Name of gene(s) containing variant">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_ID,Number=.,Type=String,Description="Entrez/LocusLink gene identifiers of genes containing variant">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_LOCATION,Number=.,Type=String,Description="Location of variant in gene(s)">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_FUNCTION,Number=.,Type=String,Description="Functional classification of variant for each gene and transcript">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_FUNCTION_DETAILS,Number=.,Type=String,Description="Functional details of variant for each gene and transcript">'])

  if cv:
    metadata['INFO'].append(['##INFO=<ID=COMMON_SCORE,Number=1,Type=Float,Description="Common score: maximum allele frequency in any population for rarest allele">'])
    metadata['INFO'].append(['##INFO=<ID=FUNCTION_SCORE,Number=1,Type=Int,Description="Function score: reported as functional variant in OMIM, dbSNP, or COSMIC">'])
    metadata['INFO'].append(['##INFO=<ID=INEXACT_VARIANTS,Number=.,Type=String,Description="Inexact variant matche">'])

  if options.kaviar:
    metadata['INFO'].append(['##INFO=<ID=KAVIAR,Number=.,Type=String,Description="Samples or datasets from Kaviar in which variant was found">'])

  for meta in vcf.metadata_order:
    for m in metadata[meta]:
      out.write('\t'.join(m))
      out.write('\n')

  out.write('#%s\n' % ('\t'.join(vcf.header)))

  for v in vcf:
    update_vcf_annotation(v, vs, cv, kaviar, options)

    # FORMAT: chrom start end names ref var filter info format genos
    row = [ v.chrom, str(v.end), ','.join(v.names) or '.', v.ref, ','.join(v.var), v.qual,
                                 ';'.join(sorted(v.filter)) or '.',
                                 ';'.join(v.info)] + [ ':'.join(g) for g in v.genos ]

    out.write('\t'.join(row))
    out.write('\n')


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('-f', '--format',   metavar='NAME', default='GLU',
                      help='File format (VCF, GLU)')
  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--cgfvariants',   metavar='NAME',
                      help='CGFvariant database annotation')
  parser.add_argument('--commonscore', metavar='T', type=float, default=0.05,
                      help='Annotate all variants with common score > T')
  parser.add_argument('--kaviar',   metavar='NAME',
                        help='Kaviar annotation (optional)')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()
  format  = options.format.upper()

  if format=='GLU':
    vs = VariantAnnotator(options.genedb, options.reference)
    cv = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
    filename = options.variants
    variants = table_reader(filename,hyphen=sys.stdin)
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
      #ref  = row[3]
      var   = row[4].split(',')[0]

      evidence = list(vs.classify(chrom, start, end, var, nsonly=True))
      evidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e[7] ]

      for e in evidence:
        out.writerow(e)

  elif format=='VCF':  # *** VCF, SNPs only ***
    annotate_vcf(options)
  else:
    raise ValueError('Unknown or Unsupported format specified: %s' % options.format)
