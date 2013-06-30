# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import re
import sys

from   operator                     import itemgetter
from   itertools                    import imap, count, izip, groupby
from   collections                  import defaultdict, namedtuple, OrderedDict

from   pysam                        import Fastafile

from   Bio.Seq                      import Seq
from   Bio.Data.CodonTable          import TranslationError

from   glu.lib.progressbar          import progress_loop
from   glu.lib.recordtype           import recordtype

from   glu.lib.seqlib.edits         import reduce_match
from   glu.lib.seqlib.intervaltree  import IntervalTree

from   glu.lib.genedb               import open_genedb


cyto_re = re.compile('(\d+|X|Y)(?:([p|q])(?:(\d+)(.\d+)?)?)?$')


BandRecord   = namedtuple('BandRecord',   'name chrom start end color')
GeneRecord   = namedtuple('GeneRecord',   'id name symbol geneid mRNA protein canonical chrom strand '
                                          'txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds '
                                          'category isRefSeq startcomplete endComplete nonsenseMediatedDecay '
                                          'strangeSplice atacIntrons selenocysteine')
GeneFeature  = namedtuple('GeneFeature',  'chrom start end cds_index exon_num type')
GeneEvidence = recordtype('GeneEvidence', 'chrom cytoband ref_start ref_end intersect gene details '
                                          'func func_class func_type ref_nuc var_nuc ref_aa var_aa')


def find_all(s, ss):
  start = 0
  n = len(ss)
  while True:
    start = s.find(ss, start)
    if start == -1: return
    yield start
    start += n


def common_suffix(s1,s2):
  min_len = min(len(s1), len(s2) )
  suffix  = []

  for c1,c2 in izip(reversed(s1),reversed(s2)):
    if c1==c2:
      suffix.append(c1)
    else:
      break

  suffix.reverse()
  return ''.join(suffix)


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
    left,right    = 'UTR5','UTR3'
    intron_offset = -1
  else:
    # Otherwise, genes transcribed on the reverse genomic strand are
    # numbered in increasing order and is bordered by the 3' UTR to the
    # "left" and 5' UTR to the "right" in genome orientation.
    exon_nums     = xrange(len(starts),0,-1)
    left,right    = 'UTR3','UTR5'
    intron_offset = 0

  if gene.category=='noncoding' or cds_start==cds_end:
    left = right = 'UTR'

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


def get_cytobands(con):
  cur = con.cursor()
  cur.execute('SELECT band,"chr"||chrom,start,stop,color FROM cytoband;')
  return imap(BandRecord._make, iter(cur))


def split_cytoband(band):
  m = cyto_re.match(str(band))
  return m.groups() if m is not None else None


def cytoband_name(bands):
  if not bands:
    return '?'
  elif len(bands)==1:
    return bands[0].name

  bands = [ split_cytoband(b.name) for b in bands if b ]

  if not bands:
    return '?'

  ref   = bands[0]
  match = 0

  for i in range(4):
    if not all(b[i] is not None and b[i]==ref[i] for b in bands):
      break
    match = i+1

  if not match:
    return '?'

  return ''.join(ref[:match])


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

    self.band_map = band_map = defaultdict(IntervalTree)
    for band in get_cytobands(self.con):
      band_map[band.chrom].insert(band.start,band.end,band)
      if band.chrom.startswith('chr') and band.chrom[3:] not in band_map:
        band_map[band.chrom[3:]] = band_map[band.chrom]

    trans = get_transcripts(self.con)
    trans = progress_loop(trans, label='Loading transcripts: ', units='transcripts')

    self.feature_map = feature_map = defaultdict(IntervalTree)
    for gene in trans:
      feature_map[gene.chrom].insert(gene.txStart,gene.txEnd,gene)

      if 0: # DEBUG
        parts = self.decode_gene(gene)
        for part in parts:
          if part.type not in ('intron','UTR5','UTR3','UTR') and '_' not in part.chrom:
            print '\t'.join(map(str,[part.chrom,part.start,part.end,gene.symbol]))

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


  def annotate(self, chrom, ref_start, ref_end, variant, nsonly=False):
    variant = variant.replace('-','')

    ref_nuc = self.reference.fetch(chrom,ref_start,ref_end).upper()
    var_nuc = variant.upper()

    evidence = []
    for feature in self.feature_map[chrom].find(ref_start, ref_end):
      evidence.extend( self.classify_feature(feature.value, ref_start, ref_end, ref_nuc, var_nuc) )

    #ns = any('NON-SYNONYMOUS' in e[3] for e in evidence)
    #if nsonly and not ns:
    #  return []

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
        evidence.append( ['UPSTREAM_GENE',gene,'',False,'','',ref_nuc,var_nuc,'',''] )

      for gene in three_prime:
        evidence.append( ['DOWNSTREAM_GENE',gene,'',False,'','',ref_nuc,var_nuc,'',''] )

    if not evidence:
      evidence.append( ['intergenic','','',False,'','',ref_nuc,var_nuc,'',''] )

    evidence = group_evidence(evidence)
    cytoband = cytoband_name(self.band_map[chrom].find_values(ref_start,ref_end))
    context  = [ chrom,cytoband,ref_start,ref_end ]

    if 0: # evidence:
      print
      for e in evidence:
        values = context+e
        for f,v in zip(GeneEvidence.__slots__,values):
          print '%15s = %s' % (f,v)
        print

    evidence = [ GeneEvidence._make(context+e) for e in evidence ]

    return evidence


  def classify_feature(self, gene, ref_start, ref_end, ref_nuc, var_nuc):
    gene_parts = self.decode_gene(gene)

    intersect = defaultdict(list)
    for part in gene_parts.find_values(ref_start, ref_end):
      intersect[part.type].append(part)

    evidence = []

    parts    = set(intersect)
    mut_type = set()

    for splice in gene_parts.find_values(ref_start-5,ref_end+5):
      if splice.type=='CDS' or 'UTR' in splice.type:
        if (0<splice.start-ref_end<=5) or (0<ref_start-splice.end<=5):
          mut_type.add('POSSIBLE_INTRONIC_SPLICE_VARIANT')

    parts    = ','.join(sorted(parts))
    mut_type = ','.join(sorted(mut_type))

    if len(intersect)==1 and len(intersect['CDS'])==1:
      e = self.classify_exonic_variant(gene, gene_parts, intersect['CDS'][0],
                                       ref_start, ref_end, ref_nuc, var_nuc)
      evidence.append(e)
    elif len(intersect['CDS']):
      evidence.append([parts,gene,'',True,'NON-SYNONYMOUS',mut_type,ref_nuc,var_nuc,'',''])
    elif mut_type:
      evidence.append([parts,gene,'',True,'PREDICTED-DISRUPT-TRANSCRIPT',mut_type,ref_nuc,var_nuc,'',''])
    elif len(intersect['UTR5'])+len(intersect['UTR3']):
      evidence.append([parts,gene,'',False,'UNKNOWN-UTR',mut_type,ref_nuc,var_nuc,'',''])
    elif len(intersect['intron']):
      evidence.append([parts,gene,'',False,'UNKNOWN-INTRONIC',mut_type,ref_nuc,var_nuc,'',''])
    else:
      evidence.append([parts,gene,'',False,'UNKNOWN-INTERGENIC',mut_type,ref_nuc,var_nuc,'',''])

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
      result += [False,'SYNONYMOUS','REFERENCE',ref_nuc,var_nuc,'','']
      return result

    ref_frame  = len(ref_nuc)%3
    var_frame  = len(var_nuc)%3
    frameshift = (len(ref_nuc)-len(var_nuc))%3

    if 0:
      print '  REF_FRAME: %d' % ref_frame
      print '  VAR_FRAME: %d' % var_frame

    mut_type = []

    if len(ref_nuc)==len(var_nuc):
      mut_type.append('SUBSTITUTION')
    elif len(ref_nuc)>len(var_nuc):
      mut_type.append('DELETION')
    else:
      mut_type.append('INSERTION')

    if exon_start<5:
      mut_type.append('POSSIBLE-SPLICE5')
    if cds.end-exon_end<5:
      mut_type.append('POSSIBLE-SPLICE3')

    if ref_frame!=var_frame:
      mut_type.append('FRAMESHIFT')
      mut_type = ','.join(sorted(mut_type))
      result += [True,'NON-SYNONYMOUS',mut_type,ref_nuc,var_nuc,'','']
      return result

    # FIXME: Request 100 bases beyond end of transcription
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
      mut_type.append('INVALID_TRANSLATION')
      mut_type = ','.join(sorted(mut_type))
      result += [True,'PRESUMED_NON-SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,'','']
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

      #assert len(ref_aa)

      result[-1] += ':aa=%d' % (aa_position+1)
      result += [False,'SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,str(ref_aa),str(ref_aa)]
      return result

    # Classify non-synonymous change by comparing AA sequences

    # Make sure ref protein doesn't appear to have spurious stops

    r = ref_cds_aa.rstrip('*')
    v = var_cds_aa.rstrip('*')

    ref_stop = r.find('*')
    var_stop = v.find('*')

    if ref_stop==-1:
      if var_stop!=-1 and not v.startswith(r):
        mut_type.append('PREMATURE_STOP')
      elif ref_cds_aa[-1]=='*' and var_cds_aa[-1]!='*':
        mut_type.append('LOSS_OF_STOP')

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
    result     += [True,'NON-SYNONYMOUS',mut_type,ref_cds_nuc,var_cds_nuc,str(ref_aa),str(var_aa)]

    return result
