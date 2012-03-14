from   __future__ import division

import sys
import csv

from   operator           import itemgetter
from   itertools          import groupby
from   collections        import namedtuple

from   glu.lib.utils      import unique, is_str
from   glu.lib.fileutils  import autofile, hyphen
from   glu.lib.recordtype import recordtype


def tryint(s):
  try:
    return int(s)
  except ValueError:
    return s


def cga_base_reader(filename,hyfile=None,**kwargs):
  '''
  #ASSEMBLY_ID    GS00325-DNA_G01_37_1100-ASM
  #COSMIC COSMIC v48
  #DBSNP_BUILD    dbSNP build 131
  #FORMAT_VERSION 1.5
  #GENERATED_AT   2010-Dec-03 18:29:56.140933
  #GENERATED_BY   callannotate
  #GENE_ANNOTATIONS       NCBI build 37.1
  #GENOME_REFERENCE       NCBI build 37
  #SAMPLE GS00325-DNA_G01
  #SOFTWARE_VERSION       1.10.0.22
  #TYPE   GENE-ANNOTATION
  #PFAM_DATE      August 13, 2010

  >index  locus   allele  chromosome      begin   end     varType reference       call    xRef    geneId  mrnaAcc pro
  teinAcc symbol  orientation     component       componentIndex  hasCodingRegion impact  nucleotidePos   proteinPosannotationRefSequence    sampleSequence  genomeRefSequence       pfam
  '''

  bufsize = kwargs.get('bufsize',-1)

  cgafile = csv.reader(hyphen(autofile(filename,bufsize=bufsize),hyfile),dialect='excel-tab')

  extra = kwargs.get('extra',None)
  attrs = {}
  header = None

  while 1:
    row = next(cgafile)

    if not row:
      continue
    elif row[0].startswith('>'):
      row[0] = row[0][1:]
      header = row
      break
    elif row[0].startswith('#'):
      attrs[row[0][1:]] = row[1] if len(row)>1 else True

  if not header:
    raise ValueError('Unable to determine CGA header in file %s' % filename)

  if extra:
    if is_str(extra):
      extra = [extra]
    header.extend(extra)

  Row = recordtype('CGARow', header)

  if extra:
    e = [None]*len(extra)
    records = (Row(*(row+e)) for row in cgafile)
  else:
    records = (Row(*row)     for row in cgafile)

  return attrs,header,records


def cga_coverage_reader(records,**kwargs):
  '''
  ['chromosome', 'position', 'uniqueSequenceCoverage', 'weightSumSequenceCoverage', 'gcCorrectedCvg', 'avgNormalizedCoverage']
  CGARow(chromosome='chr1', position='60000', uniqueSequenceCoverage='3.744', weightSumSequenceCoverage='60.316', gcCorrectedCvg='63.508', avgNormalizedCoverage='54.6')
  '''

  str_cache = {'':None}
  cache_str = str_cache.setdefault

  for rec in records:
    rec.chromosome                = cache_str(rec.chromosome,rec.chromosome)
    if hasattr(rec,'position'):
      rec.position                = int(rec.position)
    else:
      rec.begin                   = int(rec.begin)
      rec.end                     = int(rec.end)
    rec.uniqueSequenceCoverage    = float(rec.uniqueSequenceCoverage)
    rec.weightSumSequenceCoverage = float(rec.weightSumSequenceCoverage)
    rec.gcCorrectedCvg            = float(rec.gcCorrectedCvg)
    rec.avgNormalizedCoverage     = float(rec.avgNormalizedCoverage)

    yield rec


def cga_dbsnpvar_reader(records,**kwargs):
  '''
  TYPE=DBSNP-TO-CGI
  >dbSnpId alleles chromosome begin end reference found exactMatch loci zygosity
           varTypeA hapA scoreA chromosomeA beginA endA
           varTypeB hapB scoreB chromosomeB beginB endB
  '''

  str_cache = {'':None}
  cache_str = str_cache.setdefault
  foundmap  = {'Y':True,'N':False,'?':None}

  def clean_allele(a):
    a = a.replace('-','')
    return cache_str(a,a)

  for rec in records:
    rec.dbSnpId                   = rec.dbSnpId.split(':')[1]
    rec.alleles                   = [ clean_allele(a) for a in rec.alleles.split('/') ]
    rec.chromosome                = cache_str(rec.chromosome,rec.chromosome)
    rec.begin                     = int(rec.begin)
    rec.end                       = int(rec.end)
    rec.reference                 = cache_str(rec.reference,rec.reference)
    rec.found                     = foundmap[rec.found]
    rec.exactMatch                = foundmap[rec.exactMatch]
    rec.loci                      = [ int(l) for l in rec.loci.split(',') ] if rec.loci else None
    rec.zygosity                  = cache_str(rec.zygosity,rec.zygosity)
    rec.varTypeA                  = cache_str(rec.varTypeA,rec.varTypeA)
    rec.scoreA                    = int(rec.scoreA) if rec.scoreA else None
    rec.chromosomeA               = cache_str(rec.chromosomeA,rec.chromosomeA)
    rec.beginA                    = int(rec.beginA)
    rec.endA                      = int(rec.endA)
    rec.varTypeB                  = cache_str(rec.varTypeB,rec.varTypeB)
    rec.scoreB                    = int(rec.scoreB) if rec.scoreB else None
    rec.chromosomeB               = cache_str(rec.chromosomeB,rec.chromosomeB)
    rec.beginB                    = int(rec.beginB) if rec.beginB else None
    rec.endB                      = int(rec.endB) if rec.endB else None

    yield rec


def cga_cnvsegments_reader(records,**kwargs):
  '''
  >chr    begin   end     avgNormalizedCvg        relativeCvg     calledPloidy    calledCNVType   ploidyScore     CNVTypeScore    overlappingGene knownCNV        repeats
  chr1    10000   177417  59.0    1.09    N       hypervariable   0       0
  '''
  str_cache = {'':None}
  cache_str = str_cache.setdefault

  for rec in records:
    rec.chr                       = cache_str(rec.chr,rec.chr)
    rec.begin                     = int(rec.begin)
    rec.end                       = int(rec.end)
    rec.avgNormalizedCvg          = float(rec.avgNormalizedCvg)
    rec.relativeCvg               = float(rec.relativeCvg)
    rec.calledPloidy              = int(rec.calledPloidy) if rec.calledPloidy!='N' else None
    rec.calledCNVType             = cache_str(rec.calledCNVType,rec.calledCNVType)
    rec.ploidyScore               = float(rec.ploidyScore)
    rec.CNVTypeScore              = float(rec.CNVTypeScore) if rec.CNVTypeScore!='NA' else None

    yield rec


def cga_tumor_cnvsegments_reader(records,**kwargs):
  '''
  >chr    begin   end     avgNormalizedCvg        relativeCvg     calledLevel     calledCNVType   levelScore
  CNVTypeScore
  chr1    10000   177417  61.2    0.93    0.951   NA      194     NA
  chr1    227417  267719  61.2    0.93    0.951   NA      194     NA
  chr1    317719  471368  61.2    0.93    0.951   NA      194     NA
  '''
  str_cache = {'':None}
  cache_str = str_cache.setdefault

  for rec in records:
    rec.chr                       = cache_str(rec.chr,rec.chr)
    rec.begin                     = int(rec.begin)
    rec.end                       = int(rec.end)
    rec.avgNormalizedCvg          = float(rec.avgNormalizedCvg)
    rec.relativeCvg               = float(rec.relativeCvg)
    rec.calledLevel               = float(rec.calledLevel) if rec.calledLevel!='N' else None
    rec.calledCNVType             = cache_str(rec.calledCNVType,rec.calledCNVType)
    rec.levelScore                = float(rec.levelScore)
    rec.CNVTypeScore              = float(rec.CNVTypeScore) if rec.CNVTypeScore!='NA' else None

    yield rec


def cga_junctions_reader(records,**kwargs):
  '''
  >Id LeftChr  LeftPosition  LeftStrand  LeftLength
      RightChr RightPosition RightStrand RightLength
      StrandConsistent Interchromosomal Distance
      DiscordantMatePairAlignments JunctionSequenceResolved
      TransitionSequence TransitionLength
      LeftRepeatClassification RightRepeatClassification LeftGenes RightGenes
      XRef DeletedTransposableElement KnownUnderrepresentedRepeat FrequencyInBaselineGenomeSet AssembledSequence
  4415	chr1	869376	+	247	chr1	870279	+	361	Y	N	903	13	Y	TACAGGT	7	Self chain;Tandem period 28	Self chain;Tandem period 28	NM_152486	NM_152486	rs70949527 (chr1:869471-870312)			0.75	gttacaggtgggcaggggaggcggctgcgtTACAGGTgggcgggggaggcggctgcgttacaggtgg	
  '''
  str_cache = {'':None}
  cache_str = str_cache.setdefault

  for rec in records:
    rec.Id                           = int(rec.Id)
    rec.LeftChr                      = cache_str(rec.LeftChr,rec.LeftChr)
    rec.LeftPosition                 = int(rec.LeftPosition)
    rec.LeftLength                   = int(rec.LeftLength)
    rec.RightChr                     = cache_str(rec.RightChr,rec.RightChr)
    rec.RightPosition                = int(rec.RightPosition)
    rec.RightLength                  = int(rec.RightLength)
    rec.FrequencyInBaselineGenomeSet = float(rec.FrequencyInBaselineGenomeSet)

    yield rec


def cga_mastervar_reader(records,**kwargs):
  '''
  CGARow(locus='1', ploidy='2', chromosome='chr1', begin='0', end='10000', zygosity='no-call', varType='no-ref',
         reference='=', allele1Seq='?', allele2Seq='?', allele1Score='', allele2Score='', allele1HapLink='',
         allele2HapLink='', xRef='', evidenceIntervalId='', allele1ReadCount='', allele2ReadCount='',
         referenceAlleleReadCount='', totalReadCount='', allele1Gene='', allele2Gene='', pfam='',
         miRBaseId='', repeatMasker='', segDupOverlap='', relativeCoverage='', calledPloidy='')

  for TUMORS:
    locus ploidy chromosome begin end zygosity varType reference allele1Seq
    allele2Seq allele1Score allele2Score allele1HapLink allele2HapLink xRef
    evidenceIntervalId allele1Re adCount allele2ReadCount
    referenceAlleleReadCount totalReadCount allele1Gene allele2Gene pfam
    miRBaseId repeatMasker segDupOverlap relativeCoverage
  '''
  generecord = recordtype('GeneRecord', 'geneId mrnaAcc symbol component impact')

  str_cache = {'':None}
  cache_str = str_cache.setdefault

  def parseGene(genetext):
    if not genetext:
      return None

    records = []

    for g in genetext.split(';'):
      gene = generecord(*g.split(':'))

      gene.geneId    = int(gene.geneId)
      gene.mrnaAcc   = gene.mrnaAcc or None
      gene.symbol    = gene.symbol  or None
      gene.component = cache_str(gene.component,gene.component)
      gene.impact    = cache_str(gene.impact,gene.impact)

      records.append(gene)

    return records or None

  skip_ref = kwargs.get('skip_ref',False)

  for record in records:
    if skip_ref and (record.zygosity=='no-call' or record.varType=='ref'):
      continue

    record.locus                    = int(record.locus)
    record.ploidy                   = int(record.ploidy)
    record.chromosome               = cache_str(record.chromosome,record.chromosome)
    record.begin                    = int(record.begin)
    record.end                      = int(record.end)
    record.zygosity                 = cache_str(record.zygosity,record.zygosity)
    record.varType                  = cache_str(record.varType,record.varType)
    record.reference                = cache_str(record.reference,record.reference)
    record.allele1Seq               = cache_str(record.allele1Seq,record.allele1Seq)
    record.allele2Seq               = cache_str(record.allele2Seq,record.allele2Seq)
    record.allele1Score             = int(record.allele1Score) if record.allele1Score else None
    record.allele2Score             = int(record.allele2Score) if record.allele2Score else None
    record.allele1HapLink           = record.allele1HapLink or None
    record.allele2HapLink           = record.allele2HapLink or None
    record.xRef                     = record.xRef or None
    record.evidenceIntervalId       = int(record.evidenceIntervalId) if record.evidenceIntervalId else None
    record.allele1ReadCount         = int(record.allele1ReadCount) if record.allele1ReadCount else None
    record.allele2ReadCount         = int(record.allele2ReadCount) if record.allele2ReadCount else None
    record.referenceAlleleReadCount = int(record.referenceAlleleReadCount) if record.referenceAlleleReadCount else None
    record.totalReadCount           = int(record.totalReadCount) if record.totalReadCount else None
    record.pfam                     = record.pfam or None
    record.miRBaseId                = record.miRBaseId or None
    record.repeatMasker             = record.repeatMasker or None
    record.segDupOverlap            = record.segDupOverlap or None
    record.relativeCoverage         = float(record.relativeCoverage) if record.relativeCoverage else None

    if hasattr(record,'calledPloidy'):
      record.calledPloidy           = cache_str(record.calledPloidy,record.calledPloidy) or None

    if record.allele1Gene==record.allele2Gene:
      record.allele1Gene = record.allele2Gene = parseGene(record.allele1Gene)
    else:
      record.allele1Gene            = parseGene(record.allele1Gene)
      record.allele2Gene            = parseGene(record.allele2Gene)

    # 653635:NR_024540.1:WASH5P:INTRON:UNKNOWN-INC

    yield record


def cga_mastervar_reader_v2(records,**kwargs):
  '''
  CGARow(locus='1', ploidy='2', chromosome='chr1', begin='0', end='10000', zygosity='no-call', varType='no-ref',
         reference='=', allele1Seq='?', allele2Seq='?', allele1Score='', allele2Score='', allele1HapLink='',
         allele2HapLink='', xRef='', evidenceIntervalId='', allele1ReadCount='', allele2ReadCount='',
         referenceAlleleReadCount='', totalReadCount='', allele1Gene='', allele2Gene='', pfam='',
         miRBaseId='', repeatMasker='', segDupOverlap='', relativeCoverage='', calledPloidy='')

  for TUMORS:
    locus ploidy chromosome begin end zygosity varType reference allele1Seq
    allele2Seq allele1Score allele2Score allele1HapLink allele2HapLink xRef
    evidenceIntervalId allele1Re adCount allele2ReadCount
    referenceAlleleReadCount totalReadCount allele1Gene allele2Gene pfam
    miRBaseId repeatMasker segDupOverlap relativeCoverage
  '''
  generecord = recordtype('GeneRecord', 'geneId mrnaAcc symbol component impact')

  str_cache = {'':None}
  cache_str = str_cache.setdefault

  skip_ref      = kwargs.get('skip_ref',False)
  skip_lq       = kwargs.get('skip_lq',False)
  cache_strings = kwargs.get('cache_strings',False)

  def parseGene(genetext):
    if not genetext:
      return None

    records = []

    for g in genetext.split(';'):
      gene = generecord(*g.split(':'))

      gene.geneId    = int(gene.geneId)
      gene.mrnaAcc   = gene.mrnaAcc or None
      gene.symbol    = gene.symbol  or None

      if cache_strings:
        gene.component = cache_str(gene.component,gene.component)
        gene.impact    = cache_str(gene.impact,gene.impact)

      records.append(gene)

    return records

  missing = ('','N')

  for record in records:
    if skip_ref and (record.zygosity=='no-call' or record.varType=='ref'):
      continue

    if skip_lq and skip_ref:
      if record.allele1Seq==record.reference and record.allele1VarQuality=='VQHIGH' \
                                             and record.allele2VarQuality=='VQLOW':
        continue

      if record.allele2Seq==record.reference and record.allele2VarQuality=='VQHIGH' \
                                             and record.allele1VarQuality=='VQLOW':
        continue

    record.locus                    = int(record.locus)
    record.ploidy                   = int(record.ploidy)
    record.begin                    = int(record.begin)
    record.end                      = int(record.end)
    record.allele1VarScoreVAF       = int(record.allele1VarScoreVAF) if record.allele1VarScoreVAF else None
    record.allele2VarScoreVAF       = int(record.allele2VarScoreVAF) if record.allele2VarScoreVAF else None
    record.allele1VarScoreEAF       = int(record.allele1VarScoreEAF) if record.allele1VarScoreEAF else None
    record.allele2VarScoreEAF       = int(record.allele2VarScoreEAF) if record.allele2VarScoreEAF else None
    record.allele1HapLink           = record.allele1HapLink or None
    record.allele2HapLink           = record.allele2HapLink or None
    record.allele1XRef              = record.allele1XRef    or None
    record.allele2XRef              = record.allele2XRef    or None
    record.evidenceIntervalId       = int(record.evidenceIntervalId) if record.evidenceIntervalId else None
    record.allele1ReadCount         = int(record.allele1ReadCount) if record.allele1ReadCount else None
    record.allele2ReadCount         = int(record.allele2ReadCount) if record.allele2ReadCount else None
    record.referenceAlleleReadCount = int(record.referenceAlleleReadCount) if record.referenceAlleleReadCount else None
    record.totalReadCount           = int(record.totalReadCount) if record.totalReadCount else None
    record.pfam                     = record.pfam           or None
    record.miRBaseId                = record.miRBaseId      or None
    record.repeatMasker             = record.repeatMasker   or None
    record.segDupOverlap            = record.segDupOverlap  or None
    record.relativeCoverageDiploid  = float(record.relativeCoverageDiploid) if record.relativeCoverageDiploid not in missing else None
    record.calledPloidy             = int(record.calledPloidy) if record.calledPloidy not in missing else None

    if cache_strings:
      record.chromosome               = cache_str(record.chromosome,record.chromosome)
      record.zygosity                 = cache_str(record.zygosity,record.zygosity)
      record.varType                  = cache_str(record.varType,record.varType)
      record.reference                = cache_str(record.reference,record.reference)
      record.allele1Seq               = cache_str(record.allele1Seq,record.allele1Seq)
      record.allele2Seq               = cache_str(record.allele2Seq,record.allele2Seq)
      record.allele1VarQuality        = cache_str(record.allele1VarQuality,record.allele1VarQuality)
      record.allele2VarQuality        = cache_str(record.allele2VarQuality,record.allele2VarQuality)

    if hasattr(record,'relativeCoverageNondiploid'):
      record.relativeCoverageNondiploid = float(record.relativeCoverageNondiploid) if record.relativeCoverageNondiploid not in missing else None

    if hasattr(record,'calledLevel'):
      record.calledLevel            = float(record.calledLevel) if record.calledLevel not in missing else None

    if record.allele1Gene==record.allele2Gene:
      record.allele1Gene = record.allele2Gene = parseGene(record.allele1Gene)
    else:
      record.allele1Gene            = parseGene(record.allele1Gene)
      record.allele2Gene            = parseGene(record.allele2Gene)

    if (record.allele1Seq != record.reference == record.allele2Seq) or \
       (record.allele2Seq and not record.allele1Seq):
      record.allele1Seq,record.allele2Seq                 = record.allele2Seq,record.allele1Seq
      record.allele1ReadCount,record.allele2ReadCount     = record.allele2ReadCount,record.allele1ReadCount
      record.allele1VarScoreVAF,record.allele2VarScoreVAF = record.allele2VarScoreVAF,record.allele1VarScoreVAF
      record.allele1VarScoreEAF,record.allele2VarScoreEAF = record.allele2VarScoreEAF,record.allele1VarScoreEAF
      record.allele1HapLink,record.allele2HapLink         = record.allele2HapLink,record.allele1HapLink
      record.allele1XRef,record.allele2XRef               = record.allele2XRef,record.allele1XRef
      record.allele1Gene,record.allele2Gene               = record.allele2Gene,record.allele1Gene

    if skip_lq:
      if record.allele1VarQuality=='VQLOW':
        record.allele1Seq  ='?'
        record.allele1Gene = None

      if record.allele2VarQuality=='VQLOW':
        record.allele2Seq  ='?'
        record.allele2Gene = None

    # 653635:NR_024540.1:WASH5P:INTRON:UNKNOWN-INC

    yield record


def cga_gene_reader(records,**kwargs):
  Variant = namedtuple('Variant', 'index locus chromosome begin end reference varType geneId '
                                  'mrnaAcc proteinAcc symbol orientation component componentIndex '
                                  'hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence '
                                  'pfam alleles')
  Allele  = namedtuple('Allele',  'allele_index varType nucSeq aaSeq impact xRef')

  varType_nocall = set(['no-call','no-call-rc','no-call-ri'])
  impact_nocall  = set(['UNKNOWN-VNC','UNKNOWN-INC','UNKNOWN-TR'])

  #impact_scores = {'NO-CHANGE'   :
  #                 'COMPATIBLE'  :
  #                 'MISSENSE'    :
  #                 'NONSENSE'    :
  #                 'NONSTOP'     :
  #                 'DELETE'      :
  #                 'INSERT'      :
  #                 'DELETE+'     :
  #                 'INSERT+'     :
  #                 'FRAMESHIFT'  :
  #                 'MISSTART'    :
  #                 'DISRUPT'     :
  #                 'UNKNOWN-VNC' :
  #                 'UNKNOWN-INC' :
  #                 'UNKNOWN-TR'  :

  def _allele_priority(allele):
    if allele.nucSeq==v.reference:
      return 0
    elif not allele.nucSeq.strip('N'):
      return 2
    else:
      return 1

  # Require index, locus, chromosome, begin, end, reference
  for index,vrecs in groupby(records,itemgetter(0,1,3,4,5)):
    vrecs   = list(vrecs)
    alleles = unique(Allele(   int(v.allele),
                            intern(v.varType),
                            intern(v.call),
                            intern(v.sampleSequence) or None,
                            intern(v.impact),
                            intern(v.xRef) or None)           for v in vrecs)

    v       = vrecs[0]
    alleles = sorted(alleles,key=_allele_priority)

    varType = set(a.varType for a in alleles) - varType_nocall

    if len(varType)>1:
      varType.discard('ref')

    if not varType:
      varType = None
    elif len(varType)==1:
      varType = intern(varType.pop())
    else:
      varType = 'complex'

    impact  = set(a.impact for a in alleles) - impact_nocall

    if len(impact)>1:
      impact.discard('NO-CHANGE')

    if not impact:
      impact = None
    elif len(impact)==1:
      impact = intern(impact.pop())
    else:
      impact = 'COMPLEX'

    variant = Variant(   int(v.index),
                         int(v.locus),
                      intern(v.chromosome),
                         int(v.begin),
                         int(v.end),
                      intern(v.reference),
                             varType,
                         int(v.geneId),
                             v.mrnaAcc if v.mrnaAcc!='-' else None,
                             v.proteinAcc,
                             v.symbol,
                      intern(v.orientation),
                      intern(v.component),
                         int(v.componentIndex) if v.componentIndex else None,
                             v.hasCodingRegion=='Y',
                             impact,            # impact
                      tryint(v.nucleotidePos) if v.nucleotidePos else None,
                      tryint(v.proteinPos)    if v.proteinPos    else None,
                      intern(v.annotationRefSequence),
                             v.pfam,
                             alleles)

    yield variant


def is_functional(v):
  return v.impact not in (None,'NO-CHANGE','COMPATIBLE')


def cga_reader(filename,hyphen=None,**kwargs):
  attrs,header,records = cga_base_reader(filename,hyphen,**kwargs)

  version = attrs.get('FORMAT_VERSION')

  if version is None:
    version = attrs.get('SOFTWARE_VERSION')

  v2      = version.startswith('2.')
  type    = attrs.get('TYPE')

  if type=='DEPTH-OF-COVERAGE':
    records = cga_coverage_reader(records,**kwargs)
  elif type=='VAR-OLPL':
    if v2:
      records = cga_mastervar_reader_v2(records,**kwargs)
    else:
      records = cga_mastervar_reader(records,**kwargs)
  elif type=='CNV-SEGMENTS':
    records = cga_cnvsegments_reader(records,**kwargs)
  elif type=='TUMOR-CNV-SEGMENTS':
    records = cga_tumor_cnvsegments_reader(records,**kwargs)
  elif type=='JUNCTIONS':
    records = cga_junctions_reader(records,**kwargs)
  elif type=='DBSNP-TO-CGI':
    records = cga_dbsnpvar_reader(records,**kwargs)
  elif type=='GENE-ANNOTATION':
    records = cga_gene_reader(records,**kwargs)
  else:
    raise ValueError('Unknown CGA type: %s' % type)

  return attrs,header,records


def main():
  for rec in cga_reader(sys.argv[1]):
    print rec


if __name__=='__main__':
  main()
