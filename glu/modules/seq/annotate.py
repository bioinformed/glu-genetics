# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                     import itemgetter

import pysam

from   glu.lib.utils                import unique
from   glu.lib.fileutils            import list_reader, autofile, hyphen
from   glu.lib.progressbar          import progress_loop

from   glu.lib.genedb.queries       import query_segdups, query_repeats

from   glu.lib.seqlib.vcf           import VCFReader, VCFWriter
from   glu.lib.seqlib.cga           import cga_reader
from   glu.lib.seqlib.vannotator    import VariantAnnotator
from   glu.lib.seqlib.cgfvariants   import CGFVariants
from   glu.lib.seqlib.kaviar        import kaviar_reader
from   glu.lib.seqlib.refvariants   import ReferenceVariants


class OrderedReader(object):
  def __init__(self, data, references, get_contig=itemgetter(0), get_loc=itemgetter(1)):
    self.data       = data
    self.refmap     = dict( (r,i) for i,r in enumerate(references) )
    self.current    = next(self.data,None)
    self.get_contig = get_contig
    self.get_loc    = get_loc

  def get(self, chrom, loc):
    current = self.current
    if current is None:
      return None

    try:
      data = self.data
      get_contig = self.get_contig

      current_chrom = get_contig(current)

      if chrom!=current_chrom:
        refmap   = self.refmap
        chromidx = refmap[chrom]

        # Stream is past current contig, do nothing...
        if chromidx<refmap[current_chrom]:
          return None

        while chromidx>refmap[current_chrom]:
          last_chrom = current_chrom
          while last_chrom==current_chrom:
            current       = next(data)
            current_chrom = get_contig(current)

      assert chrom==current_chrom

      get_loc     = self.get_loc
      current_loc = get_loc(current)

      while current_loc<loc:
        current       = next(data)
        current_chrom = get_contig(current)
        current_loc   = get_loc(current)

        if current_chrom!=chrom or current_loc>loc:
          self.current = current
          return None

      if current_loc>loc:
        self.current = current
        return None

      assert current_loc==loc
      results = []

      while 1:
        results.append(current)

        self.current = current = next(data,None)

        if current is None:
          break

        current_chrom = get_contig(current)
        current_loc   = get_loc(current)

        if current_chrom!=chrom or current_loc!=loc:
          break

      return results

    except StopIteration:
      self.current = None
      return None


def make_infomap(info):
  infomap = {}
  for inf in info:
    if '=' in inf:
      key,value = inf.split('=',1)
    else:
      key,value = inf,''

    infomap[key] = value
  return infomap


def polyphen2_code(c):
  if c=='?':
    return '.'
  elif len(c)==1:
    return c
  elif 'D' in c:
    return'D'
  elif 'd' in c:
    return'd'
  elif '?' in c:
    return'.'
  elif 'b' in c:
    return'b'
  else:
    return'.'


def aa_change(e):
  return '%s->%s' % (e.ref_aa,e.var_aa) if e.ref_aa or e.var_aa else ''


def update_vcf_annotation(v, vs, cv, esp, kaviar, refvars, polyphen2, options):
  new_info = []

  if vs:
    #a = vs.annotate(v.chrom, v.start, v.end, v.reference, v.var[0])
    #a = vs.annotate(v.chrom, v.start, v.end, v.var)
    #print '!!!',a
    #return v

    # FIXME: Order genes and evidence consistently
    evidence   = list(vs.annotate(v.chrom, v.start, v.end, v.var[0], nsonly=False))
    #v.names   = sorted(set(str(v) for e in evidence for v in e.varid_exact)|set(v.names))
    cytoband   = sorted(set(e.cytoband    for e in evidence if e.cytoband))
    genes      = sorted(set(e.gene.symbol for e in evidence if e.gene and e.gene.symbol))
    geneids    = sorted(set(e.gene.geneid for e in evidence if e.gene and e.gene.geneid))
    location   = sorted(set(e.intersect   for e in evidence if e.intersect))
    function   = sorted(set(e.func_type or e.func_class  for e in evidence if e.func))
    nsevidence = [ e for e in evidence if e.func ]
    nsinfo     = ( '%s:%s:%s:%s' % (e.gene.symbol,e.func_type,e.details,aa_change(e))
                   for e in nsevidence )
    nsinfo     = list(unique(nsinfo))

    segdups    = query_segdups(vs.con, v.chrom, v.start, v.end)
    repeats    = query_repeats(vs.con, v.chrom, v.start, v.end)

    while 'PASS' in v.filter:
      v.filter.remove('PASS')

    if not genes:
      v.filter.append('Intergenic')

    if not nsevidence:
      v.filter.append('NPF')

    if cytoband:
      new_info.append('CYTOBAND=%s' % (','.join(cytoband)))

    if genes:
      new_info.append('GENE_NAME=%s'             % (','.join(genes   )))

    if geneids:
      new_info.append('GENE_ID=%s'               % (':'.join(str(g) for g in geneids)))

    if location:
      new_info.append('GENE_LOCATION=%s'         % (','.join(location)))

    if function:
      new_info.append('GENE_FUNCTION=%s'         % (','.join(function)))

    if nsinfo:
      new_info.append('GENE_FUNCTION_DETAILS=%s' % (','.join(nsinfo  )))

    if segdups:
      segdup_info = ('chr%s:%s-%s:%.2f' % (s.other_chrom,s.other_start,s.other_stop,s.matchFraction*100) for s in segdups)
      segdup_info = ','.join(segdup_info)
      v.filter.append('SegDup')
      new_info.append('SEGDUP_COUNT=%d' % len(segdups))
      new_info.append('SEGDUP_INFO=%s' % segdup_info)

    if repeats:
      repeat_info = ('%s:%s:%s' % (r.repeatName,r.repeatClass,r.repeatFamily) for r in repeats)
      repeat_info = ','.join(repeat_info)
      v.filter.append('Repeat')
      new_info.append('REPEAT_COUNT=%d' % len(repeats))
      new_info.append('REPEAT_INFO=%s' % repeat_info)

  if cv:
    cvinfo  = cv.score_and_classify(v.chrom,v.start,v.end,[v.ref,v.var[0]])
    if cvinfo.exact_vars:
      v.names = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)|set(v.names))

    if cvinfo.exact_vars or cvinfo.common_score>0:
      new_info.append('COMMON_SCORE=%.2f' % cvinfo.common_score)

    if cvinfo.function_info:
      function_info = ','.join(cvinfo.function_info)
      new_info.append('FUNCTION_INFO=%s' % function_info)

    if cvinfo.inexact_vars:
      inexact = ','.join(sorted(set(cvinfo.inexact_vars)))
      new_info.append('INEXACT_VARIANTS=%s' % inexact)

    if cvinfo.common_score>options.commonscore:
      v.filter.append('Common')

  if esp:
    chrom = v.chrom
    vars  = set(v.var)
    if chrom.startswith('chr'):
      chrom = chrom[3:]

    try:
      esp_info = esp.fetch(chrom,v.start,v.end)
      esp_info = [ e for e in esp_info if e.start==v.start and e.end==v.end and vars.intersection(e.var) ]
    except ValueError:
      esp_info = None

    if esp_info:
      gtc = maf_ea = maf_aa = 0

      for einfo in esp_info:
        infomap  = make_infomap(einfo.info)
        if 'MAF' in infomap:
          mafs     = map(float,infomap['MAF'].split(','))
          maf_ea   = max(maf_ea,mafs[0]/100)
          maf_aa   = max(maf_aa,mafs[1]/100)

        if 'GTC' in infomap:
          gtcs     = map(int,infomap['GTC'].split(','))
          gtc      = max(gtc,gtcs[0]+gtcs[1])

      maf = max(maf_ea,maf_aa)

      v.filter.append('ESP')
      if maf>options.commonscore:
        v.filter.append('ESPCommon')

      new_info.append('ESP_COUNT=%d' % gtc)
      new_info.append('ESP_MAF_EA=%.4f' % maf_ea)
      new_info.append('ESP_MAF_AA=%.4f' % maf_aa)

  if kaviar:
    kinfo = kaviar.get(v.chrom,v.start) or []
    ktext = [ stuff.replace(';',',') for chrom,loc,mallele,maf,allele,stuff in kinfo if allele in v.var ]

    if ktext:
      v.filter.append('Kaviar')

      mallele = kinfo[0][2]
      maf     = kinfo[0][3]

      if mallele and mallele in v.var:
        new_info.append('KAVIAR_MAF=%.2f' % maf)

        if maf>options.commonscore:
          v.filter.append('KaviarCommon')

      new_info.append('KAVIAR_NAMES=%s' % ','.join(ktext))


  if refvars:
    ingroup,outgroup = refvars.get(v.chrom,v.start,v.end,v.var) if refvars else ([],[])

    new_info.append('REFVAR_INGROUP_COUNT=%d'  % len(ingroup))
    new_info.append('REFVAR_OUTGROUP_COUNT=%d' % len(outgroup))

    if ingroup:
      ingroup  = ','.join(ingroup)
      new_info.append('REFVAR_INGROUP_NAMES=%s'  % ingroup)

    if outgroup:
      outgroup = ','.join(outgroup)
      new_info.append('REFVAR_OUTGROUP_NAMES=%s' % outgroup)
      v.filter.append('RefVar')

  if vs and polyphen2 and v.end-v.start==1 and 'CDS' in location:
    pmap  = {'A':0,'C':1,'G':2,'T':3}

    try:
      pvars = [ p.rstrip().split('\t') for p in polyphen2.fetch(v.chrom,v.start,v.end) ]
    except ValueError:
      pvars = None

    if pvars:
      hdivs = []
      hvars = []

      for a in v.var:
        if a in pmap:
          i = pmap[a]

          hdiv = hvar = ''
          for p in pvars:
            hdiv += p[3][i]
            hvar += p[4][i]

          hdiv = polyphen2_code(hdiv)
          hvar = polyphen2_code(hvar)
        else:
          hdiv = '.'
          hvar = '.'

        assert hdiv!='r' and hvar!='r'

        hdivs.append(hdiv)
        hvars.append(hvar)

      if hdivs.count('.')!=len(hdivs):
        new_info.append('POLYPHEN2_HDIV=%s' % (','.join(hdivs)))
      if hvars.count('.')!=len(hvars):
        new_info.append('POLYPHEN2_HVAR=%s' % (','.join(hvars)))

  if 'tgp' in v.names:
    v.names.remove('tgp')
    v.filter.append('1000G')

  if not v.ref or v.var==['']:
    v.filter.append('Indel')

  # Remove any old fields that have been replaced by a new field
  new_info_fields = set(f.split('=',1)[0] for f in new_info)
  v.info          = new_info+[ f for f in v.info if f.split('=',1)[0] not in new_info_fields ]

  v.filter = v.filter or ['PASS']

  return v


def annotate_vcf(options):
  vs       = VariantAnnotator(options.genedb, options.reference)
  vcf      = VCFReader(options.variants,sys.stdin)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  esp      = VCFReader(options.esp) if options.esp else None
  polyphen2= pysam.Tabixfile(options.polyphen2) if options.polyphen2 else None

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars  = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  metadata = vcf.metadata

  metadata.setdefault('FILTER',[])
  metadata.setdefault('INFO',[])

  metadata['FILTER'].append('##FILTER=<ID=PartiallyCalled,Description="Variant is not called for one or more samples">')
  metadata['FILTER'].append('##FILTER=<ID=Indel,Description="Variant is an insertion or deletion">')
  metadata['FILTER'].append('##FILTER=<ID=Intergenic,Description="Variant not in or near a gene">')
  metadata['FILTER'].append('##FILTER=<ID=NPF,Description="Variant is not predicted to alter a protein">')
  metadata['FILTER'].append('##FILTER=<ID=SegDup,Description="Variant occurs in a segmentally duplicated region">')
  metadata['FILTER'].append('##FILTER=<ID=Repeat,Description="Variant occurs in a repetitive or low-complexity region">')

  metadata['INFO'].append('##INFO=<ID=CYTOBAND,Number=.,Type=String,Description="Name of cytoband(s) containing variant">')
  metadata['INFO'].append('##INFO=<ID=GENE_NAME,Number=.,Type=String,Description="Name of gene(s) containing variant">')
  metadata['INFO'].append('##INFO=<ID=GENE_ID,Number=.,Type=String,Description="Entrez/LocusLink gene identifiers of genes containing variant">')
  metadata['INFO'].append('##INFO=<ID=GENE_LOCATION,Number=.,Type=String,Description="Location of variant in gene(s)">')
  metadata['INFO'].append('##INFO=<ID=GENE_FUNCTION,Number=.,Type=String,Description="Functional classification of variant for each gene and transcript">')
  metadata['INFO'].append('##INFO=<ID=GENE_FUNCTION_DETAILS,Number=.,Type=String,Description="Functional details of variant for each gene and transcript">')
  metadata['INFO'].append('##INFO=<ID=SEGDUP_COUNT,Number=1,Type=Integer,Description="Number of segmental duplications that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=SEGDUP_INFO,Number=.,Type=String,Description="Details of segmental duplications that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=REPEAT_COUNT,Number=1,Type=Integer,Description="Number of repetitive or low complexity elements that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=REPEAT_INFO,Number=.,Type=String,Description="Details of repetitive or low complexity elements that overlap variant locus">')

  if cv:
    metadata['FILTER'].append('##FILTER=<ID=Common,Description="Variant is likely common with common score>%f">' % options.commonscore)
    metadata['FILTER'].append('##FILTER=<ID=1000G,Description="Variant was reported by 1000 Genomes project">')
    metadata['INFO'  ].append('##INFO=<ID=COMMON_SCORE,Number=1,Type=Float,Description="Common score: maximum allele frequency in any population for rarest allele">')
    metadata['INFO'  ].append('##INFO=<ID=FUNCTION_INFO,Number=.,Type=String,Description="Annotated as function by OMIM, dbSNP, or COSMIC">')
    metadata['INFO'  ].append('##INFO=<ID=INEXACT_VARIANTS,Number=.,Type=String,Description="Observed variant matches inexactly: it has different alleles or overlaps observed">')

  if esp:
    metadata['FILTER'].append('##FILTER=<ID=ESP,Description="Variant appears in the UW ESP database">')
    metadata['FILTER'].append('##FILTER=<ID=ESPCommon,Description="Variant appears in the UW ESP database and appears to be common with MAF>%f">' % options.commonscore)
    metadata['INFO'  ].append('##INFO=<ID=ESP_COUNT,Number=1,Type=Integer,Description="Count of individuals with one or more variant alleles in UW ESP database">')
    metadata['INFO'  ].append('##INFO=<ID=ESP_MAF_EA,Number=1,Type=Float,Description="Minor allele frequency in European-Americans according to UW ESP database">')
    metadata['INFO'  ].append('##INFO=<ID=ESP_MAF_AA,Number=1,Type=Float,Description="Minor allele frequency in African-Americans according to UW ESP database">')

  if kaviar:
    metadata['FILTER'].append('##FILTER=<ID=Kaviar,Description="Variant appears in the Kaviar database">')
    metadata['FILTER'].append('##FILTER=<ID=KaviarCommon,Description="Variant appears in the Kaviar database and appears to be common with MAF>%f">' % options.commonscore)
    metadata['INFO'  ].append('##INFO=<ID=KAVIAR_MAF,Number=1,Type=Float,Description="Minor allele frequency according to Kaviar database">')
    metadata['INFO'  ].append('##INFO=<ID=KAVIAR_NAMES,Number=.,Type=String,Description="Samples or datasets from Kaviar in which variant was found">')

  if refvars:
    metadata['FILTER'].append('##FILTER=<ID=RefVar,Description="Variant appears in the local reference variant list">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_INGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in intra-group samples">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_OUTGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in extra-group samples">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_INGROUP_NAMES,Number=.,Type=String,Description="Intra-group samples in which Variant is present">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_OUTGROUP_NAMES,Number=.,Type=String,Description="Extra-group samples in which Variant is present">')

  if polyphen2:
    metadata['INFO'  ].append('##INFO=<ID=POLYPHEN2_HDIV,Number=.,Type=String,Description="Polyphen2 HDIV prediction code for each variant SNV allele (b=benign, d=possibly damaging, D=probably damaging, .=unknown)">')
    metadata['INFO'  ].append('##INFO=<ID=POLYPHEN2_HVAR,Number=.,Type=String,Description="Polyphen2 HVAR prediction code for each variant SNV allele (b=benign, d=possibly damaging, D=probably damaging, .=unknown)">')

  out = VCFWriter(options.output, metadata, vcf.samples, options.reference)

  for v in vcf:
    update_vcf_annotation(v, vs, cv, esp, kaviar, refvars, polyphen2, options)

    out.write_locus(v)


def valid_allele(a):
  return a is not None and a!='=' and 'N' not in a and '?' not in a


def annotate_mastervar(options):
  vs       = VariantAnnotator(options.genedb, options.reference)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  out      = autofile(hyphen(options.output,sys.stdout),'w')

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars    = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  attrs,header,vars = cga_reader(options.variants,sys.stdin)

  for var in vars:
    varid_exact     = []
    varid_inexact   = []
    common_score    = 0
    function_score  = 0
    kaviar_maf      = 0
    kaviar_details  = []
    geneinfo        = []

    allele1         = var.allele1Seq if valid_allele(var.allele1Seq) else None
    allele2         = var.allele2Seq if valid_allele(var.allele2Seq) else None

    alleles = []
    if allele1 is not None:
      alleles.append(allele1)
    if allele2 is not None:
      alleles.append(allele2)

    if len(alleles)!=2:
      continue

    if vs:
      evidence1     = list(vs.classify(var.chromosome, var.begin, var.end, allele1)) if allele1 is not None else []
      evidence2     = list(vs.classify(var.chromosome, var.begin, var.end, allele2)) if allele2 is not None else []
      evidence      = evidence1+evidence2

      geneinfo      = [ (e.gene.symbol, e.gene.geneid,
                         e.intersect,e.func_class,e.func_type,e.details,e.ref_aa,e.var_aa) for e in evidence if e.gene ]

      geneinfo      = list(unique(geneinfo))

    if cv:
      cvinfo        = cv.score_and_classify(var.chromosome,var.begin,var.end,alleles)
      if cvinfo.exact_vars:
        varid_exact = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)|set(varid_exact))
      if cvinfo.inexact_vars:
        varid_inexact = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.inexact_vars)|set(varid_inexact))

      #cvinfo.common_score)
      #cvinfo.function_score)

    if kaviar and var.end-var.begin==1:
      kinfo = kaviar.get(var.chromosome,var.begin) or []
      kaviar_details += [ stuff.replace(';',',') for chrom,loc,mallele,maf,allele,stuff in kinfo if allele in (allele1,allele2) ]

      if ktext:
        mallele = kinfo[0][2]
        maf     = kinfo[0][3]

        if mallele and mallele in (allele1,allele2):
          kaviar_maf = maf

    #ingroup,outgroup = refvars.get(var.chrom,var.start,var.end,var.var) if refvars else ([],[])

    #ingroup  = ','.join(ingroup)  if  ingroup else ''
    #outgroup = ','.join(outgroup) if outgroup else ''

    print var
    print '    varid_exact=',varid_exact
    print '  varid_inexact=',varid_inexact
    if geneinfo:
      print '      gene info=',geneinfo[0]
      for g in geneinfo[1:]:
        print '                ',g
    print


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('-f', '--format',   metavar='NAME', default='VCF',
                      help='File format (VCF)')
  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--cgfvariants',   metavar='NAME',
                      help='CGFvariant database annotation')
  parser.add_argument('--commonscore', metavar='T', type=float, default=0.05,
                      help='Annotate all variants with common score > T')
  parser.add_argument('--esp',   metavar='NAME',
                        help='UW Exome Sequencing Project (ESP) annotation (optional)')
  parser.add_argument('--kaviar',   metavar='NAME',
                        help='Kaviar annotation (optional)')
  parser.add_argument('--refvars',   metavar='NAME',
                        help='Reference variant list')
  parser.add_argument('--polyphen2',   metavar='NAME',
                        help='Polyphen2 exome annotation (optional)')
  parser.add_argument('--refingroup',   metavar='NAME',
                        help='List of subjects defined to be intra-group reference variants')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()
  format  = options.format.upper()

  if format=='VCF':  # *** VCF, SNPs only ***
    annotate_vcf(options)
  elif format=='MASTERVAR':
    annotate_mastervar(options)
  else:
    raise ValueError('Unknown or Unsupported format specified: %s' % options.format)
