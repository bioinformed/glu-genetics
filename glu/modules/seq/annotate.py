# -*- coding: utf-8 -*-

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                     import itemgetter

from   glu.lib.fileutils            import table_writer, table_reader, list_reader, autofile, hyphen
from   glu.lib.progressbar          import progress_loop

from   glu.lib.seqlib.vcf           import VCFReader
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


def update_vcf_annotation(v, vs, cv, kaviar, refvars, options):
  if vs:
    # FIXME: Order genes and evidence consistently
    evidence   = list(vs.classify(v.chrom, v.start, v.end, v.var[0], nsonly=False))
    #v.names   = sorted(set(str(v) for e in evidence for v in e.varid_exact)|set(v.names))
    cytoband   = sorted(set(e.cytoband    for e in evidence if e.cytoband))
    genes      = sorted(set(e.gene.symbol for e in evidence if e.gene and e.gene.symbol))
    geneids    = sorted(set(e.gene.geneid for e in evidence if e.gene and e.gene.geneid))
    location   = sorted(set(e.intersect   for e in evidence if e.intersect))
    function   = sorted(set(e.func_type or e.func_class  for e in evidence if e.func_type or e.func_class))
    nsevidence = [ e for e in evidence if 'NON-SYNONYMOUS' in e.func_class ]
    nsinfo     = [ '%s:%s:%s->%s' % (e.gene.symbol,e.details,e.ref_aa,e.var_aa)
                   for e in nsevidence ]
    if not genes:
      v.filter.append('NonGenic')
    elif not nsevidence:
      v.filter.append('Synonymous')

    new_info = ['CYTOBAND=%s'              % (','.join(cytoband)),
                'GENE_NAME=%s'             % (','.join(genes   )),
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
      else:
        new_info.append('KAVIAR_MAF=')

      new_info.append('KAVIAR_NAMES=%s' % ','.join(ktext))

    else:
      new_info.append('KAVIAR_MAF=0')
      new_info.append('KAVIAR_NAMES=')

  ingroup,outgroup = refvars.get(v.chrom,v.start,v.end,v.var) if refvars else ([],[])

  new_info.append('REFVARS_INGROUP_COUNT=%d'  % len(ingroup))
  new_info.append('REFVARS_OUTGROUP_COUNT=%d' % len(outgroup))

  ingroup  = ','.join(ingroup)  if  ingroup else ''
  outgroup = ','.join(outgroup) if outgroup else ''

  new_info.append('REFVARS_INGROUP_NAMES=%s'  % ingroup)
  new_info.append('REFVARS_OUTGROUP_NAMES=%s' % outgroup)

  if outgroup:
    v.filter.append('RefVar')

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

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars    = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  metadata = vcf.metadata
  metadata['FILTER'].append(['##FILTER=<PartiallyCalled,Description="Variant is not called for one or more samples">'])
  metadata['FILTER'].append(['##FILTER=<NonGenic,Description="Variant not in or near a gene">'])
  metadata['FILTER'].append(['##FILTER=<Synonymous,Description="Variant does not alter an amino-acid">'])

  if cv:
    metadata['FILTER'].append(['##FILTER=<Common,Description="Variant is likely common with common score>%f">' % options.commonscore])
    metadata['FILTER'].append(['##FILTER=<1000G,Description="Variant was reported by 1000 Genomes project">'])

  if kaviar:
    metadata['FILTER'].append(['##FILTER=<KaviarCommon,Description="Variant appears in the Kaviar database and appears to be common with MAF>%f">' % options.commonscore])
    metadata['FILTER'].append(['##FILTER=<Kaviar,Description="Variant appears in the Kaviar database">'])

  if refvars:
    metadata['FILTER'].append(['##FILTER=<RefVar,Description="Variant appears in the local reference variant list">'])

  #metadata['FILTER'].append(['##FILTER=<NotDominant,Description="Variant does not fit dominant heritibility model">'])
  #metadata['FILTER'].append(['##FILTER=<NotRecessive,Description="Variant does not fit recessive heritibility model">'])

  metadata['INFO'].append(['##INFO=<ID=CYTOBAND,Number=.,Type=String,Description="Name of cytoband(s) containing variant">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_NAME,Number=.,Type=String,Description="Name of gene(s) containing variant">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_ID,Number=.,Type=String,Description="Entrez/LocusLink gene identifiers of genes containing variant">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_LOCATION,Number=.,Type=String,Description="Location of variant in gene(s)">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_FUNCTION,Number=.,Type=String,Description="Functional classification of variant for each gene and transcript">'])
  metadata['INFO'].append(['##INFO=<ID=GENE_FUNCTION_DETAILS,Number=.,Type=String,Description="Functional details of variant for each gene and transcript">'])

  if cv:
    metadata['INFO'].append(['##INFO=<ID=COMMON_SCORE,Number=1,Type=Float,Description="Common score: maximum allele frequency in any population for rarest allele">'])
    metadata['INFO'].append(['##INFO=<ID=FUNCTION_SCORE,Number=1,Type=Int,Description="Function score: reported as functional variant in OMIM, dbSNP, or COSMIC">'])
    metadata['INFO'].append(['##INFO=<ID=INEXACT_VARIANTS,Number=.,Type=String,Description="Inexact variant matche">'])

  if kaviar:
    metadata['INFO'].append(['##INFO=<ID=KAVIAR_MAF,Number=1,Type=Float,Description="Minor allele frequency accourding to Kaviar database">'])
    metadata['INFO'].append(['##INFO=<ID=KAVIAR_NAMES,Number=.,Type=String,Description="Samples or datasets from Kaviar in which variant was found">'])

  if refvars:
    metadata['INFO'].append(['##INFO=<ID=REFVAR_INGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in intra-group samples">'])
    metadata['INFO'].append(['##INFO=<ID=REFVAR_OUTGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in extra-group samples">'])
    metadata['INFO'].append(['##INFO=<ID=REFVAR_INGROUP_NAMES,Number=.,Type=String,Description="Intra-group samples in which Variant is present">'])
    metadata['INFO'].append(['##INFO=<ID=REFVAR_OUTGROUP_NAMES,Number=.,Type=String,Description="Extra-group samples in which Variant is present">'])

  for meta in vcf.metadata_order:
    for m in metadata[meta]:
      out.write('\t'.join(m))
      out.write('\n')

  out.write('#%s\n' % ('\t'.join(vcf.header)))

  for v in vcf:
    update_vcf_annotation(v, vs, cv, kaviar, refvars, options)

    # FORMAT: chrom start end names ref var filter info format genos
    row = [ v.chrom, str(v.end), ','.join(v.names) or '.', v.ref, ','.join(v.var), v.qual,
                                 ';'.join(sorted(v.filter)) or '.',
                                 ';'.join(v.info), v.format ] + [ ':'.join(g) for g in v.genos ]

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
  parser.add_argument('--refvars',   metavar='NAME',
                        help='Reference variant list')
  parser.add_argument('--refingroup',   metavar='NAME',
                        help='List of subjects defined to be intra-group reference variants')
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
