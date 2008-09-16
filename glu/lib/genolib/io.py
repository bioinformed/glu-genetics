# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'GLU genotype data input/output objects'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


from   glu.lib.utils             import is_str
from   glu.lib.fileutils         import namefile, guess_format, parse_augmented_filename, get_arg

from   glu.lib.genolib.streams   import GenotripleStream, GenomatrixStream
from   glu.lib.genolib.transform import GenoTransform
from   glu.lib.genolib.merge     import get_genomerger
from   glu.lib.genolib.locus     import load_genome, Genome
from   glu.lib.genolib.phenos    import load_phenome, Phenome
from   glu.lib.genolib.reprs     import get_genorepr

# FIXME: Format support should ultimately be pluggable with a registration protocol
from   glu.lib.genolib.formats   import *


INPUT_FORMATS  = ['ldat','sdat','tdat','trip','genotriple',
                  'prettybase','pb',
                  'lbat','sbat','tbat',
                  'ped','tped','bed',
                  'hapmap','mach','merlin', 'wtccc-raw',
                  'eigensoft','smartpca']

OUTPUT_FORMATS = ['ldat','sdat','tdat','trip','genotriple',
                  'prettybase','pb',
                  'lbat','sbat','tbat',
                  'ped','tped','bed',
                  'mach','merlin','structure','phase','wtccc',
                  'eigensoft','smartpca']


def guess_informat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  return guess_format(filename, INPUT_FORMATS)


def guess_informat_list(filenames):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  formats = set( guess_informat(f) for f in filenames )
  formats.discard(None)
  if len(formats) == 1:
    return formats.pop()
  return None


def guess_outformat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  return guess_format(filename, OUTPUT_FORMATS)


def geno_options(group,input=False,output=False,merge=False,filter=False,transform=False):
  if input:
    group.add_option('-f', '--informat', dest='informat', metavar='NAME',
                      help='Input genotype format')
    group.add_option('-g', '--ingenorepr', dest='ingenorepr', metavar='REP',
                      help='Input genotype representation')

  if output:
    group.add_option('-F','--outformat',  dest='outformat', metavar='NAME',
                      help='Output genotype format')
    group.add_option('-G', '--outgenorepr', dest='outgenorepr', metavar='REP',
                      help='Output genotype representation')

  if input or output:
    group.add_option('-l', '--loci', dest='loci', metavar='FILE',
                      help='Locus description file and options')
    group.add_option('-p', '--pedigree', dest='pedigree', metavar='FILE',
                      help='Pedigree description file and options')

  if merge:
    group.add_option('--merge', dest='merge', metavar='METHOD:T', default='unanimous',
                      help='Genotype merge algorithm and optional consensus threshold used to form a '
                           'consensus genotypes. Values=unique,unanimous,vote,ordered.  Value may be '
                           'optionally followed by a colon and a threshold.  Default=unanimous')
    group.add_option('--samplemerge', dest='samplemerge', metavar='FILE',
                    help='Sample concordance statistics output to FILE (optional)')
    group.add_option('--locusmerge',  dest='locusmerge',  metavar='FILE',
                    help='Locus concordance statistics output to FILE (optional)')

  if filter:
    group.add_option('--filtermissing', action='store_true', dest='filtermissing',
                      help='Filters out the samples or loci with missing genotypes')

    group.add_option('--includesamples', dest='includesamples', metavar='FILE',
                      help='List of samples to include')
    group.add_option('--includeloci', dest='includeloci', metavar='FILE',
                      help='List of loci to include')

    group.add_option('--excludesamples', dest='excludesamples', metavar='FILE',
                      help='List of samples to exclude')
    group.add_option('--excludeloci', dest='excludeloci', metavar='FILE',
                      help='List of loci to exclude')

  if transform:
    group.add_option('--renamesamples', dest='renamesamples', metavar='FILE',
                      help='Rename samples from a file containing rows of original name, tab, new name')
    group.add_option('--renameloci', dest='renameloci', metavar='FILE',
                      help='Rename loci from a file containing rows of original name, tab, new name')
    group.add_option('--renamealleles', dest='renamealleles', metavar='FILE',
                      help='Rename alleles based on file of locus name, tab, old alleles (comma separated), '
                           'tab, new alleles (comma separated)')
    group.add_option('--ordersamples', dest='ordersamples', metavar='FILE',
                      help='Order samples based on the order of names in FILE')
    group.add_option('--orderloci', dest='orderloci', metavar='FILE',
                      help='Order loci based on the order of names in FILE')


def load_genostream(filename, transform=None, extra_args=None, **kwargs):
  '''
  Load genomatrix file depending on matrix format and return a GenotripleMatrix object

  @param filename: a file name or file object
  @type  filename: str or file object
  @param   format: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                   lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type    format: str
  @param  genorepr: string or representation object for text genotypes. Default is None
  @type   genorepr: str, UnphasedMarkerRepresentation or similar object
  @param   unique: flag indicating if repeated row or column elements do not exist. Default is None
  @type    unique: bool
  @para    genome: Genome metadata or filename
  @type    genome: str or Genome instance
  @para   phenome: Pedigree and phenotype metadata or filename
  @type   phenome: str or Phenome instance
  @return        : loaded genomatrix stream
  @rtype         : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\tCC\\tCT\\tTT\\n")
  >>> ldat = load_genostream(data,format='ldat',genorepr='snp')
  >>> ldat.columns
  ('s1', 's2', 's3')
  >>> for row in ldat:
  ...   print row
  ('l1', [('A', 'A'), ('A', 'G'), ('G', 'G')])
  ('l2', [('C', 'C'), ('C', 'T'), ('T', 'T')])
  >>> ldat.loci
  >>> ldat.unique
  True

  >>> from StringIO import StringIO
  >>> data = StringIO('s1\\tl1\\tAA\\ns1\\tl2\\tGG\\ns2\\tl1\\tAG\\ns2\\tl2\\tCC\\n')
  >>> triples = load_genostream(data,format='tdat',genorepr='snp')
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', ('A', 'G'))
  ('s2', 'l2', ('C', 'C'))
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  hyphen  = get_arg(args, ['hyphen'])
  format  = get_arg(args, ['format'])
  genome  = get_arg(args, ['genome'])
  phenome = get_arg(args, ['phenome'])

  if filename == '-':
    if hyphen is None:
      raise ValueError("loading genotypes from '-' is not supported")
    if is_str(hyphen):
      raise ValueError('a file object must be supplied for hyphen redirection')
    filename = hyphen

  if format is None:
    format = guess_informat(filename)

  # FIXME: Check genome is None behavior
  if genome is None:
    genome = Genome()
  elif is_str(genome):
    genome = load_genome(genome,extra_args=args)

  if phenome is None:
    phenome = Phenome()
  elif is_str(phenome):
    phenome = load_phenome(phenome,extra_args=args)

  kwargs = dict(genome=genome,phenome=phenome)

  if format == 'hapmap':
    genos = load_hapmap(filename,extra_args=args,**kwargs)
  elif format == 'ldat':
    genos = load_genomatrix_text(filename,format,extra_args=args,**kwargs)
  elif format == 'sdat':
    genos = load_genomatrix_text(filename,format,extra_args=args,**kwargs)
  elif format == 'lbat':
    genos = load_genomatrix_binary(filename,'ldat',extra_args=args,**kwargs)
  elif format == 'sbat':
    genos = load_genomatrix_binary(filename,'sdat',extra_args=args,**kwargs)
  elif format in ('tdat','trip','genotriple'):
    genos = load_genotriples_text(filename,extra_args=args,**kwargs)
  elif format in ('pb','prettybase'):
    genos = load_prettybase(filename,extra_args=args,**kwargs)
  elif format=='tbat':
    genos = load_genotriples_binary(filename,extra_args=args,**kwargs)
  elif format in ('plink_ped','ped'):
    genos = load_plink_ped(filename,extra_args=args,**kwargs)
  elif format in ('plink_tped','tped'):
    genos = load_plink_tped(filename,extra_args=args,**kwargs)
  elif format in ('plink_bed','bed'):
    genos = load_plink_bed(filename,extra_args=args,**kwargs)
  elif format in ('merlin','mach'):
    genos = load_merlin(filename,extra_args=args,**kwargs)
  elif format in ('eigensoft','smartpca'):
    genos = load_eigensoft_smartpca(filename,extra_args=args,**kwargs)
  elif format == 'wtccc-raw':
    genos = load_wtccc_raw(filename,extra_args=args,**kwargs)
  elif not format:
    raise ValueError("Input file format for '%s' must be specified" % namefile(filename))
  else:
    raise NotImplementedError("File format '%s' is not supported" % format)

  # Apply any requested transformations
  # FIXME: The length checking code is there because transformed may not be
  #        smart enough to handle null transformations as a matter of
  #        course.  This can be relaxed once that path is rechecked.
  transform_args = len(args)
  more_transform = GenoTransform.from_kwargs(args)

  if len(args) != transform_args:
    if transform:
      transform = GenoTransform.from_object(transform).merge(more_transform)
    else:
      transform = more_transform

  if transform:
    genos = genos.transformed(transform)

  # Kludge until extra_args is smarter about default values.
  #   This is needed in the short term for genorepr=None for formats that do
  #   not use genoreprs.
  for k,v in args.items():
    if v is None:
      del args[k]

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  return genos


def save_genostream(filename, genos, extra_args=None, **kwargs):
  '''
  Write genotype data to file

  @param  filename: a file name or file object
  @type   filename: str or file object
  @param     genos: genomatrix/genotriple stream
  @type      genos: sequence
  @param    format: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                    lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type     format: str
  @param  genorepr: string or representation object for text genotypes. Default is None
  @type   genorepr: str, UnphasedMarkerRepresentation or similar object
  @param mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type  mergefunc: callable
  @param  compress: flag indicating if a compressed format is desired. Default is True
  @type   compress: bool
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  hyphen    = get_arg(args, ['hyphen'])
  format    = get_arg(args, ['format'])

  if filename == '-':
    if hyphen is None:
      raise ValueError("saving genotypes to '-' is not supported")
    if is_str(hyphen):
      raise ValueError('a file object must be supplied for hyphen redirection')
    filename = hyphen

  if format is None:
    format = guess_outformat(filename)

  if format == 'ldat':
    save_genomatrix_text(filename, genos, format='ldat', extra_args=args)
  elif format == 'sdat':
    save_genomatrix_text(filename, genos, format='sdat', extra_args=args)
  elif format in ('tdat','trip','genotriple'):
    save_genotriples_text(filename, genos, extra_args=args)
  elif format in ('pb','prettybase'):
    genos = save_prettybase(filename, genos, extra_args=args)
  elif format == 'lbat':
    save_genomatrix_binary(filename, genos, format='ldat', extra_args=args)
  elif format == 'sbat':
    save_genomatrix_binary(filename, genos, format='sdat', extra_args=args)
  elif format == 'tbat':
    save_genotriples_binary(filename, genos, extra_args=args)
  elif format in ('plink_ped','ped'):
    save_plink_ped(filename, genos, extra_args=args)
  elif format in ('plink_tped','tped'):
    save_plink_tped(filename, genos, extra_args=args)
  elif format in ('plink_bed','bed'):
    save_plink_bed(filename, genos, extra_args=args)
  elif format in ('merlin','mach'):
    save_merlin(filename, genos, extra_args=args)
  elif format == 'structure':
    save_structure(filename, genos, extra_args=args)
  elif format == 'phase':
    save_phase(filename, genos, extra_args=args)
  elif format == 'wtccc':
    save_wtccc(filename, genos, extra_args=args)
  elif format in ('eigensoft','smartpca'):
    genos = save_eigensoft_smartpca(filename, genos, extra_args=args)
  elif not format:
    raise ValueError("Output file format for '%s' must be specified" % namefile(filename))
  else:
    raise NotImplementedError("File format '%s' is not supported" % format)


def transform_files(infiles,informat,ingenorepr,
                    outfile,outformat,outgenorepr,
                    transform=None,genome=None,phenome=None,
                    mergefunc=None,
                    inhyphen=None,outhyphen=None):
  '''
  A driver for transforming multiple genodata files into different formats
  (ldat, sdat, trip, or genotriples), representations (...) and, depending
  on the presence and attributes of the transform object, performing
  operations on samples and loci such as exclude, include, and rename.  The
  results are then saved to the specified output file.

  @param     infiles: list of input file names or file objects
  @type      infiles: str or file objects
  @param    informat: format of input file: hapmap,ldat,sdat,trip,genotriple,prettybase,pb,
                      lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type     informat: str
  @param  ingenorepr: internal genotype representation for the input
  @type   ingenorepr: UnphasedMarkerRepresentation or similar object
  @param    outfiles: output file name or file object
  @type     outfiles: str or file object
  @param   outformat: format of output file: ldat,sdat,trip,genotriple,prettybase,pb,
                      lbat,sbat,tbat. Default is None, which attempts to autodetect the file type.
  @type    outformat: str
  @param outgenorepr: internal genotype representation for the output
  @type  outgenorepr: UnphasedMarkerRepresentation or similar object
  @param   transform: transformation object (optional)
  @type    transform: GenoTransform object
  @para       genome: Genome metadata or filename
  @type       genome: str or Genome instance
  @para      phenome: Pedigree and phenotype metadata or filename
  @type      phenome: str or Phenome instance
  @param   mergefunc: function to merge multiple genotypes into a consensus genotype. Default is None
  @type    mergefunc: callable

  >>> from StringIO import StringIO
  >>> data = StringIO("ldat\\ts1\\ts2\\ts3\\nl1\\tAA\\tAG\\tGG\\nl2\\t\\tCT\\tTT\\n")
  >>> out  = StringIO()
  >>> transform_files([data],'ldat','snp',out,'tdat','marker')
  >>> print out.getvalue() # doctest: +NORMALIZE_WHITESPACE
  s1  l1      A/A
  s2  l1      A/G
  s3  l1      G/G
  s1  l2
  s2  l2      C/T
  s3  l2      T/T
  '''
  if informat is None:
    informat = guess_informat_list(infiles)

  if is_str(ingenorepr):
    ingenorepr = get_genorepr(ingenorepr)

  if is_str(outgenorepr):
    outgenorepr = get_genorepr(outgenorepr)

  if not outgenorepr:
    outgenorepr = ingenorepr

  if genome is None:
    genome = Genome()
  elif is_str(genome):
    genome = load_genome(genome)

  if phenome is None:
    phenome = Phenome()
  elif is_str(phenome):
    phenome = load_phenome(phenome)

  if is_str(mergefunc):
    mergefunc = get_genomerger(mergefunc)

  genos = [ load_genostream(f,format=informat,genorepr=ingenorepr,
                              genome=genome,phenome=phenome,hyphen=inhyphen)
                            .transformed(transform) for f in infiles ]
  n = len(genos)

  if outformat is None:
    outformat = guess_outformat(outfile)

  # Guess output format based on input format if it is unique
  if outformat is None:
    outformat = informat

  # FIXME: Refactor preferred output classes
  if outformat in ('ldat','lbat','plink_tped','tped','plink_bed','bed','wtccc','eigensoft','smartpca'):
    genos = GenomatrixStream.from_streams(genos,'ldat',mergefunc=mergefunc)
  elif outformat in ('sdat','sbat','plink_ped','ped','plink_bed_ind','merlin','mach','structure','phase'):
    genos = GenomatrixStream.from_streams(genos,'sdat',mergefunc=mergefunc)
  elif outformat in ('tdat','trip','genotriple','pb','prettybase','tbat'):
    genos = GenotripleStream.from_streams(genos,mergefunc=mergefunc)
  elif not outformat:
    raise ValueError("Output file format for '%s' must be specified" % namefile(outfile))
  else:
    raise NotImplementedError("File format '%s' is not supported" % outformat)

  # Order again after merging, if necessary
  if n>1 and (transform.loci.order or transform.samples.order):
    genos = genos.transformed(order_loci=transform.loci.order,
                              order_samples=transform.samples.order)

  save_genostream(outfile,genos,format=outformat,genorepr=outgenorepr,hyphen=outhyphen)


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
