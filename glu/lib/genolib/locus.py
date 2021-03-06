# -*- coding: utf-8 -*-

__abstract__  = 'locus model input and output'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   collections               import defaultdict

from   glu.lib.fileutils         import namefile,table_reader,get_arg,parse_augmented_filename

from   glu.lib.genolib.genoarray import build_model


class Nothing(object): pass


STRAND_UNKNOWN,STRAND_FORWARD,STRAND_REVERSE = None,'+','-'

STRAND_MAP = {'+':STRAND_FORWARD, 'F':STRAND_FORWARD, 'FORWARD':STRAND_FORWARD,
              '-':STRAND_REVERSE, 'R':STRAND_REVERSE, 'REVERSE':STRAND_REVERSE,
               '':STRAND_UNKNOWN, 'U':STRAND_UNKNOWN, 'UNKNOWN':STRAND_UNKNOWN,
              '?':STRAND_UNKNOWN}


# FIXME: Try various cases
def _get_header(filename,hmap,header,required=False):
  index = hmap.get(header,[])

  if not index:
    if required:
      raise ValueError('Cannot find %s header in locus file %s' % (header,namefile(filename)))
    return None

  elif len(index) != 1:
    raise ValueError('Ambiguous %s header in locus file %s' % (header,namefile(filename)))

  return index[0]


class Locus(object):
  '''
  Locus metadata
  '''
  __slots__ = ('name','model','chromosome','location','strand')

  def __init__(self, name, model=None, chromosome=None, location=None, strand=STRAND_UNKNOWN):
    self.name       = name
    self.model      = model
    self.chromosome = chromosome
    self.location   = location
    self.strand     = strand


class Genome(object):
  '''
  Genome metadata storage, default genotype model factory and collection of
  locus descriptors

  FIXME: docstring
  '''

  def __init__(self, default_model=None, default_fixed=False, max_alleles=None, alleles=None, loci=None):
    '''
    FIXME: docstring
    '''
    if not alleles and not max_alleles:
      max_alleles = 2

    self.loci          = dict( (locus.name,locus) for locus in loci or [])
    self.max_alleles   = max_alleles or 0
    self.alleles       = alleles or []
    self.default_model = default_model
    self.default_fixed = default_fixed

  def set_locus(self, name, model=Nothing, chromosome=Nothing,
                            location=Nothing, strand=Nothing):
    '''
    FIXME: docstring
    '''
    locus = self.loci.get(name)
    if locus is None:
      locus = self.loci[name] = Locus(name)

    if model is not Nothing:
      locus.model = model
    if chromosome is not Nothing:
      locus.chromosome = chromosome
    if location is not Nothing:
      locus.location = location
    if strand is not Nothing:
      locus.strand = strand

  def merge_locus(self, name, model=None, chromosome=None,
                              location=None, strand=STRAND_UNKNOWN,
                              warn=False):
    '''
    FIXME: docstring
    '''
    locus = self.loci.get(name)

    # Fastpath: No merge needed
    if locus is None:
      self.loci[name] = Locus(name,model,chromosome,location,strand)
      return

    # Slowpath: Must check for incompatibilities
    if model is not None:
      if locus.model is None:
        locus.model = model
      elif locus.model is not model:
        msg = 'Locus record %s merging incompatible models' % name
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if chromosome is not None:
      if locus.chromosome is None:
        locus.chromosome = chromosome
      elif locus.chromosome != chromosome:
        msg = 'Locus record %s merging incompatible chromosomes (%s != %s)' % (name,locus.chromosome,chromosome)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if location is not None:
      if locus.location is None:
        locus.location = location
      elif locus.location != location:
        msg = 'Locus record %s merging incompatible locations (%s != %s)' % (name,locus.location,location)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

    if strand not in (STRAND_UNKNOWN,Nothing):
      if locus.strand is STRAND_UNKNOWN:
        locus.strand = strand
      elif locus.strand != strand:
        msg = 'Locus record %s merging incompatible strands (%s != %s)' % (name,locus.strand,strand)
        if warn:
          sys.stderr.write('[WARNING] %s\n' % msg)
        else:
          raise ValueError(msg)

  def get_locus(self, name):
    '''
    FIXME: docstring
    '''
    locus = self.loci.get(name)

    if locus is None:
      locus = self.loci[name] = Locus(name)

    return locus

  def get_locus_model(self, name):
    '''
    FIXME: docstring
    '''
    locus = self.get_locus(name)

    if locus.model is None:
      if self.default_model is not None:
        locus.model = self.default_model
      else:
        locus.model = build_model(self.alleles, max_alleles=self.max_alleles)

    return locus

  def get_model(self, name):
    '''
    FIXME: docstring
    '''
    return self.get_locus_model(name).model


def load_locus_records(filename,extra_args=None,**kwargs):
  '''
  Load a locus file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      alleles: alleles
  @param      alleles: sequence of objects
  @param  max_alleles: maximum number of alleles allowed. Default is 2
  @type   max_alleles: int or None
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @return            : default_max_alleles, default_alleles, sequence of
                       locus models of the form: locus name, max_alleles, alleles

  >>> from StringIO import StringIO
  >>> l = """LOCUS\\tMAX_ALLELES\\tALLELES\\tCHROMOSOME\\tLOCATION
  ... l1\\t2\\tA,C\\t1\\t0
  ... l2
  ... l3\\t\\tA,G\\t1
  ... l4\\t45\\t\\tX\\t1234
  ... l5\\t45\\tA,B,C,D,E,F,G,H,I,J,K\\tM\\t5541
  ... """
  >>> max_alleles,alleles,loci = load_locus_records(StringIO(l))
  >>> max_alleles
  2
  >>> alleles
  []
  >>> for lname,max_alleles,alleles,chromosome,location,strand in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 2 [] None None
  l3 2 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> max_alleles,alleles,loci = load_locus_records(StringIO(l),max_alleles=4,alleles=['A','C','G','T'])
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles,chromosome,location,strand in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 4 ['A', 'C', 'G', 'T'] None None
  l3 4 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> import tempfile
  >>> f = tempfile.NamedTemporaryFile()
  >>> f.write(l)
  >>> f.flush()
  >>> max_alleles,alleles,loci = load_locus_records(f.name + ':max_alleles=4:alleles=A,C,G,T')
  >>> for lname,max_alleles,alleles,chromosome,location,strand in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  l1 2 ['A', 'C'] 1 0
  l2 4 ['A', 'C', 'G', 'T'] None None
  l3 4 ['A', 'G'] 1 None
  l4 45 [] X 1234
  l5 45 ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] M 5541

  >>> max_alleles,alleles,loci = load_locus_records(':max_alleles=4:alleles=A,C,G,T')
  >>> max_alleles
  4
  >>> alleles
  ['A', 'C', 'G', 'T']
  >>> for lname,max_alleles,alleles,chromosome,location in loci:
  ...   print lname,max_alleles,alleles,chromosome,location
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  if isinstance(filename,basestring) and filename.startswith(':'):
    filename = parse_augmented_filename(filename,args)
    assert not filename
  elif filename:
    rows = table_reader(filename,want_header=True,extra_args=args)
    header = rows.next()

  default_alleles     = get_arg(args, ['alleles'], None)
  default_max_alleles = int(get_arg(args, ['max_alleles'], 2)) or None

  if default_alleles is None:
    default_alleles = []
  elif isinstance(default_alleles,basestring):
    default_alleles = [ a.strip() for a in default_alleles.split(',') ]

  if default_max_alleles is not None and default_alleles is not None:
    default_max_alleles = max(default_max_alleles,len(default_alleles))

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if not filename:
    return default_max_alleles,default_alleles,[]

  hmap = defaultdict(list)
  for i,h in enumerate(header):
    hmap[h].append(i)

  # FIXME: Add aliases -- POSITION,MARKERNAME,NAME,CHR
  lname_index       = _get_header(filename,hmap,'LOCUS',required=True)
  max_alleles_index = _get_header(filename,hmap,'MAX_ALLELES')
  alleles_index     = _get_header(filename,hmap,'ALLELES')
  chromosome_index  = _get_header(filename,hmap,'CHROMOSOME')
  location_index    = _get_header(filename,hmap,'LOCATION')
  strand_index      = _get_header(filename,hmap,'STRAND')

  def _loci():
    for i,row in enumerate(rows):
      if not row:
        continue

      n = len(row)

      if lname_index >= n:
        lname = ''
      else:
        lname = intern(row[lname_index].strip())

      if not lname:
        raise ValueError('Invalid locus file record %d in %s (blank locus name)' % (i+1,namefile(filename)))

      informative = False

      max_alleles = default_max_alleles
      if max_alleles_index is not None and max_alleles_index<n:
        m = row[max_alleles_index].strip()
        if m:
          max_alleles = int(m)
          informative = max_alleles>0

      alleles = None
      if alleles_index is not None and alleles_index<n:
        alleles = ( a.strip() for a in row[alleles_index].split(',') )
        alleles = [ intern(a) for a in alleles if a ]

      if not alleles and informative:
        alleles = []
      elif not alleles:
        alleles = default_alleles

      max_alleles = max(max_alleles or 0,len(alleles)) or None

      chromosome = None
      if chromosome_index is not None and chromosome_index<n:
        chromosome = intern(row[chromosome_index].strip()) or None

      location = None
      if location_index is not None and location_index<n:
        loc = row[location_index].strip() or None
        if loc:
          location = int(loc)

      strand = None
      if strand_index is not None and strand_index<n:
        strand = intern(row[strand_index].strip()) or None

      yield lname,max_alleles,alleles,chromosome,location,strand

  return default_max_alleles,default_alleles,_loci()


def populate_genome(genome,loci):
  '''
  Return the default model and a sequence of Locus objects from an augmented
  locus description file

  FIXME: Add docstring and doctests
  '''
  # FIXME: Model merge must be more careful about stepping on existing
  #        models, when they are actually compatible.  The current interface
  #        is required identical model objects, so problems will result if
  #        populate genome is called on a non-empty instance.

  for lname,max_alleles,alleles,chromosome,location,strand in loci:
    model = build_model(alleles,max_alleles=max_alleles)
    genome.merge_locus(lname, model, chromosome, location, strand)


def load_genome(filename,**kwargs):
  '''
  Return the default model and a sequence of Locus objects from an augmented
  locus description file

  FIXME: Add docstring and doctests
  '''
  default_max_alleles,default_alleles,loci = load_locus_records(filename,**kwargs)

  if default_alleles:
    default_model = build_model(default_alleles,max_alleles=default_max_alleles)
  else:
    default_model = None

  genome = Genome(default_model=default_model, default_fixed=bool(default_alleles),
                  max_alleles=default_max_alleles, alleles=default_alleles)

  populate_genome(genome,loci)

  return genome


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
