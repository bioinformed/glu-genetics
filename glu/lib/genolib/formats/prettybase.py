# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'PrettyBase genotype format input/output objects'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['PrettybaseWriter', 'save_prettybase', 'load_prettybase']

__genoformats__ = [
  #     LOADER             SAVER              WRITER        PFORMAT     ALIAS      EXTS
  ('load_prettybase', 'save_prettybase', 'PrettybaseWriter', 'trip', 'prettybase', 'pb') ]


import re

from   glu.lib.fileutils         import autofile,namefile,parse_augmented_filename,get_arg, \
                                        get_csv_dialect,trybool,list_reader

from   glu.lib.genolib.streams   import GenotripleStream


def load_prettybase(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load genotype triples from file

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       unique: assume rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @rtype             : GenotripleStream

  >>> from StringIO import StringIO
  >>> data = StringIO('l1 s1 A A\\nl2 s1 G G\\nl1 s2 N N\\nl2 s2 C C\\n')
  >>> triples = load_prettybase(data,'pb')
  >>> for triple in triples:
  ...   print triple
  ('s1', 'l1', ('A', 'A'))
  ('s1', 'l2', ('G', 'G'))
  ('s2', 'l1', (None, None))
  ('s2', 'l2', ('C', 'C'))
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename = parse_augmented_filename(filename,args)

  unique   = trybool(get_arg(args, ['unique'], False))
  order    = get_arg(args, ['order'])
  samples  = get_arg(args, ['samples'])
  loci     = get_arg(args, ['loci'])
  dialect  = get_csv_dialect(args)

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  if samples:
    samples = set(list_reader(samples,**dialect))

  if loci:
    loci = set(list_reader(loci,**dialect))

  re_spaces = re.compile('[\t ,]+')
  gfile = autofile(filename)

  def _load():
    # Micro-optimization
    split        = re_spaces.split
    local_intern = intern
    local_strip  = str.strip
    amap         = {'N':None,'n':None}

    for line_num,line in enumerate(gfile):
      row = split(local_strip(line))
      if not row:
        continue
      elif len(row) != 4:
        raise ValueError('Invalid prettybase row on line %d of %s' % (line_num+1,namefile(filename)))

      locus  = local_intern(local_strip(row[0]))
      sample = local_intern(local_strip(row[1]))
      a1,a2  = row[2],row[3]
      geno   = amap.get(a1,a1),amap.get(a2,a2)

      yield sample,locus,geno

  return GenotripleStream.from_tuples(_load(),genome=genome,phenome=phenome,
                                              samples=samples,loci=loci,
                                              unique=unique,order=order)


class PrettybaseWriter(object):
  '''
  Object to write genotype triple data to a Prettybase format file

  Genotype triple files must be supplied as and are output to whitespace
  delimited ASCII files as a sequence of four items:

    1. Locus name
    2. Sample name
    3. Allele 1, N for missing
    4. Allele 2, N for missing

  All rows output have exactly these four columns and no file header is
  output. Sample and locus names are arbitrary and user-specified strings.

  >>> triples = [('s1','l1',('C','T')), ('s1','l2',(None,None)),
  ...            ('s1','l3',('A','A')), ('s2','l2', ('C','C'))]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> with PrettybaseWriter(o,'pb',None,triples.genome,triples.phenome) as w:
  ...   triples = iter(triples)
  ...   w.writerow(*triples.next())
  ...   w.writerow(*triples.next())
  ...   w.writerows(triples)
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 s1 C T
  l2 s1 N N
  l3 s1 A A
  l2 s2 C C
  '''
  def __init__(self,filename,format,header,genome,phenome,extra_args=None,**kwargs):
    '''
    @param     filename: file name or file object
    @type      filename: str or file object
    '''
    if extra_args is None:
      args = kwargs
    else:
      args = extra_args
      args.update(kwargs)

    filename = parse_augmented_filename(filename,args)

    if extra_args is None and args:
      raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

    self.out = autofile(filename,'w')

  def writerow(self, sample, locus, geno):
    '''
    Write a genotype triple (sample,locus,genotype)

    @param sample: sample identifier
    @type  sample: str
    @param  locus: locus identifier
    @type   locus: str
    @param   geno: genotypes internal representation
    @type    geno: genotype representation
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    out.write( ' '.join( [locus,sample,geno[0] or 'N',geno[1] or 'N'] ) )
    out.write('\n')

  def writerows(self, triples):
    '''
    Write a genotype sequence of triples (sample,locus,genotype)

    @param  triples: sequence of (sample,locus,genotype)
    @type   triples: sequence of (str,str,genotype representation)
    '''
    out = self.out
    if out is None:
      raise IOError('Cannot write to closed writer object')

    write = out.write
    join  = ' '.join

    for sample,locus,geno in triples:
      write( join( [locus,sample,geno[0] or 'N',geno[1] or 'N'] ) )
      write('\n')

  def close(self):
    '''
    Close the writer.

    A closed writer cannot be used for further I/O operations and will
    result in an error if called more than once.
    '''
    if self.out is None:
      raise IOError('Writer object already closed')
    self.out = None

  def __enter__(self):
    '''
    Context enter function
    '''
    return self

  def __exit__(self, *exc_info):
    '''
    Context exit function that closes the writer upon exit
    '''
    self.close()


def save_prettybase(filename,genos,format,extra_args=None,**kwargs):
  '''
  Write the genotype triple data to file.

  @param     filename: file name or file object
  @type      filename: str or file object
  @param      triples: genotype triple data
  @type       triples: sequence

  >>> triples = [ ('s1', 'l1',  ('C','T')),
  ...             ('s1', 'l2', (None,None)),
  ...             ('s1', 'l3',  ('A','A')) ]
  >>> triples = GenotripleStream.from_tuples(triples)
  >>> from cStringIO import StringIO
  >>> o = StringIO()
  >>> save_prettybase(o,triples,'pb')
  >>> print o.getvalue() # doctest: +NORMALIZE_WHITESPACE
  l1 s1 C T
  l2 s1 N N
  l3 s1 A A
  '''
  with PrettybaseWriter(filename,format,None,genos.genome,genos.phenome,extra_args=extra_args,**kwargs) as w:
    w.writerows(genos.as_genotriples())


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
