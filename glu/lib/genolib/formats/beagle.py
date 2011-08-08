# -*- coding: utf-8 -*-

from __future__ import with_statement

__abstract__  = 'Beagle genotype format output object'
__copyright__ = 'Copyright (c) 2007-2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['load_beagle_dosage']

  #LOADER        SAVER          WRITER       PFORMAT  ALIAS    EXTS
__genoformats__ = []


from   glu.lib.fileutils         import autofile,parse_augmented_filename,get_arg,guess_related_file,namefile


def load_beagle_samples(filename):
  sfile = autofile(filename)
  header = sfile.next().split()

  if header[:3] != ['ID_1','ID_2','missing']:
    raise ValueError('Invalid header found for Beagle sample file %s' % namefile(filename))

  header2 = sfile.next().split()

  if header2[:3] != ['0','0','0']:
    raise ValueError('Invalid header2 found for Beagle sample file %s' % namefile(filename))

  for line_num,line in enumerate(sfile):
    fields = line.split()

    sample = ''

    if not fields:
      pass
    elif len(fields)==1:
      sample = fields[0]
    elif fields[0] and fields[1] and fields[0]!=fields[1]:
      sample = '%s:%s' % (fields[0],fields[1])
    else:
      sample = fields[0] or fields[1]

    if not sample:
      raise ValueError('Invalid Beagle sample data on line %d of %s' % (line_num+1,namefile(filename)))

    yield sample


def load_beagle_dosage(filename,format,genome=None,phenome=None,extra_args=None,**kwargs):
  '''
  Load a Beagle file and return a NumPy array of dosage values

  See http://www.stats.ox.ac.uk/%7Emarchini/software/gwas/file_format.html

  @param     filename: file name or file object
  @type      filename: str or file object
  @param       format: text string expected in the first header field to
                       indicate data format, if specified
  @type        format: string
  @param     genorepr: function to convert list genotype strings to desired
                       internal representation
  @type      genorepr: unary function
  @param       genome: genome descriptor
  @type        genome: Genome instance
  @param       unique: rows and columns are uniquely labeled (default is True)
  @type        unique: bool
  @param   extra_args: optional dictionary to store extraneous arguments, instead of
                       raising an error.
  @type    extra_args: dict
  @rtype             : GenomatrixStream

  >>> from StringIO import StringIO
  >>> data = StringIO(
  ... 'marker alleleA alleleB\\n'
  ... 'l1 A G 1 0 0 0 1 0 0 0 1\\n'
  ... 'l2 C G 0 0 0 0 1 0 0 0 0\\n'
  ... 'l3 C G 0 0 0 0 1 0 0 0 0\\n'
  ... 'l4 C T 0 .5 .5 1 0 0 0 1 0\\n')
  >>> samples = StringIO('ID_1 ID_2 missing\\n0 0 0 0 0 0\\nf1 s1 0\\nf1 s2 0\\nf1 s3 0\\n')
  >>> samples,genos = load_beagle_dosage(data,'beagle',samples=samples)
  >>> samples
  ('f1:s1', 'f1:s2', 'f1:s3')
  >>> for row in genos:
  ...   print row
  ('l1', None, None, 'A', 'G', array([ 0. ,  0.5,  1. ]))
  ('l2', None, None, 'C', 'G', array([ nan,  0.5,  nan]))
  ('l3', None, None, 'C', 'G', array([ nan,  0.5,  nan]))
  ('l4', None, None, 'C', 'T', array([ 0.75,  0.  ,  0.5 ]))
  '''
  if extra_args is None:
    args = kwargs
  else:
    args = extra_args
    args.update(kwargs)

  filename  = parse_augmented_filename(filename,args)

  format    = get_arg(args, ['format'])
  samples   = get_arg(args, ['samples','s']) or guess_related_file(filename,['lst'])

  if samples is None:
    raise ValueError('Sample file must be specified when loading Beagle files')

  if extra_args is None and args:
    raise ValueError('Unexpected filename arguments: %s' % ','.join(sorted(args)))

  samples = tuple(load_beagle_samples(samples))
  gfile   = autofile(filename)

  header  = next(gfile).rstrip('\n').split(' ')

  assert header[:3]==['marker','alleleA','alleleB']

  def _load_beagle_dosage():
    import numpy as np

    m           = len(samples)
    n           = 3+m*3
    dosagevals  = np.array([0,0.5,1],dtype=float)

    for line_num,row in enumerate(gfile):
      fields = row.rstrip('\n').split()

      if len(fields) != n:
        raise ValueError('Invalid Beagle row on line %d of %s (found %d fields, expected %d)' \
                   % (line_num+2,namefile(filename),len(fields),n))

      lname   = fields[0]
      a,b     = fields[1:3]

      genos  = np.array(fields[3:],dtype=float).reshape( (-1,3) )
      dosage = (genos*dosagevals).sum(axis=1)/genos.sum(axis=1)

      yield lname,None,None,a,b,dosage

  genos = _load_beagle_dosage()

  return samples,genos


def test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  test()
