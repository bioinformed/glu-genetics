# -*- coding: utf-8 -*-

__abstract__  = 'file related utility functions'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import os

from   glu.lib.utils              import is_str
from   glu.lib.fileutils.parser   import parse_augmented_filename


__all__ = ['autofile','namefile','hyphen','guess_format','related_file','guess_related_file']


COMPRESSED_SUFFIXES = {'gz':'gzip', 'Z'  :'gzip',
                       'bz':'bzip2','bz2':'bzip2'}


def spawn_compressor(binary, filename, mode, bufsize=-1):
  '''
  Spawn a subprocess to run a compressor like gzip or bzip2 and connect to
  input/output pipes

  @param    binary: executable name
  @type     binary: str
  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @param   bufsize: buffering mode and size.  0=unbuffered, 1=linebuffered, >1 buffer size,
                    -1 default buffering (default)
  @type    bufsize: int
  @return         : file object to read from or write to
  @rtype          : file object
  '''
  if os.environ.get('GLU_NOSPAWN'):
    raise OSError('Spawning external processes disabled by GLU_NOSPAWN')

  from subprocess import Popen,PIPE

  if 'w' in mode:
    out = file(filename,mode)
    cmd = [binary,'-c']
    f   = Popen(cmd, stdin=PIPE, stdout=out, bufsize=bufsize).stdin
  else:
    cmd = [binary,'-d','-c',filename]
    f   = Popen(cmd, stdout=PIPE, universal_newlines='U' in mode, bufsize=bufsize).stdout

  return f


def hyphen(filename,defaultfile,args=None):
  '''
  Return the default if input is '-', otherwise itself

  @param  defaultfile: default file name or file object
  @type   defaultfile: str or file object
  @param     filename: file name or file object
  @type      filename: str or file object
  @param         args: augmented argument output dictionary or None
  @type          args: dict or None
  @return            : file name or file object to read from or write to
  @rtype             : str or file object
  '''
  if not is_str(filename):
    return filename

  if args is None:
    args = {}

  if parse_augmented_filename(filename,args) == '-':
    return defaultfile

  return filename


def namefile(filething):
  '''
  Return a human-comprehensible name for the file or file object provided.
  Recognizes file objects with the 'name' attribute, including sys.stdin,
  sys.stdout, and sys.stderr.

  @param  filething: file or filename
  @type   filething: file object or string
  @return:           human-comprehensible file name

  >>> import sys
  >>> namefile(file('/dev/null'))
  '/dev/null'
  >>> namefile(sys.stdin)
  '<stdin>'
  >>> namefile(sys.stderr)
  '<stderr>'
  >>> namefile('/dev/null')
  '/dev/null'
  '''
  if is_str(filething):
    return filething
  elif getattr(filething, 'name', None) is not None:
    return filething.name
  else:
    return repr(filething)


_autofile_errors = (ImportError,ValueError,OSError)


def compressed_filename(filename):
  '''
  Determine if the input file is in or needs to be in compressed format

  @param  filename: file name or file object
  @type   filename: str or file object
  @return         : compression format, if compressed, otherwise an empty string
  @rtype          : str

  >>> compressed_filename('subjects.sdat')
  ''
  >>> compressed_filename('subjects.sdat.gz')
  'gzip'
  >>> compressed_filename('../subjects.sdat')
  ''
  >>> compressed_filename('../subjects.sdat.gz')
  'gzip'
  >>> compressed_filename('../subjects.sdat.bz2')
  'bzip2'
  '''
  if not is_str(filename):
    return ''

  filename = os.path.expanduser(filename)
  parts    = os.path.basename(filename).split('.')
  ext      = parts[-1] if parts else ''
  return COMPRESSED_SUFFIXES.get(ext,'')


def autofile(filename, mode='r', bufsize=-1):
  '''
  Return a file object in the correct compressed format as specified, which
  is ready to read from or write to

  @param  filename: file name or file object
  @type   filename: str or file object
  @param      mode: determine whether the file objects should be opened for input or output,
                    either 'w' or 'r'
  @type       mode: str
  @return         : file object to read from or write to
  @rtype          : file object
  '''
  # Pass non-string filename objects back as file-objects
  if not is_str(filename):
    return filename

  filename = os.path.expanduser(filename)
  comp     = compressed_filename(filename)

  if not comp:
    f = file(filename, mode)
  elif comp == 'gzip':
    try:
      f = spawn_compressor(os.environ.get('GLU_GZIP','gzip'), filename, mode, bufsize=bufsize)
    except _autofile_errors:
      import gzip
      f = gzip.GzipFile(filename, mode)
  elif comp == 'bzip2':
    try:
      f = spawn_compressor(os.environ.get('GLU_BZIP2','bzip2'), filename, mode, bufsize=bufsize)
    except _autofile_errors:
      import bz2
      f = bz2.BZ2File(filename, mode, buffering=max(0,bufsize))

  return f


def guess_format(filename, formats, args=None):
  '''
  Return the guessed format indicated by the filename itself

  @param  filename: file name
  @type   filename: str
  @param    format: any expected file format
  @type     format: list of strs
  @return         : file format guessed from the filename
  @rtype          : str

  >>> f = ['sdat','ldat']
  >>> guess_format('subjects.sdat', f)
  'sdat'
  >>> guess_format('subjects.sdat.gz', f)
  'sdat'
  >>> guess_format('../subjects.sdat', f)
  'sdat'
  >>> guess_format('../subjects.sdat.gz', f)
  'sdat'
  >>> guess_format('../subjects.sdat.gz:format=ldat', f)
  'ldat'
  '''

  if not is_str(filename):
    return None

  # Parse to remove augmented arguments
  if args is None:
    args = {}

  filename = parse_augmented_filename(filename,args)

  if 'format' in args:
    return args['format']

  parts = os.path.basename(filename).split('.')

  if len(parts)>1 and parts[-1] in COMPRESSED_SUFFIXES:
    parts.pop()

  if parts and parts[-1] in formats:
    return parts[-1]

  return None


def related_file(filename,extension):
  '''
  Return a filename with the extension provided

  @param   filename: base filename with extension
  @type    filename: str
  @param extensions: new extension (without a '.')
  @type  extensions: str
  @return          : new filename
  @rtype           : str

  >>> related_file('foo.ldat', 'sdat')
  'foo.sdat'
  >>> related_file('foo', 'sdat')
  'foo.sdat'
  '''
  if not is_str(filename):
    return None

  parts = os.path.basename(filename).split('.')

  if len(parts)>1 and parts[-1] in COMPRESSED_SUFFIXES:
    parts.pop()

  if not parts:
    raise ValueError('invalid filename')

  if len(parts)>1:
    parts.pop()

  parts.append(extension)

  return '.'.join(parts)


def guess_related_file(filename,extensions):
  '''
  Find a related file with the same name except different prefix.  Only
  files that exist are returned.

  @param   filename: base filename with extension
  @type    filename: str
  @param extensions: list of alternative extensions (without a '.')
  @type  extensions: list of str
  @return          : new filename
  @rtype           : str

  >>> guess_related_file('fileutils.dat',['py']) # doctest:+SKIP
  'fileutils.py'
  '''
  if not is_str(filename):
    return None

  prefix,ext = os.path.splitext(filename)

  if not prefix:
    raise ValueError('invalid filename')

  for new_ext in extensions:
    testfile = '%s.%s' % (prefix,new_ext)
    if os.path.isfile(testfile):
      return testfile

  return None


def _test():
  import doctest
  return doctest.testmod()


if __name__ == '__main__':
  _test()
