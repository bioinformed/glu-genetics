# -*- coding: utf-8 -*-

__abstract__  = 'GLU genotype input/output formats'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

__all__       = ['guess_informat','guess_informat_list','guess_outformat',
                 'get_genostream_loader','get_genostream_saver','get_genostream_writer',
                 'genostream_preferred_format']


from   glu.lib.fileutils import guess_format


########################################################################################################
# Internal data structures

INPUT_EXTS   = set()
OUTPUT_EXTS  = set()
LOADERS      = {}
SAVERS       = {}
WRITERS      = {}
PREF_FORMATS = {}


########################################################################################################
# Public API
def guess_informat(filename):
  '''
  @param filename: a file name or file object
  @type  filename: str or file object
  '''
  return guess_format(filename, INPUT_EXTS)


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
  return guess_format(filename, OUTPUT_EXTS)


def get_genostream_loader(format):
  try:
    return LOADERS[format]
  except KeyError:
    raise NotImplementedError("File format '%s' is not supported" % format)


def get_genostream_saver(format):
  try:
    return SAVERS[format]
  except KeyError:
    raise NotImplementedError("File format '%s' is not supported" % format)


def get_genostream_writer(format):
  try:
    return WRITERS[format]
  except KeyError:
    raise NotImplementedError("File format '%s' is not supported" % format)


def genostream_preferred_format(format):
  return PREF_FORMATS.get(format)


########################################################################################################
# Format discovery function

def discover_modules():
  import pkgutil

  path     = __path__
  prefix   = __name__ + '.'

  for (i,name,ispkg) in pkgutil.walk_packages(path,prefix=prefix):
    try:
      module = __import__(name,fromlist=['*'])
    except:
      continue

    # Skip all but explicitly indexable modules and packages
    if not getattr(module,'__genoformats__',False):
      continue

    # Strip module prefix
    if name.startswith(prefix):
      name = name[len(prefix):]

    formats = getattr(module,'__genoformats__')

    yield name,module


def discover_formats():
  from   glu.lib.utils import is_str

  for name,module in discover_modules():
    formats = module.__genoformats__

    for loader,saver,writer,pformat,aliases,extensions in formats:
      loader = getattr(module,loader,None) if loader else None
      saver  = getattr(module,saver, None) if saver  else None
      writer = getattr(module,writer,None) if writer else None

      aliases = aliases or []
      if is_str(aliases):
        aliases = [aliases]

      extensions = extensions or []
      if is_str(extensions):
        extensions = [extensions]

      allnames = set(aliases) | set(extensions)

      if loader:
        for name in extensions:
          INPUT_EXTS.add(name)
        for name in allnames:
          if name in LOADERS:
            raise ValueError('Conflicting genotype loader for format %s' % name)
          LOADERS[name] = loader

      if saver:
        for name in extensions:
          OUTPUT_EXTS.add(name)
        for name in allnames:
          if name in SAVERS:
            raise ValueError('Conflicting genotype saver for format %s' % name)
          SAVERS[name] = saver

      if writer:
        for name in extensions:
          OUTPUT_EXTS.add(name)
        for name in allnames:
          if name in WRITERS:
            raise ValueError('Conflicting genotype writer for format %s' % name)
          WRITERS[name] = writer


discover_formats()
