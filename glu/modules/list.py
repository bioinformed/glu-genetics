# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Print a list of GLU modules'
__index__     = 'General modules'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import pydoc
import pkgutil
import textwrap

try:
  from cStringIO import StringIO
except ImportError:
  from StringIO import StringIO

import glu.modules


def find_glu_modules():
  path     = glu.modules.__path__
  prefix   = glu.modules.__name__ + '.'
  descrs   = {():glu.modules.__abstract__}
  orders   = {():10}
  modules  = []

  for (i,name,ispkg) in pkgutil.walk_packages(path,prefix=prefix):
    try:
      module = __import__(name,fromlist=['main'])
    except:
      continue

    # Skip all but explicitly indexable modules and packages
    if not getattr(module,'__gluindex__',False):
      continue

    # Strip module prefix
    if name.startswith(prefix):
      name = name[len(prefix):]

    # Construct package path and determine base package
    parts = tuple(name.split('.'))
    base  = parts[:-1]

    # Do not index modules unless their base is indexable
    if base not in descrs:
      continue

    # Get attributes
    ismodule = hasattr(module,'main')
    index    = getattr(module,'__index__',None)
    order    = getattr(module,'__order__',50)

    # Use abstract text for packages if no index name is available
    if not index and ispkg and not ismodule:
      index = getattr(module,'__abstract__',None)

    # Index location is provided by the module
    if index:
      orders[index] = min(orders.get(index,100),order)
      descrs[parts] = index

    # Get index location from the base package
    else:
      index = descrs[base]

    # Add only indexed modules
    if ismodule and index:
      order = getattr(module,'__order__', 50)
      orders[parts] = min(orders.get(parts,100),order)
      modules.append( (index,parts,name,module) )

  modules = sorted((orders.get(m[0],50),orders.get(m[1],50))+m for m in modules)
  modules = [ (m[2],)+m[4:] for m in modules ]
  return modules


def main():
  modules = find_glu_modules()

  # FIXME: Figure out what to do when width is too large
  width   = max(len(m[1]) for m in modules)+4

  out = StringIO()
  out.write('Welcome to GLU: The Genotype Library and Utilities\n')

  inheading = None
  for heading,name,module in modules:
    if heading != inheading:
      inheading = heading
      out.write('\n%s:\n\n' % heading)

    abstract = getattr(module,'__abstract__','')

    initial_width = max(0,width-len(name)-4)
    text = textwrap.wrap(abstract, 78, initial_indent=' '*initial_width,
                                       subsequent_indent=' '*width)
    out.write('  %s  %s\n' % (name,'\n'.join(text)))

  pydoc.pager(out.getvalue())
