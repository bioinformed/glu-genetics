# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Print a list of GLU modules'
__index__     = 'General modules'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys
import pydoc
import pkgutil
import textwrap

try:
  from cStringIO import StringIO
except ImportError:
  from StringIO import StringIO

import glu
import glu.modules


def main():
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

    # Add only module indexed modules
    if ismodule and index:
      order = getattr(module,'__order__', 50)
      orders[parts] = min(orders.get(parts,100),order)
      modules.append( (index,parts,name,module) )

  # FIXME: Figure out what to do when width is too large
  width   = max(len(m[2]) for m in modules)+4
  modules = sorted((orders.get(m[0],50),orders.get(m[1],50))+m for m in modules)

  out = StringIO()
  out.write('Welcome to GLU: The Genotype Library and Utilities\n')

  inindex = None
  for o11,o2,index,parts,name,module in modules:
    if index != inindex:
      inindex = index
      out.write('\n%s:\n\n' % index)

    abstract = getattr(module,'__abstract__','')

    initial_width = max(0,width-len(name)-4)
    text = textwrap.wrap(abstract, 78, initial_indent=' '*initial_width,
                                       subsequent_indent=' '*width)
    out.write('  %s  %s\n' % (name,'\n'.join(text)))

  pydoc.pager(out.getvalue())
