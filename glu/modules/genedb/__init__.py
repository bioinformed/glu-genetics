# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Genomic annotation database'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import os
import sys
import sqlite3


DEFAULT_GENEDB = ['genedb_ncbi36.3_dbsnp129',
                  'genedb_ncbi36.3_dbsnp128',
                  'genedb_ncbi36.3_huge',
                  'genedb_ncbi36.3_hapmap23a'
                  'genedb_ncbi36.3_gwas '
                  'genedb_ncbi36.3_tiny',
                  'genome36.3']

DEFAULT_PATHS  = [os.path.join(os.path.dirname(__file__),'data'),
                  '/usr/local/share/genedb',
                  '/usr/share/genedb']


def search_file(filename, path):
  '''
  Find a file given a search path
  '''
  if os.path.isfile(filename):
    return os.path.abspath(filename)

  if not path:
    return None

  for pdir in path:
    if pdir and os.path.isdir(pdir):
      name = os.path.join(pdir, filename)
      if os.path.isfile(name):
        return os.path.abspath(name)

  return None


def open_genedb(dbname=None,path=None):
  '''
  Open a genedb instance based on an optional database name and search path.
  If not specified, a series of standard database names and paths will be
  explored.
  '''
  path  = os.environ.get('GLU_GENEDB_PATH','').split(os.pathsep)
  path += DEFAULT_PATHS

  dbnames = []
  if dbname:
    dbnames += [dbname]
  else:
    if 'GLU_GENEDB' in os.environ:
      dbnames += [os.environ['GLU_GENEDB']]
    dbnames += DEFAULT_GENEDB

  for name in dbnames:
    if not name.endswith('.db'):
      name += '.db'

    filename = search_file(name, path)

    if filename:
      sys.stderr.write('[INFO] Opening genedb: %s\n' % filename)
      return sqlite3.connect(filename)

  raise IOError('Cannot open genedb')
