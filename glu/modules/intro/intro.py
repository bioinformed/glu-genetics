# -*- coding: utf-8 -*-

__gluindex__  = False
__abstract__  = 'Brief introduction to GLU for new users'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

intro = '''\
Welcome to GLU: The Genotype Library and Utilities.

GLU provides tools for the management of large amounts of SNP genotype data
and programs to check its quality and to test for association between SNP
markers with continuous or discrete trait phenotypes.

Run the following modules for more information:

  intro.overview   : More detailed overview of GLU
  intro.formats    : Basic information on GLU file formats
  intro.quickstart : Just the basics for how to get started
  list             : List major modules
'''

def main():
  import pydoc
  pydoc.pager(intro)
