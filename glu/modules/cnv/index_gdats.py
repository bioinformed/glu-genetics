# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Create an index of one or more GDAT files'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import sys

from   glu.modules.cnv.gdat import GDATIndex, GDATFile


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('index',            help='Index file to create or update')
  parser.add_argument('gdat',  nargs='+', help='GDAT file(s)')

  return parser


def main():
  parser    = option_parser()
  options   = parser.parse_args()
  indexfile = GDATIndex(options.index)

  for filename in options.gdat:
    gdat     = GDATFile(filename)
    manifest = gdat.attrs['ManifestName']

    print 'GDAT: %s (snps=%d, samples=%d, manifest=%s)' % (filename,len(gdat.snps),len(gdat.samples),manifest)

    indexfile.index(gdat)


if __name__=='__main__':
  main()
