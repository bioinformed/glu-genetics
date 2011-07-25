# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Dump the contents of an Illumina BPM Manifest file to a text file'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   glu.lib.fileutils import table_writer
from   glu.lib.illumina  import IlluminaManifest


def option_parser():
  import optparse

  usage = 'usage: %prog [options] manifest.bpm'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output manifest information')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help(sys.stderr)
    sys.exit(2)

  manifest = IlluminaManifest(args[0])
  rows     = iter(manifest)
  header   = rows.next()

  out = table_writer(options.output,hyphen=sys.stdout)
  out.writerow(header)
  out.writerows(rows)


if __name__=='__main__':
  main()
