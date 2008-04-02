# -*- coding: utf-8 -*-
'''
File:          completion.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs completion analysis on genotype data

Requires:      Python 2.5, glu

Revision:      $Id:
'''

__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import csv

from   itertools         import izip
from   textwrap          import fill
from   collections       import defaultdict

from   glu.lib.utils     import percent
from   glu.lib.fileutils import autofile,hyphen
from   glu.lib.genolib   import load_genostream
from   glu.lib.sections  import save_section, SectionWriter, save_metadata_section


class MatrixResults(object):
  def __init__(self,rowresults,colresults,completed,total_all,total_inf):
    self.rowresults = rowresults
    self.colresults = colresults
    self.completed  = completed
    self.total_all  = total_all
    self.total_inf  = total_inf


class SliceResults(object):
  def __init__(self,label,completion,empty,dropped,total):
    self.label      = label
    self.completion = completion
    self.empty      = empty
    self.dropped    = dropped
    self.total      = total


def completion_summary(outfile,results):
  for x in (results.rowresults,results.colresults):
    outfile.write('%-7s Total: %7d, Empty: %7d, Dropped: %7d, Informative: %7d\n' % \
               (x.label,x.total,len(x.empty),x.dropped,len(x.completion)-len(x.empty)))
  outfile.write('\n')
  outfile.write('GLOBAL GENOTYPE COMPLETION RATE FOR ALL DATA:         %10d / %10d = %5.3f%%\n' %
                 (results.completed,results.total_all,percent(results.completed,results.total_all)))
  outfile.write('GLOBAL GENOTYPE COMPLETION RATE FOR INFORMATIVE DATA: %10d / %10d = %5.3f%%\n' %
                 (results.completed,results.total_inf,percent(results.completed,results.total_inf)))

  e = ' '*5
  empty = e+'(empty)'
  outfile.write('\n')
  for x in (results.rowresults,results.colresults):
    outfile.write('%s with no data:\n' % x.label)
    outfile.write(fill(', '.join(x.empty), initial_indent=e,subsequent_indent=e) or empty)
    outfile.write('\n\n')


def completion_output(outfile, results1, results2):
  data = sorted( (n,k) for k,(n,m) in results1.completion.iteritems() )
  outfile.write('MISSING GENOTYPES BY %s\n' % results1.label.upper())
  outfile.write('\n')
  outfile.write('                                        Informative                  All\n')
  outfile.write('  Rank  %-25s        N / Total      %%           N / Total      %%  \n' % results1.label)
  outfile.write('  ----  -------------------------  -----------------  ------  -----------------  ------\n')

  total_inf = len(results2.completion)-len(results2.empty)
  total_all = total_inf+len(results2.empty)+results2.dropped

  for r,(n,l) in enumerate(data):
    data = (r+1,l,n,total_inf,percent(n,total_inf),
                  n,total_all,percent(n,total_all))
    outfile.write('  %4d  %-25s  %7d / %7d  %6.2f  %7d / %7d  %6.2f\n' % data)
  outfile.write('\n')


def completion(genos, rowlabel, collabel, droppedrows=0, droppedcols=0):
  print >> sys.stderr, 'Computing completion rates...',

  rowcomp  = defaultdict(lambda: [0,0])
  colcomp  = defaultdict(lambda: [0,0])
  colcomps = [ colcomp[colkey] for colkey in genos.columns ]

  for rowkey,rgenos in genos:
    rcomp = rowcomp[rowkey]

    for ccomp,geno in izip(colcomps,rgenos):
      # Dirty python trick where not missing == 0, missing == 1
      missing = not geno
      rcomp[missing] += 1
      ccomp[missing] += 1

  print >> sys.stderr, 'Done.'

  # Informative = mixture of missing and non-missing genotypes for each row
  #               and column
  # Empty       = Observed whole rows or columns of missing genotypes
  # Dropped     = Empty rows or columns not included in the observed data
  #
  # ==================|=============|------------
  # ||                |             |           |
  # ||  Informative   |<-EmptyCols->|           |
  # ||                |             |<-Dropped->|
  # ||----------------|-------------|   Cols    |
  # ||      ^         |             |           |
  # ||   EmptyRows    |  EmptyBoth  |           |
  # ||      v         |             |           |
  # ||----------------|-------------|-----------|
  # |                 ^             |           |
  # |             DroppedRows       |  Dropped  |
  # |                 v             |   Both    |
  # |-----------------|-------------|-----------|

  emptyrows = set(key for key,counts in rowcomp.iteritems() if not counts[0])
  emptycols = set(key for key,counts in colcomp.iteritems() if not counts[0])

  missingrows = len(emptyrows) + droppedrows
  missingcols = len(emptycols) + droppedcols

  totalrows = len(rowcomp) + droppedrows
  totalcols = len(colcomp) + droppedcols

  empty   = len(emptyrows)*len(colcomp) + len(emptycols)*len(rowcomp) - len(emptyrows)*len(emptycols)
  missing = totalrows*droppedcols + totalcols*droppedrows - droppedcols*droppedrows

  if len(rowcomp) < len(colcomp):
    comp = rowcomp
  else:
    comp = colcomp

  inf,noninf = map(sum,izip(*comp.itervalues()))
  genos_inf = inf + noninf - empty
  genos_all = inf + noninf + missing

  rowresults = SliceResults(rowlabel,rowcomp,emptyrows,droppedrows,totalrows)
  colresults = SliceResults(collabel,colcomp,emptycols,droppedcols,totalcols)

  return MatrixResults(rowresults,colresults,inf,genos_all,genos_inf)


def save_completion_summary(sw, results):
 data = [['slice',       results.label],
         ['total',       results.total],
         ['empty',       ','.join(results.empty)],
         ['dropped',     results.dropped],
         ['informative', len(results.completion)-len(results.empty)]]

 save_section(sw, 'summary', data)


def save_completion_results(sw, results):
 data  = [['slice', results.label], ['id', 'complete']]
 data += [[name,complete] for name,(complete,total) in results.completion.iteritems()]
 save_section(sw, 'data', data)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-f', '--format', dest='format', metavar='NAME',
                    help='Format of the input data. Values=sdat,ldat')
  parser.add_option('-g', '--genorepr',        dest='genorepr',        metavar='REP',
                    help='Input genotype representations')
  parser.add_option('-l', '--loci', dest='loci', metavar='FILE',
                    help='Locus description file and options')
  parser.add_option('-o', '--output', dest='output', metavar='FILE', default='-',
                    help='Output of completion report')
  parser.add_option('-r', '--droppedrows', dest='droppedrows', metavar='N', type='int', default=0,
                    help='Number of rows that where dropped from the dataset previously.  Used to compute overall completion.')
  parser.add_option('-c', '--droppedcols', dest='droppedcols', metavar='N', type='int', default=0,
                    help='Number of columns that where dropped from the dataset previously.  Used to compute overall completion.')
  parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE',
                    help='Generate machine readable tabular output of results')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  genos = load_genostream(args[0],format=options.format,genorepr=options.genorepr,
                                  genome=options.loci)

  if genos.format not in ('sdat','ldat'):
    genos = genos.as_ldat()

  if genos.format=='sdat':
    rowlabel='Sample'
    collabel='Locus'
  elif genos.format=='ldat':
    rowlabel='Locus'
    collabel='Sample'

  outfile = autofile(hyphen(options.output, sys.stdout),'w')

  results = completion(genos,rowlabel,collabel,options.droppedrows,options.droppedcols)

  print >> sys.stderr, 'Writing completion output...',

  completion_summary(outfile,results)
  completion_output(outfile, results.rowresults, results.colresults)
  completion_output(outfile, results.colresults, results.rowresults)

  if options.tabularoutput:
    sw = SectionWriter(options.tabularoutput)
    save_metadata_section(sw, analysis='completion', analysis_version='0.1')
    save_completion_summary(sw,results.rowresults)
    save_completion_summary(sw,results.colresults)
    save_completion_results(sw,results.rowresults)
    save_completion_results(sw,results.colresults)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
