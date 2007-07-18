# -*- coding: utf-8 -*-
'''
File:          completion.py

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       2006-06-29

Abstract:      Performs completion analysis on genotype data

Requires:      Python 2.5, biozilla

Revision:      $Id:
'''

__version__ = '0.2'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import sys
import csv
from   itertools           import islice, izip
from   textwrap            import fill
from   biozilla.utils      import percent
from   biozilla.genodata   import load_genomatrix
from   biozilla.sections   import save_section, SectionWriter, save_metadata_section


# FIXME: Use an output stream (-o)
def completion_summary(rowlabel,collabel,results,rowresults,colresults):
  # FIXME: This is beginning to look like a closure
  inf,genos_all,genos_inf = results
  rowcomp,emptyrows,droppedrows,totalrows = rowresults
  colcomp,emptycols,droppedcols,totalcols = colresults

  print '%s     Total: %7d, Empty: %7d, Dropped: %7d, Informative: %7d' % \
             (rowlabel,totalrows,len(emptyrows),droppedrows,len(rowcomp)-len(emptyrows))
  print '%s  Total: %7d, Empty: %7d, Dropped: %7d, Informative: %7d' % \
             (collabel,totalcols,len(emptycols),droppedcols,len(colcomp)-len(emptycols))
  print
  print 'GLOBAL GENOTYPE COMPLETION RATE FOR ALL DATA:         %10d / %10d = %5.3f%%' % (inf,genos_all,percent(inf,genos_all))
  print 'GLOBAL GENOTYPE COMPLETION RATE FOR INFORMATIVE DATA: %10d / %10d = %5.3f%%' % (inf,genos_inf,percent(inf,genos_inf))
  e = ' '*5
  empty = e+'(empty)'
  print
  print '%s with no data:' % rowlabel
  print fill(', '.join(emptyrows), initial_indent=e,subsequent_indent=e) or empty
  print
  print '%s with no data:' % collabel
  print fill(', '.join(emptycols), initial_indent=e,subsequent_indent=e) or empty
  print


# FIXME: Use an output stream (-o)
def completion_output(name, data, total_inf, total_uninf):
  data = sorted( (n,k) for k,(n,m) in data.iteritems() )
  print 'MISSING GENOTYPES BY %s' % name.upper()
  print
  print '                                        Informative                  All'
  print '  Rank  %-25s        N / Total      %%           N / Total      %%  ' % name
  print '  ----  -------------------------  -----------------  ------  -----------------  ------'

  total_all = total_inf+total_uninf
  for r,(n,l) in enumerate(data):
    data = (r+1,l,n,total_inf,percent(n,total_inf),n,total_all,percent(n,total_all))
    print '  %4d  %-25s  %7d / %7d  %6.2f  %7d / %7d  %6.2f' % data
  print


# FIXME: Use an output stream (-o)
def completion_output2(name, data, m):
  data = sorted( (v,k) for k,v in data.iteritems() )
  e = ' '*43
  rank = 1

  print 'MISSING GENOTYPES BY %s' % name.upper()
  print
  print '  Ranks              # Missing        %%    Missing %s names' % name
  print '  -------------  -----------------  -----  -------------------------------------------------'
  for items in data:
    k = len(items)
    if k > 9:
      it = items[:9] + ['...']
    else:
      it = items
    text = fill(', '.join(it), initial_indent=e,subsequent_indent=e,width=99).lstrip(' ')
    p = percent(n,m)
    print '  %6d-%6d  %7d / %7d  %4.1f%%  %s' % (rank,rank+k-1,n,m,p,text)
    rank += k

  print


def completion(genos, droppedrows=0, droppedcols=0):
  header = list(genos.next())

  print >> sys.stderr, 'Computing completion rates...',

  rowcomp = {}
  colcomp = {}

  colcomps = [ colcomp.setdefault(colkey,[0,0]) for colkey in header ]

  for row in genos:
    rowkey,rgenos = row
    items = izip(colcomps,rgenos)

    rcomp = rowcomp.setdefault(rowkey,[0,0])
    for ccomp,geno in items:
      # Dirty python trick where not missing == 0, missing == 1
      g = geno==0
      rcomp[g] += 1
      ccomp[g] += 1

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

  # FIXME: This is beginning to look like a closure
  return (inf,genos_all,genos_inf),                 \
         (rowcomp,emptyrows,droppedrows,totalrows), \
         (colcomp,emptycols,droppedcols,totalcols)

def save_completion_summary(sw, results, section_type):
 section='summary'
 comp,empty,dropped,total = results
 data = [['slice',        section_type],
         ['total',       total],
         ['empty',       ','.join(empty)],
         ['dropped',     dropped],
         ['informative', len(comp)-len(empty)]]

 save_section(sw, section, data)


def save_completion_results(sw, results, section_type):
 section='data'
 data = [[name,complete] for name,(complete,total) in results.iteritems()]
 data = [['slice', section_type], ['id', 'complete']] + data
 save_section(sw, section, data)


def option_parser():
  import optparse

  usage = 'usage: %prog [options] file'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-o', '--output', dest='output', metavar='FILE',
                    help='Output of duplicate check report')
  parser.add_option('-l', '--rowlimit', dest='rowlimit', metavar='N', type='int', default=0,
                    help='Limit the number of rows considered to N for testing purposes (default=0 for unlimited)')
  parser.add_option('-L', '--collimit', dest='collimit', metavar='N', type='int', default=0,
                    help='Limit the number of columns considered to N for testing purposes (default=0 for unlimited)')
  parser.add_option('-r', '--droppedrows', dest='droppedrows', metavar='N', type='int', default=0,
                    help='Number of rows that where dropped from the dataset previously.  Used to compute overall completion.')
  parser.add_option('-c', '--droppedcols', dest='droppedcols', metavar='N', type='int', default=0,
                    help='Number of columns that where dropped from the dataset previously.  Used to compute overall completion.')
  parser.add_option('-f', '--format', dest='format', metavar='NAME', type='string', default='',
                    help='Format of the input data. Values=sdat,ldat')
  parser.add_option('--tabularoutput', dest='tabularoutput', metavar='FILE', type='string',
                    help='Generate machine readable tabular output of results')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if len(args) != 1:
    parser.print_help()
    return

  format,genos = load_genomatrix(args[0],format=options.format,limit=options.collimit or None)

  format = format or options.format

  if format=='sdat':
    rowlabel='Sample'
    collabel='Loci'
  elif format=='ldat':
    rowlabel='Loci'
    collabel='Samples'
  else:
    rowlabel='Rows'
    collabel='Columns'

  if options.rowlimit:
    genos = islice(genos, options.rowlimit+1)

  results,rowresults,colresults = completion(genos, options.droppedrows, options.droppedcols)

  print >> sys.stderr, 'Writing completion output...',

  completion_summary(rowlabel,collabel,results,rowresults,colresults)

  # FIXME: This is beginning to look like a closure
  rowcomp,emptyrows,droppedrows,totalrows = rowresults
  colcomp,emptycols,droppedcols,totalcols = colresults

  completion_output(rowlabel, rowcomp, len(colcomp)-len(emptycols), len(emptycols)+droppedcols)
  completion_output(collabel, colcomp, len(rowcomp)-len(emptyrows), len(emptyrows)+droppedrows)

  if options.tabularoutput:
    sw = SectionWriter(options.tabularoutput)
    save_metadata_section(sw, analysis='completion', analysis_version=__version__)
    save_completion_summary(sw,rowresults,rowlabel)
    save_completion_summary(sw,colresults,collabel)
    save_completion_results(sw,rowcomp,rowlabel)
    save_completion_results(sw,colcomp,collabel)

  print >> sys.stderr, 'Done.'


if __name__ == '__main__':
  main()
