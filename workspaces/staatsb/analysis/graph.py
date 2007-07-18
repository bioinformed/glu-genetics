# -*- coding: utf-8 -*-
'''
File:          graph.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-11-07

Requires:      Python 2.5, biozilla

Revision:      $Id: $
'''
__version__   = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'
REPORTING     = 'Reporting v%s(BETA)' % __version__

import  os
import  csv
import  sys
import  time

from    operator                  import itemgetter
from    contextlib                import contextmanager
from    itertools                 import islice,chain,groupby,izip,repeat

from    reportlab.platypus        import Paragraph
from    reportlab.pdfgen          import canvas
from    reportlab.lib             import colors, enums as raenums
from    reportlab.platypus        import Spacer

from    biozilla.utils            import autofile, percent
from    biozilla.sections         import read_sections,filter_sections,index_sections
from    biozilla.stats            import quantile

import  styles

# TODO: list each instead of *
from    reportlib                 import *
from    graphics                  import draw_lineplot


# TODO: find a better way of accessing the first page from the doc

# page template indexes in the document
FIRST_PAGE=0
LATER_PAGE=1


class GraphPageStyle(object):
  # Why does this object need to be smart enough to infer the title?  ie,
  # add title as an argument, drop meta until it is needed
  def __init__(self, meta, options):
    self.title = 'Unknown Graph'
    if meta['analysis'] == 'completion':
      self.title = 'Genotype Completion'

  def __call__(self, canvas, doc):
    titleframe      = doc.pageTemplates[FIRST_PAGE].title
    subtitleframe   = doc.pageTemplates[FIRST_PAGE].subtitle
    contentframe    = doc.pageTemplates[FIRST_PAGE].content

    titleStyle      = styles.Style(fontsize=24, fontname='Helvetica-Bold', fontcolor=colors.white)
    titleStyle      = titleStyle.paragraphstyle('titleStyle')
    subtitleStyle   = styles.Style(fontsize=10, fontname='Helvetica-Bold', fontcolor=colors.white)
    subtitleStyle   = subtitleStyle.paragraphstyle('subtitleStyle', alignment=styles.TA_RIGHT)

    titleframe.add(styles.p(self.title, style=titleStyle), canvas)
    subtitle =  [styles.p(REPORTING, style=subtitleStyle),styles.p(time.asctime(), style=subtitleStyle)]
    subtitleframe.addFromList(subtitle, canvas)


# Need docstring and doctests
def check_sections(indexed_sections, expected_sections):
  for section,mins,maxs in expected_sections:
    found = indexed_sections.get(section, [])

    if not (mins <= len(found) <= maxs):
      raise ValueError, 'Unexpected number of section bodies found for: ' + section

    if maxs == 1:
      indexed_sections[section] = found[0] if found else None

  return indexed_sections


class CompletionData(object):
  def __init__(self, summary, completion, empty, dropped):
    self.summary      = summary     # { 'tape':'samples', 'completed': 6, 'total': 9 }
    self.completion   = completion  # [ ('id', [members#, empty#, dropped#, completed#, total#]) ]
    self.empty        = empty       # [ 'rs3781', 'rs7410', 'rs27180', 'rs10028']
    self.dropped      = dropped     # [ 'rs3781', 'rs7410', 'rs27180', 'rs10028']


class CompletionModel(object):
  def __init__(self, metadata, summary):
    self.metadata = metadata
    self.summary  = summary
    self.data     = []

  def append(self, data):
    '''
    Append a CompletionData object to the list of CompletionData objects.
    '''
    self.data.append(data)

  def find_types(self, ctype):
    return [ c for c in self.data if c.summary.get('type') == ctype ]

  def find_type(self, ctype):
    found = self.find_types(ctype)
    if len(found) != 1:
      raise KeyError, 'Found more than one type of data for ' + ctype

    return found[0]

  def __iter__(self):
    return iter(self.data)

  def __contains__(self, ctype):
    return any(self.find_types(ctype))


def identity(x):
  return x


def process_sections(sections, intfields=[], floatfields=[]):
  key,ctype = section[0]
  assert key=='type'
  header = section[1]

  fieldmap = defaultdict(lambda: identity)

  for f in intfields:
    fieldmap[f] = int

  for f in floatfields:
    fieldmap[f] = float

  section = [ [ fieldmap[h](v) for h,v in izip(header,row) ] for row in islice(section,2,None) ]

  return ctype,section


def parse_completion(completion):
  data      = {}
  summaries = {}
  groups    = {}
  empty     = {}
  dropped   = {}

  for section in completion['data']:
    ctype,section = process_sections(section, ('members#','empty#','dropped#','completed#','total#'))
    data[ctype] = section

  for section in completion['summary']:
    ctype,section = process_sections(section)
    summaries[ctype] = dict(section)

  for section in completion['group']:
    ctype,section = process_sections(section, ('id','members#','empty#','dropped#','completed#','total#'))
    groups[ctype] = section

  for section in completion['empty']:
    key,value = section[0]
    assert key=='type'
    empty[value] = section[1:]

  for section in completion['dropped']:
    key,value = section[0]
    assert key=='type'
    dropped[value] = section[1:]

  assert len(summaries) >= len(data) == len(empty) == len(dropped)

  meta = dict(completion['header'])

  model = CompletionModel(meta, summaries['global'])

  # Append to model, not list
  for k in data:
    model.append( CompletionData(summaries[k], data[k], empty[k], dropped[k]) )

  return model


def parse_completion(file):
  data     = csv.reader(autofile(file), dialect='excel-tab')
  sections = index_sections(read_sections(data))

  expected_sections = [('header',1,1)]
  data              = check_sections(sections, expected_sections)
  meta              = dict(data['header'])

  datatype = meta['analysis']

  if datatype != 'completion':
    raise NotImplementedError, datatype + ' data type is not yet supported.'

  expected_sections = [('summary',  3,3),
                       ('group',    2,2),
                       ('data',     2,2),
                       ('empty',    2,2),
                       ('dropped',  2,2)]
  data  = check_sections(data, expected_sections)
  model = parse_completion(data)

  return meta,model


def create_content(model, page, options):
  spacewidth  = 1
  spaceheight = 100
  content     = []
  frame       = page.content

  if options.graphtype == 'quantile':
    style       = create_quantile_style(options.graphstyle)
    samplecolor = style.data[0].fontcolor
    locuscolor  = style.data[1].fontcolor
    # title could be an option
    text        = 'Genotype Assay Completion by Quantiles for <font color=%s>Samples</font> and <font color=%s>Loci</font>' % (samplecolor,locuscolor)
    content.append(create_graph_title(text, style))
    content.append(create_quantile_graph(model, style, frame, options))

    # content.append(Spacer(spacewidth, spaceheight))

  return content


def create_quantile_style(graphstyle):
  if graphstyle == 'nci':
    from styles import ncicolors

    line1 = styles.Style(strokewidth=1, strokecolor=ncicolors['darkred'],   fontcolor=ncicolors['darkred'])
    line2 = styles.Style(strokewidth=1, strokecolor=ncicolors['darkblue'],  fontcolor=ncicolors['darkblue'])
    axis  = styles.Style(fontsize=8,    fontcolor=ncicolors['black'])
    title = styles.Style(fontsize=14,   fontcolor=ncicolors['black'])
    style = styles.GraphStyle([line1,line2], axis, axis, title)
  else:
    raise NotImplementedError, graphstyle + ' style type is not yet supported.'

  return style


def create_graph_title(text, style):
  titlestyle = style.title.paragraphstyle('Title', alignment=raenums.TA_CENTER)
  return Paragraph(text, style=titlestyle)


def create_quantile_graph(model, style, frame, options, step_count=1000):
  width   = frame.width-(2*SPACER)
  height  = frame.height/3

  quantiles       = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1]
  quantile_lables = ['%.02f' % q for q in quantiles ]

  xsteps     = [float(i)/step_count for i in xrange(step_count)]

  xaxis      = Axis(0, step_count, 1, quantile_labels, name='Quantiles')
  yaxis      = Axis(0,100,10,name='Percent Completion')

  # [ ('id', [members#, empty#, dropped#, completed#, total#]) ]
  samcomp = sorted(values[3] for ind,values in model.find_type('samples').completion)
  loccomp = sorted(values[3] for ind,values in model.find_type('loci').completion)

  nsam    = loci.summary['total']
  nloc    = samples.summary['total']
  samperc = [ percent(quantile(samcomp,x,presorted=True),nsam) for x in xsteps ]
  locperc = [ percent(quantile(loccomp,x,presorted=True),nloc) for x in xsteps ]

  samtable = ['Sample']+[ '%6.2f' % percent(quantile(samcomp,x,presorted=True),nsam) for x in quantiles ]
  loctable = ['Locus'] +[ '%6.2f' % percent(quantile(loccomp,x,presorted=True),nloc) for x in quantiles ]
  table = zip(samtable,loctable)

  # order matters, since plot lines are color coded
  points = enumerate(samperc),enumerate(locperc)
  return draw_lineplot(points, xaxis, yaxis, style, width, height, table)


def option_parser():
  import optparse

  usage   = 'usage: %prog [options] [args]'
  parser  = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-o', '--output', dest='output', metavar='NAME', type='string',
                    help='File output name.')
  parser.add_option('--graphtype', dest='graphtype', metavar='NAME', type='string', default='cumulative',
                    help='Type of graph desired for plotting data. Default = cumulative')
  parser.add_option('--xaxis', dest='xaxis', metavar='NAME', type='string', default=None,
                    help='A comma separated string of values to alter the xaxis dependent on the specified graph.')
  parser.add_option('--yaxis', dest='yaxis', metavar='NAME', type='string', default=None,
                    help='A comma separated string of values to alter the yaxis dependent on the specified graph.')

  parser.add_option('--graphstyle', dest='graphstyle', metavar='NAME', type='string', default='nci',
                    help='Style for the graph from the styles library. Default = nci')
  parser.add_option('--tabletype', dest='tabletype', metavar='NAME', type='string', default='cumulative',
                    help='Type of table desired for compiling the data. Default = cumulative')
  parser.add_option('--tablestyle', dest='tablestyle', metavar='NAME', type='string', default='nci',
                    help='Style for the graph from the styles library. Default = nci')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not args:
    print >> sys.stderr, 'Must provide input files'
    return

  if not options.output:
    print >> sys.stderr, 'Must provide output file name'
    return

  # Create models from input file.
  # The header section (meta) in file dictates which model to use.
  meta,model = parse_completion(args[0])

  # initialize pdf document
  doc   = DocTemplate(options.output + '.pdf')
  first = FirstPageGraphTemplate(id='First')
  later = LaterPageGraphTemplate(id='Later')

  # parameters to GraphPageStyle could come from tabular output and options
  info  = GraphPageStyle(meta, options)

  # Create all the graphs and tables for insertion into pdf
  content = create_content(model, first, options)

  doc.build(content, first, later, onFirstPage=info)


if __name__ == '__main__':
  main()
