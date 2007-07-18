# -*- coding: utf-8 -*-
'''
File:          summaryreport.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-07-07

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import  os
import  sys
import  csv
import  time

from    itertools            import islice,chain,groupby
from    operator             import itemgetter

from    reportlab.lib.styles import getSampleStyleSheet
from    reportlab.rl_config  import defaultPageSize

from    glu.lib.utils        import autofile, percent
from    glu.lib.sections     import read_sections,filter_sections,index_sections

from    reportlib            import *
from    styles               import *
from    graphics             import *


MISSING_INPUT = 'Analysis not performed. Requires data from investigator to complete.'

def summ_concordance(concordance, frame):
  '''
  Summarize concordance results.
  '''
  out = [p("Concordance", style=HeaderStyle3)]

  if concordance:
    raise NotImplementedError,'Summary report do not yet support concordance data'
  out.append( p(MISSING_INPUT) )

  return out


def summ_hwp(hwp, frame):
  '''
  Summarize hwp results.
  '''
  out = [p("Hardy–Weinberg Proportions", style=HeaderStyle3)]

  if hwp:
    raise NotImplementedError,'Summary report do not yet support Hardy-Weinberg summary'
  out.append( p(MISSING_INPUT) )

  return out


def summ_dupcheck(dupcheck, frame):
  '''
  Summarize dupcheck results.
  '''
  out = [p("Duplicate Checking", style=HeaderStyle3)]

  if dupcheck:
    raise NotImplementedError,'Summary report do not yet support dupcheck summary'
  else:
    out.append(p(MISSING_INPUT))

  return out


def summ_future():
  '''
  Place holder for future summarization of additional results.
  '''
  return [p("Contingency Association Analysis", style=HeaderStyle3),
          p(MISSING_INPUT)]


def generate_completion_table(samplecounts,sampletotal,locuscounts,locustotal,quantiles):
  '''
  Creates a styled table for the quantile counts and percents defined for completion by sample and locus.
  '''

  samplepercents  = [percent(n,sampletotal) for n in samplecounts]
  locuspercents   = [percent(n,locustotal)  for n in locuscounts]

  samplecount_quantile,samplepercent_quantile = compute_centiles(samplecounts,locustotal,quantiles)
  locuscount_quantile,locuspercent_quantile   = compute_centiles(locuscounts,sampletotal,quantiles)

  tabledata = [['','min','5th','med','95th','max',''],
              ['Samples by Locus']  +  ['%d'     % c for c in locuscount_quantile]        + ['#'],
              ['']                  +  ['%5.2f'  % per for per in locuspercent_quantile]  + ['%'],
              ['Loci by Sample']    +  ['%d'     % c for c in samplecount_quantile]       + ['#'],
              ['']                  +  ['%5.2f'  % per for per in samplepercent_quantile] + ['%']]
  return Table(tabledata, style=compcentileTS, hAlign='LEFT')


def generate_completion_graph(samplecounts,sampletotal,locuscounts,locustotal,quantiles):
  '''
  Creates a cumulative line plot for the completion counts of sample and locus
  '''
  quantiles = quantiles[1:4]
  xmax      = 1000
  labels    = [[q*xmax for q in quantiles],['5th','med','95th']]

  plotcentiles = [i/(xmax*1.00) for i in xrange(xmax+1)]

  samplecount_centile =   [(i,quantile(sorted(samplecounts), c, presorted=True)) for i,c in enumerate(plotcentiles)]
  locuscount_centile  =   [(i,quantile(sorted(locuscounts),  c, presorted=True)) for i,c in enumerate(plotcentiles)]

  samplepercent_centile = [(i,percent(c,locustotal))  for i,c in samplecount_centile]
  locuspercent_centile  = [(i,percent(c,sampletotal)) for i,c in locuscount_centile]

  graphcolors = [colors.CMYKColor(*ncilightred),colors.CMYKColor(*ncilightblue)]
  x,y,w,h     = 0,SPACER,190,90
  return lineplot([samplepercent_centile,locuspercent_centile],labels,graphcolors,x,y,w,h)


def summ_completion(completion,frame):
  '''
  Summarize completion results.
  '''
  out = [p('Completion', style=HeaderStyle3)]

  if completion:
    sample,locus      = extract_completion_results(completion)
    sam_summ,sam_data = sample
    loc_summ,loc_data = locus

    samples       = list(islice(sam_data,2,None))
    sampletotal   = len(samples)
    samplecounts  = [int(s[1]) for s in samples[:-1]]
    loci          = list(islice(loc_data,2,None))
    locustotal    = len(loci)
    locuscounts   = [int(l[1]) for l in loci[:-1]]

    quantiles = [0,.05,.5,.95,1]
    table     = generate_completion_table(samplecounts,sampletotal,locuscounts,locustotal,quantiles)
    graph     = generate_completion_graph(samplecounts,sampletotal,locuscounts,locustotal,quantiles)

    out.append(Table([[table,graph]], style=completionTS))
  else:
    out.append(p(MISSING_INPUT))

  return out

def analysisreport(page, completion, dupcheck, concordance, hwp):
  '''
  Drives the summary for each type of analysis and combines the resulting lists for
  inclusion into the reportlab's platypus machinery.

  To add additional summaries, create a summary method to be driven here and include
  the result list to the end of the return list below. Your summary method must create
  a list of reportlab flowables.
  '''
  analysisframe = page.analysis

  div           = [line(0,0,analysisframe.width-SPACER,0)]

  completion    = summ_completion(completion,analysisframe)
  concordance   = summ_concordance(concordance, analysisframe)
  hwp           = summ_hwp(hwp,analysisframe)
  dupcheck      = summ_dupcheck(dupcheck,analysisframe)
  future        = summ_future()

  return  completion  + div + \
          concordance + div + \
          hwp         + div + \
          dupcheck    + div + \
          future      #+ div + \
          # testframeflow()


def testframeflow():
  '''
  For testing the frame division and treatment of flowables.
  To be removed soon.
  '''
  test = [Spacer(1,0.2*inch)]
  for i in range(100):
    bogustext = ("This is Paragraph number %s.  " % i) *20
    txt = p(bogustext, ParaStyle)
    test.append(txt)
    test.append(Spacer(1,0.2*inch))
  return test

def summ_project(project):
  '''
  Format Project information such as its name, descrtiption, and
  possible comment made by a CGF project manager.
  '''
  return [p(project.get('name'),         HeaderStyle2),
         p(project.get('description'),  MedParaStyle),
         Spacer(.05*inch, .05*inch),
         p(project.get('comment'), SmallParaStyle)]


def summ_cgfproject(project):
  '''
  Format specific project information into a styled table.
  '''
  tabledata   = [['project:',      project.get('projectid'),   ''],
                 ['study:',        project.get('studyid'),     ''],
                 ['cost center:',  project.get('cas'),         ''],
                 ['investigator:', project.get('investigator'),''],
                 ['coordinator:',  project.get('cgfcontact'),  ''],
                 ['analyst:',      project.get('cgfanalyst'),  '']]

  tablestyle  = [('ALIGN',(0,0),(-2,-1),'LEFT'),
                 ('ALIGN',(1,0),(-2,-1),'RIGHT'),
                 ('SIZE',(0,0),(-1,-1),10),
                 ('TEXTCOLOR',(0,0),(-1,-1),darkncigrey)]

  return [Table (tabledata, style=tablestyle,rowHeights=12)]


def generate_sampledata(samples, summary, total):
  '''
  Generates a list of sample groupings for each sample status
  '''
  samples     = sorted(samples, key=itemgetter(8))
  groups      = groupby(samples, itemgetter(8))
  statues     = {}
  notassayed  = 0

  for group,members in groups:
    if group!='TRUE':
      num = len(list(members))
      notassayed += num
      statues[group] = num

  failedassay = len(summary.get('empty').split(','))+int(summary.get('dropped'))
  statues['failed genotyping'] = failedassay

  informative = total-failedassay-notassayed

  data = list([[num,status] for status,num in statues.iteritems()])
  data.append([informative,'informative samples'])

  return informative,data


def generate_assaydata(assays, summary, total):
  '''
  Generates a list of locus/assay groupings for each locus/assay status
  '''
  valid       = int(summary.get('total'))
  notvalid    = total-valid
  failed      = len(summary.get('empty').split(','))+int(summary.get('dropped'))
  informative = total-failed-notvalid

  data =[[failed,'loci failed in assay'],
          [informative,'informative loci']]

  if notvalid!=0:
    data.insert(0, [notvalid,'loci failed validation'])

  return informative,data

def summ_samples_assays(data, summary, frame, canvas, datatype):
  '''
  Summarizes the samples and assays for the project based on their status and
  genotyping results.
  '''
  summary = dict(summary)
  data    = list(data)
  total   = len(data)

  if datatype == 'sample':
    informative,data = generate_sampledata(data, summary, total)
  if datatype == 'assay':
    informative,data = generate_assaydata(data, summary, total)

  datalist=[p("%d / %d Informative Loci" % (informative,total), style=HeaderStyle3),
            Spacer(.05*inch, .05*inch)]

  ts = [('TEXTCOLOR',(0,e),(-1,-1),color) for e,color in enumerate(ncicolors)]

  datalist.append(Table(data, style=TableStyle(ts+summaryTS), rowHeights=9, hAlign='LEFT'))

  datagraph = [generate_pie(data,total,ncilightcolors,frame,canvas)]

  return datalist,datagraph


def summaryreport(canvas, doc):
  '''
  Driving method for summarizing project details including assay and sample summaries.
  This summary is independent of the analysis for the project with the exception of completion.
  '''
  projectframe      = doc.pageTemplates[0].project
  cgfprojectframe   = doc.pageTemplates[0].cgfproject
  sampleframe       = doc.pageTemplates[0].sample
  samplegraphframe  = doc.pageTemplates[0].samplegraph
  assayframe        = doc.pageTemplates[0].assay
  assaygraphframe   = doc.pageTemplates[0].assaygraph

  parser = option_parser()
  options,args = parser.parse_args()

  # project
  projectdef  = csv.reader(autofile(options.projectdef), dialect='excel-tab')
  sections    = index_sections(filter_sections(read_sections(projectdef),'project'))
  project     = dict(retrieve_section(sections, 'project'))
  projectframe.addFromList(summ_project(project), canvas)

  cgfprojectframe.addFromList(summ_cgfproject(project), canvas)

  # samples
  sampledef   = csv.reader(autofile(options.sampledef), dialect='excel-tab')
  sections    = index_sections(filter_sections(read_sections(sampledef), 'samples'))
  samples     = islice(retrieve_section(sections, 'samples'),1,None)

  # assays
  assaydef    = csv.reader(autofile(options.assaydef), dialect='excel-tab')
  sections    = index_sections(filter_sections(read_sections(assaydef),'assays'))
  assays      = islice(retrieve_section(sections, 'assays'),1,None)

  # completion
  sample,locus = extract_completion_results(options.completion)
  sam_summ,sam_data = sample
  loc_summ,loc_data = locus

  samplelist,samplegraph = summ_samples_assays(samples,sam_summ,sampleframe,canvas,'sample')
  sampleframe.addFromList(samplelist, canvas)
  samplegraphframe.addFromList(samplegraph, canvas)

  assaylist,assaygraph = summ_samples_assays(assays,loc_summ,assayframe,canvas,'assay')
  assayframe.addFromList(assaylist, canvas)
  assaygraphframe.addFromList(assaygraph, canvas)


def extract_completion_results(completion):
  '''
  Retrieves completion results from the file and prepares the appropriate data for use.
  '''
  COMPLETION_SECTIONS=['summary', 'data']
  completion  = list(csv.reader(autofile(completion), dialect='excel-tab'))
  completion  = index_sections(filter_sections(read_sections(completion),COMPLETION_SECTIONS))
  sam_summ    = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='sample')
  loc_data    = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='loci')
  loc_summ    = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='loci')
  sam_data    = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='sample')
  return (sam_summ,sam_data),(loc_summ,loc_data)

def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage)

  parser.add_option('-s', '--projectsummary', dest='projectsummary', metavar='N', type='str',
                    help='option not defined')

  parser.add_option('-c', '--completion', dest='completion', metavar='FILE',
                    help='option not defined')
  parser.add_option('-d', '--dupcheck', dest='dupcheck', metavar='FILE',
                    help='option not defined')
  parser.add_option('-n', '--concordance', dest='concordance', metavar='FILE',
                    help='option not defined')
  parser.add_option('-w', '--hwp', dest='hwp', metavar='FILE',
                    help='option not defined')

  parser.add_option('-p', '--projectdef', dest='projectdef', metavar='FILE',
                    help='option not defined')
  parser.add_option('-a', '--assaydef', dest='assaydef', metavar='FILE',
                    help='option not defined')
  parser.add_option('-m', '--sampledef', dest='sampledef', metavar='FILE',
                    help='option not defined')
  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  if not options.projectdef or not options.sampledef or not options.assaydef:
    print >> sys.stderr, 'Must provide definition files'
    return

  if options.projectsummary:
    doc = SummaryDocTemplate(options.projectsummary+'.pdf')

    first = FirstSummaryPage(id='First')
    later = LaterSummaryPage(id='Later')

    report = analysisreport(first, options.completion, options.dupcheck, options.concordance, options.hwp)

    doc.build(report, first, later, onFirstPage=summaryreport)


if __name__ == '__main__':
  main()