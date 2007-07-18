# -*- coding: utf-8 -*-
'''
File:          reportlib.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-07-17

Abstract:      Supporting report functions

Requires:      Python 2.4, biozilla

Revision:      $Id: $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import  os
import  sys
import  time

from    reportlab.platypus        import BaseDocTemplate,SimpleDocTemplate,Paragraph,Frame,Image,Table,Spacer,PageTemplate
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.rl_config       import defaultPageSize
from    reportlab.lib.units       import inch
from    reportlab.graphics.shapes import Line,Rect
from    reportlab.pdfgen          import canvas

from    biozilla.utils            import percent
from    biozilla.stats            import quantile

from    graphics                  import *
from    styles                    import *


IMAGEDIR = '/home/staatsb/images/'

PAGE_HEIGHT=defaultPageSize[1]
PAGE_WIDTH=defaultPageSize[0]
SPACER=.125*inch
BORDER=SPACER*2


def compute_centiles(data, total, centiles):
  '''
  Utility method for summary reports
  '''
  quantiles   = [quantile(sorted(data),  c, presorted=True) for c in centiles]
  percentiles = [percent(q,total)  for q in quantiles]
  return quantiles,percentiles


def generate_pie(data,total,colors,frame,canvas):
  '''
  Utility method for summary reports to hide some of the implimentation
  details of creating a pie graph
  '''
  data      = [num for num,status in data]
  labels    = ['%2.1f' % percent(num,total) for num in data]
  pieheight = frame.height*.65
  x,y       = SPACER*1.8,SPACER*-1.65
  return pie(data,labels,colors,x,y,pieheight)


def retrieve_section(sections, heading, section_type=None):
  '''
  Uses the section.py library to retireve specified sections. Handles multiple
  sections with with the same heading and performs some error checks.
  '''
  section = sections.get(heading)

  if not section_type and len(section) >1:
    raise NoneUniqueSectionHeadingError, 'More than one section with the heading: ' + heading

  if not section_type:
    return section[0]

  for items in section:
    for item in items:
      if section_type in set(item):
        return items


def _doNothing(canvas, doc):
  "Dummy callback for onPage"
  pass

def make_ncibanner(x,y,width,height,canvas):
  '''
  Creates the summary report's header
  '''
  canvas.saveState()
  canvas.setFillColorCMYK(*ncired)
  canvas.rect(x,y,width,height, stroke=0, fill=1)
  canvas.setFillColor(colors.white)
  canvas.setFontSize(10)
  x,y = x+SPACER,y+(SPACER/2)
  canvas.drawString(x,y,'National Cancer Institute')
  canvas.drawRightString(x+width-BORDER,y,'Core Genotyping Facility')
  canvas.restoreState()


def make_ncifooter(x,y,width,height,canvas):
  '''
  Creates the summary report's footer
  '''
  canvas.saveState()
  canvas.setFillColorCMYK(*ncired)
  canvas.rect(x,y,width,height, stroke=0, fill=1)
  canvas.setFillColor(colors.white)
  canvas.setFontSize(10)
  x,y = x+SPACER,y+(SPACER/2)
  canvas.drawString(x,y,'CGF Project Summary')
  canvas.drawRightString(x+width-BORDER,y,str(canvas.getPageNumber()))
  canvas.restoreState()


def buildlogos(x,y,width, canvas):
  '''
  Builds the standard affiliated government logos for the CGF

  The logo images are stored in IMAGEDIR, which also may contain other images
  '''
  canvas.saveState()

  # TODO: generate higher resolution logos (partially done)
  footernci       = os.path.join(os.path.dirname(IMAGEDIR),'footer_ncis.gif')
  footerhhs       = os.path.join(os.path.dirname(IMAGEDIR),'footer_hhss.gif')
  footernih       = os.path.join(os.path.dirname(IMAGEDIR),'footer_nihs.gif')
  footerfirstgov  = os.path.join(os.path.dirname(IMAGEDIR),'footer_firstgovs.gif')

  ncigif      = Image(footernci)
  hhsgif      = Image(footerhhs)
  nihgif      = Image(footernih)
  firstgovgif = Image(footerfirstgov)
  gifwidth    = ncigif.drawWidth + hhsgif.drawWidth + nihgif.drawWidth + firstgovgif.drawWidth

  x2,y2 = width-(SPACER*4)-gifwidth,y+SPACER
  canvas.drawImage(footernci, x2, y2)

  x2 = x2+ncigif.drawWidth
  canvas.drawImage(footerhhs, x2, y2)

  x2 = x2+hhsgif.drawWidth
  canvas.drawImage(footernih, x2, y2)

  x2 = x2+nihgif.drawWidth
  canvas.drawImage(footerfirstgov, x2, y2)

  canvas.restoreState()


def backdrop(x,y,width,height,canvas,fillcolor=None,strokecolor=None):
  '''
  Generates rectangles with specified fill and stroke colors.
  '''
  canvas.saveState()

  fill = stroke = 0
  if strokecolor:
    canvas.setStrokeColor(strokecolor)
    stroke=1

  if fillcolor:
    canvas.setFillColor(fillcolor)
    fill=1

  canvas.rect(x,y,width,height, stroke=stroke, fill=fill)
  canvas.restoreState()


class LaterSummaryPage(PageTemplate) :
  '''
  All subsequent pages that follow the first page (FirstSummaryPage) of the summary report.
  '''

  def __init__(self, id):
    self.pagewidth      = PAGE_WIDTH
    self.pageheight     = PAGE_HEIGHT
    self.headerheight   = SPACER*1.5
    self.xtop,self.ytop = BORDER,self.pageheight-BORDER

    x,y,w,h = BORDER+SPACER, BORDER+SPACER+self.headerheight, self.pagewidth-(BORDER*2)-(SPACER*2), self.pageheight-(self.headerheight*2)-(BORDER*2)-(SPACER*2)
    self.analysis = Frame(x,y,w,h)
    PageTemplate.__init__(self, id, self.analysis)

  def beforeDrawPage(self, canvas, doc):
    '''
    Before the page is drawn, generate the visual page layout based on the generated frames.
    '''
    x,y,w,h = BORDER,BORDER,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2)
    backdrop(x,y,w,h,canvas,fillcolor=nciblue,strokecolor=colors.black)

    x,y,w,h = self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.headerheight
    make_ncibanner(x,y,w,h,canvas)

    x,y,w,h = self.analysis.x1,self.analysis.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.analysis.height
    backdrop(x,y,w,h,canvas,fillcolor=colors.white)

    x,y,w,h = self.analysis.x1-SPACER,self.analysis.y1-self.headerheight-SPACER,self.pagewidth-(BORDER*2),self.headerheight
    make_ncifooter(x,y,w,h,canvas)


class FirstSummaryPage(PageTemplate):
  '''
  This is the first page of the summary report.
  '''
  def __init__(self, id):
    self.pagewidth      = PAGE_WIDTH
    self.pageheight     = PAGE_HEIGHT
    self.headerheight   = SPACER*1.5
    self.xtop,self.ytop = BORDER,self.pageheight-BORDER

    self.buildframes(self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2))

    PageTemplate.__init__(self, id, self.analysis)


  def beforeDrawPage(self, canvas, doc):
    '''
    Before the page is drawn, generate the visual page layout based on the generated frames.
    '''
    x,y,w,h = BORDER,BORDER,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2)
    backdrop(x,y,w,h,canvas,fillcolor=nciblue,strokecolor=colors.black)

    x,y,w,h = self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.headerheight
    make_ncibanner(x,y,w,h,canvas)

    x,y,w,h = self.project.x1,self.project.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.project.height
    backdrop(x,y,w,h,canvas,fillcolor=colors.white)

    x,y,w,h = self.sample.x1,self.sample.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.sample.height
    backdrop(x,y,w,h,canvas,fillcolor=colors.white)

    x,y,w,h = self.analysis.x1,self.analysis.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.analysis.height
    backdrop(x,y,w,h,canvas,fillcolor=colors.white)

    buildlogos(x,y,self.pagewidth,canvas)

    x,y,w,h = self.analysis.x1-SPACER,self.analysis.y1-self.headerheight-SPACER,self.pagewidth-(BORDER*2),self.headerheight
    make_ncifooter(x,y,w,h,canvas)

    self.title.add(p('CGF Project Summary', style=titleStyle), canvas)
    self.subtitle.addFromList([ p('CGF Reporting v0.2(BETA)', style=subtitleStyle),p(time.asctime(), style=subtitleStyle)], canvas)


  def buildframes(self, x,y,width,height):
    '''
    Generates the frames for this page.
    '''
    titleheight     = SPACER*4
    subtitlewidth   = 2*inch
    x,y,w,h         = x,y-titleheight,width-subtitlewidth,titleheight
    self.title      = Frame(x, y, w, h)
    self.subtitle   = Frame(x+w,y-SPACER,subtitlewidth,h+SPACER)

    projectheight = SPACER*10
    projectwidth  = width*.65
    x,y,w,h       = x+SPACER,y-SPACER-projectheight,projectwidth,projectheight
    self.project  = Frame(x,y,w,h)

    managerwidth    = width*.35
    x1,y1,w1,h1     = x+w, y, managerwidth-(SPACER*2), h
    self.cgfproject = Frame(x1,y1,w1,h1)

    samplewidth   = (width-(SPACER*2))/4
    sampleheight  = SPACER*10
    x,y,w,h       = x,y-SPACER-sampleheight,samplewidth+(BORDER*2),sampleheight
    self.sample   = Frame(x,y,w,h)

    x1,y1,w1,h1       = x+samplewidth+(BORDER*2),y,samplewidth-(BORDER*2),h
    self.samplegraph  = Frame(x1,y1,w1,h1)
    x1,y1,w1,h1       = x1+samplewidth-(BORDER*2),y,samplewidth+(BORDER*2),h
    self.assay        = Frame(x1,y1,w1,h1)
    x1,y1,w1,h1       = x1+samplewidth+(BORDER*2),y,samplewidth-(BORDER*2),h
    self.assaygraph   = Frame(x1,y1,w1,h1)

    self.analysisheight = y-(SPACER*4)-self.headerheight
    x,y,w,h             = x,y-SPACER-self.analysisheight,width-(SPACER*2),self.analysisheight
    self.analysis       = Frame(x,y,w,h)



class SummaryDocTemplate(BaseDocTemplate):
  '''
  Custom summary DocTemplate speciicly for the CGF summary reports.
  '''
  def handle_pageBegin(self):
    '''override base method to add a change of page template after the firstpage.
    '''
    self._handle_pageBegin()
    self._handle_nextPageTemplate('Later')

  def build(self,flowables, firstPage, laterPage, onFirstPage=_doNothing, onLaterPages=_doNothing, canvasmaker=canvas.Canvas):
    self._calc()    #in case we changed margins sizes etc

    self.addPageTemplates([firstPage,laterPage])

    if onFirstPage is not _doNothing:
      self.pageTemplates[0].afterDrawPage = onFirstPage
    if onLaterPages is not _doNothing:
      self.pageTemplates[1].afterDrawPage = onLaterPages

    BaseDocTemplate.build(self,flowables, canvasmaker=canvasmaker)


def main():
  pass

if __name__ == '__main__':
  main()