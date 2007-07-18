# -*- coding: utf-8 -*-
'''
File:          reportlib.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-07-07

Requires:      Python 2.4, biozilla

Revision:      $Id: $
'''

from __future__ import with_statement

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

import  sys
import  time
import  csv
import  os

from    itertools                 import islice,chain,groupby,izip,repeat
from    operator                  import itemgetter
from    contextlib                import contextmanager

from    reportlab.platypus        import BaseDocTemplate,PageTemplate,Frame
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.rl_config       import defaultPageSize
from    reportlab.pdfgen          import canvas
from    reportlab.lib.units       import inch
from    reportlab.lib             import colors

from    biozilla.utils            import autofile, percent
from    biozilla.sections         import read_sections,filter_sections,index_sections
from    biozilla.stats            import quantile


from    styles                    import ncicolors
from    graphics                  import *


@contextmanager
def saveState(canvas):
  canvas.saveState()
  try:
    yield canvas
  finally:
    canvas.restoreState()

# TODO: use styles.PAGE_HEIGHT etc... throughout instead on import *
PAGE_HEIGHT=defaultPageSize[1]
PAGE_WIDTH=defaultPageSize[0]
SPACER=.125*inch
BORDER=SPACER*2
MISSING_INPUT = 'Analysis not performed. Requires data from investigator to complete.'


def _doNothing(canvas, doc):
  'Dummy callback for onPage'
  pass


class DocTemplate(BaseDocTemplate):
  '''
  Custom BaseDocTemplate generalized for the CGF reports.
  '''
  def handle_pageBegin(self):
    '''
    Override base method to add a change of page template after the firstpage.
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


class FirstPageGraphTemplate(PageTemplate):
  '''
  First page of graphs
  '''

  def __init__(self, id):
    self.pagewidth      = PAGE_WIDTH
    self.pageheight     = PAGE_HEIGHT
    self.headerheight   = SPACER*1.5
    self.xtop,self.ytop = BORDER,self.pageheight-BORDER

    titleheight     = SPACER*4
    subtitlewidth   = 2*inch
    x,y,w,h         = self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2)
    x,y,w,h         = x,y-titleheight,w-subtitlewidth,titleheight
    self.title      = Frame(x, y, w, h)
    self.subtitle   = Frame(x+w,y-SPACER,subtitlewidth,h+SPACER)

    x,y = BORDER+SPACER, BORDER+SPACER+self.headerheight
    w,h = self.pagewidth-(BORDER*2)-(SPACER*2), self.pageheight-(self.headerheight*2)-(BORDER*2)-(SPACER*2)-titleheight
    self.content = Frame(x,y,w,h)
    PageTemplate.__init__(self, id, self.content)

  def beforeDrawPage(self, canvas, doc):
    '''
    Before the page is drawn, generate the visual page layout based on the generated frames.
    '''
    with saveState(canvas):
      x,y,w,h = BORDER,BORDER,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2)
      canvas.setStrokeColor(colors.HexColor(ncicolors['black']))
      canvas.setFillColor(colors.HexColor(ncicolors['blue']))
      canvas.rect(x, y, w, h, stroke=1, fill=1)

      x,y,w,h = self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.headerheight
      drawbanner(canvas,x,y,w,h,bgcolor=ncicolors['red'],txtcolor=ncicolors['white'],left='National Cancer Institute',right='Core Genotyping Facility')

      x,y,w,h = self.content.x1,self.content.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.content.height
      canvas.setFillColor(colors.white)
      canvas.rect(x, y, w, h, fill=1, stroke=0)

      x,y,w,h = self.content.x1-SPACER,self.content.y1-self.headerheight-SPACER,self.pagewidth-(BORDER*2),self.headerheight
      drawbanner(canvas,x,y,w,h,bgcolor=ncicolors['red'],txtcolor=ncicolors['white'],right=str(canvas.getPageNumber()))


class LaterPageGraphTemplate(PageTemplate):
  '''
  All subsequent pages that follow the first page (FirstPageGraphTemplate) of the summary report.
  '''

  def __init__(self, id):
    self.pagewidth      = PAGE_WIDTH
    self.pageheight     = PAGE_HEIGHT
    self.headerheight   = SPACER*1.5
    self.xtop,self.ytop = BORDER,self.pageheight-BORDER

    x,y = BORDER+SPACER, BORDER+SPACER+self.headerheight
    w,h = self.pagewidth-(BORDER*2)-(SPACER*2), self.pageheight-(self.headerheight*2)-(BORDER*2)-(SPACER*2)
    self.content = Frame(x,y,w,h)
    PageTemplate.__init__(self, id, self.content)

  def beforeDrawPage(self, canvas, doc):
    '''
    Before the page is drawn, generate the visual page layout based on the generated frames.
    '''
    with saveState(canvas):
      x,y,w,h = BORDER,BORDER,self.pagewidth-(BORDER*2),self.pageheight-(BORDER*2)
      canvas.setStrokeColor(colors.HexColor(ncicolors['black']))
      canvas.setFillColor(colors.HexColor(ncicolors['blue']))
      canvas.rect(x, y, w, h, stroke=1, fill=1)

      x,y,w,h = self.xtop,self.ytop-self.headerheight,self.pagewidth-(BORDER*2),self.headerheight
      drawbanner(canvas,x,y,w,h,bgcolor=ncicolors['red'],txtcolor=ncicolors['white'],left='National Cancer Institute',right='Core Genotyping Facility')

      x,y,w,h = self.content.x1,self.content.y1,self.pagewidth-(BORDER*2)-(SPACER*2),self.content.height
      canvas.setFillColor(colors.HexColor(colors.white))
      canvas.rect(x, y, w, h, fill=1, stroke=0)

      x,y,w,h = self.content.x1-SPACER,self.content.y1-self.headerheight-SPACER,self.pagewidth-(BORDER*2),self.headerheight
      drawbanner(canvas,x,y,w,h,bgcolor=ncicolors['red'],txtcolor=ncicolors['white'],right=str(canvas.getPageNumber()))


def drawbanner(canvas,x,y,width,height,bgcolor=ncicolors['white'],txtcolor=ncicolors['black'],left='',center='',right=''):
  '''
  Creates a banner with the specified text
  '''
  with saveState(canvas):
    canvas.setFillColor(colors.HexColor(bgcolor))
    canvas.rect(x,y,width,height, stroke=0, fill=1)
    canvas.setFillColor(colors.HexColor(txtcolor))
    canvas.setFontSize(10)
    x,y   = x+SPACER, y+(SPACER/2)
    canvas.drawString(x,y,left)
    x     = x+width-BORDER
    canvas.drawRightString(x,y,right)
    x     = x/2
    canvas.drawCentredString(x,y,center)