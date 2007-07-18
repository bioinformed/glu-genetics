# -*- coding: utf-8 -*-
'''
File:          graphics.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from    reportlab.graphics.shapes import Drawing
from    reportlab.graphics.shapes import Line,Rect
from    reportlab.lib             import colors

def lineplot(data,labels,graphcolor,x,y,w,h):
  from reportlab.graphics.charts.lineplots import LinePlot

  lp = LinePlot()

  lp.x      = x
  lp.y      = y
  lp.height = h*.95
  lp.width  = w*.95
  lp.data   = data
  lp.lines[0].strokeColor = graphcolor[0]
  lp.lines[1].strokeColor = graphcolor[1]

  lp.xValueAxis.visibleGrid       = 1
  lp.xValueAxis.valueSteps        = labels[0]
  lp.xValueAxis.labelTextFormat   = labels[1]
  lp.xValueAxis.labels.fontSize   = 8
  lp.xValueAxis.labels.textAnchor = 'end'

  lp.yValueAxis.valueMin  = 80
  lp.yValueAxis.valueMax  = 100
  lp.yValueAxis.valueStep = 10

  drawing = Drawing(w, h, hAlign='RIGHT')
  drawing.add(lp)

  return drawing


def drawingat(drawing,x,y,canvas):
  from reportlab.graphics import renderPDF
  renderPDF.draw(drawing, canvas, x,y)

def pie(data,labels,piecolors,x,y,diameter):
  from reportlab.graphics.charts.piecharts import Pie

  pc = Pie()

  pc.startAngle         = 240
  pc.x                  = x
  pc.y                  = y
  pc.width              = diameter
  pc.height             = diameter
  pc.data               = data
  pc.labels             = labels
  pc.slices.strokeWidth = 0.5
  pc.slices.labelRadius = 1.4
  pc.slices.fontSize    = 8

  for i,c in enumerate(piecolors):
    pc.slices[i].fillColor = c
    pc.slices[i].fontColor = c

  d = Drawing(diameter,diameter)
  d.add(pc)

  return d

def line(x1,y1,x2,y2):
  height  = 15
  d       = Drawing(1,height, hAlign='CENTER')
  l       = Line(x1,y1+(height/2),x2,y2+(height/2), strokeColor=colors.black, strokeWidth=.5)
  d.add(l)
  return d
