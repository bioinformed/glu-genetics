from    reportlab.graphics.shapes               import Drawing,Line,Rect,Circle,String
from    reportlab.lib                           import colors
from    reportlab.graphics.charts.lineplots     import LinePlot
from    reportlab.graphics.widgets.markers      import makeMarker
from    reportlab.graphics.charts.textlabels    import Label


class Axis(object):
  def __init__(self, amin, amax, steps, name=None, labels=None):
    self.name             = name
    self.amin             = amin
    self.amax             = amax
    self.steps            = steps
    self.range            = range(amin, amax+1, steps)  # <--- wrong, since range is defined only over ints
                          # steps,text
    self.labels           = labels if labels else [self.range,map(str, self.range)] # <-- map(str,...) gives you no control
    self.labels_maxlength = max(len(l) for l in self.labels[1])


def create_graph_axis(axis, axismodel, axisstyle, visible=False):
  axis.valueMin        = axismodel.amin
  axis.valueMax        = axismodel.amax
  axis.valueStep       = axismodel.steps
  axis.labels.fontSize = axisstyle.fontsize
  axis.valueSteps,axis.labelTextFormat = axismodel.labels
  axis.visibleGrid     = bool(visible)


def create_label(text, fontsize, x, y):
  lab             = Label()
  lab.textAnchor  = 'middle'
  lab.fontSize    = fontsize
  lab.setOrigin(x,y)
  lab.setText(text)
  return lab


def draw_lineplot(data, xaxis, yaxis, style, width, height, table=None):
  '''
  Creates a line plot within a Drawing object.

    NOTE: The width and height are the inner width/height of the graph,
          need to reduce size to fit labels and tick markers

          For the axis labels, reportlab does not give a way to calculate label length
          in order to compensate for the space needed for the label. Thus a temporary
          approximation is used: length of maximum label on axis + 1/4 of the fontsize.
          Works ok with for most commonly used variations of fontsize and label length.
  '''
  from styles import Style

  drawing   = Drawing(width, height, hAlign='CENTER')
  lp        = LinePlot()

  lp.width  = width-(style.yaxis.fontsize*5)    # 3 fontsizes for label compensation and padding
  lp.height = height-(style.xaxis.fontsize*5)   # 3 fontsizes for label compensation and padding
  lp.x      = lp.x+xaxis.labels_maxlength+(style.xaxis.fontsize/4)
  lp.y      = lp.y+style.yaxis.fontsize
  lp.data   = data

  drawing.add(lp)

  # yaxis settings
  create_graph_axis(lp.yValueAxis, yaxis, style.yaxis)

  if yaxis.name:
    label = create_label(yaxis.name, style.yaxis.fontsize+3, 0, lp.height/2)
    label.angle = 90
    drawing.add(label)

  # xaxis settings
  create_graph_axis(lp.xValueAxis, xaxis, style.xaxis, visible=True)

  if xaxis.name:
    label = create_label(xaxis.name, style.xaxis.fontsize+3, lp.width/2, 0)
    drawing.add(label)

  if table:
    draw_graph_table(table, style, lp, drawing)

  # color code plot lines
  for i,l in enumerate(style.data):
    lp.lines[i].strokeColor = colors.HexColor(style.data[i].strokecolor)

  return drawing


def draw_graph_table(table, style, graph, drawing):
  steplength  = graph.width/(len(table)-1) # <--- worried about the -1
  y           = -25 # padding between xaxis and table

  for i,d in enumerate(style.data):
    x  = graph.x
    for l in table:
      drawing.add(String(x,y, l[i], fontSize=style.xaxis.fontsize+2,
                                    fillColor=colors.HexColor(style.data[i].strokecolor),
                                    textAnchor='middle'))
      x += steplength

    y -= 2*style.xaxis.fontsize # (style.xaxis.fontsize*2) is padding between lines



def lineplot(data,labels,graphcolor,x,y,w,h):
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
