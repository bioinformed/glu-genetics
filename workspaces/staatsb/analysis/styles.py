from    reportlab.lib.styles      import ParagraphStyle
from    reportlab.lib             import colors
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.lib.enums       import TA_RIGHT
from    reportlab.platypus        import TableStyle,Paragraph

stylesheet = getSampleStyleSheet()

# COLORS
basecolors  = { 'white': '#FFFFFF',
                'black': '#000000'}


ncicolors   = { 'red':      '#9C2C20',
                'darkred':  '#B72622',
                'blue':     '#6A80A3',
                'darkblue': '#2741B2'
              }
ncicolors.update(basecolors)

# My way <---
#basecolors.update(ncicolors)

# DELTE SOON
# ncipurple       = .45,.36,.10,.01
# nciblue         = .49,.32,.08,.01
# ncidarkblue     = '#336699' #'#669999'
# ncigreen        = .52,.25,.29,.01
# ncitan          = .49,.43,.49,.04
# ncigrey         = .41,.34,.31,0
# ncired          = .06,.85,.98,.18
# ncidarkred      = '#990000'
# ncicolors       = [colors.CMYKColor(*ncired),colors.CMYKColor(*ncigrey),colors.CMYKColor(*ncitan),colors.CMYKColor(*ncigreen),colors.CMYKColor(*nciblue),colors.CMYKColor(*ncipurple)]
#
# ncilightpurple  = .30,.22,.07,.0
# ncilightblue    = .43,.29,.09,0
# ncilightgreen   = .33,.1,.15,.0
# ncilighttan     = .23,.19,.26,.0
# ncilightgrey    = .26,.20,.16,.0
# ncilightred     = .24,.41,.29,.01
# ncilightcolors  = [colors.CMYKColor(*ncilightred),colors.CMYKColor(*ncilightgrey),colors.CMYKColor(*ncilighttan),colors.CMYKColor(*ncilightgreen),colors.CMYKColor(*ncilightblue),colors.CMYKColor(*ncilightpurple)]
#
# # ncidarkblue     = .70,.48,.15,.04
# darkncigrey     = .60,.51,.47,.04

class Style(object):
  '''The Style class.'''
  def __init__(self, strokewidth=0, strokecolor=None, fillcolor=None, fontsize=8, fontcolor=None, fontname=None):
    self.strokewidth  = strokewidth
    self.strokecolor  = strokecolor
    self.fillcolor    = fillcolor
    self.fontsize     = fontsize
    self.fontcolor    = fontcolor
    self.fontname     = fontname

  def paragraphstyle(self, name, **kwargs):
    return ParagraphStyle(name, parent=stylesheet['Normal'],
                          fontSize=self.fontsize,
                          fontname=self.fontname,
                          textColor=self.fontcolor, **kwargs)


class GraphStyle(object):
  def __init__(self, data, xaxis, yaxis, title):
    self.title  = title
    self.data   = data
    self.xaxis  = xaxis
    self.yaxis  = yaxis


def p(txt,style=stylesheet['Normal']):
  return Paragraph(txt, style=style)
