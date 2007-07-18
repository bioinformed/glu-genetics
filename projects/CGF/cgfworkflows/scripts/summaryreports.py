# -*- coding: utf-8 -*-
'''
File:          summaryreport.py

Authors:       Brian Staats (staatsb@mail.nih.gov), Nick Xiao (xiaon@mail.nih.gov)

Created:       2006-07-07

Requires:      Python 2.4, biozilla

Revision:      $Id: $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2006 Science Applications International Corporation ("SAIC"). All rights reserved.'

Title = "Competion Report"

Author = "Brian Staats"

URL = "http://cgf.nci.nih.gov"

email = "staatsb@mail.nih.gov"

Abstract = """Completion analysis explaination."""


import  sys
import  time
import  csv
import  os

from    reportlab.platypus        import *
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.lib.styles      import ParagraphStyle
from    reportlab.rl_config       import defaultPageSize
from    reportlab.lib.units       import inch
from    reportlab.lib             import colors
from    reportlab.graphics.shapes import Drawing,Line,Rect
from    reportlab.lib.enums       import *

from    biozilla.utils            import autofile, percent
from    biozilla.stats            import quantile
from    biozilla.sections         import read_sections,filter_sections,index_sections
from    itertools                 import islice,chain,groupby
from    operator                  import itemgetter
from    reports                   import retrieve_section


imagedir = '/home/staatsb/images/'

MISSING_INPUT = 'Analysis not performed. Requires data from investigator to complete.'

# colors colors colors
# golden colors
# saturation: light, normal, heavy
# color order red, green, blue
goldencolors = [[ (.08,.13,.1,0),(.11,.05,.15,0),(.17,.13,.03,0) ],
                [ (.01,.31,.35,0),(.17,0,.4,0),(.45,.33,0,0) ],
                [ (.01,.66,.88,.02),(.28,0,.79,0),(.82,.61,0,0) ]]


ncipurple = .45,.36,.10,.01
nciblue = .49,.32,.08,.01
ncigreen = .52,.25,.29,.01
ncitan = .49,.43,.49,.04
ncigrey = .41,.34,.31,0
ncired = .06,.85,.98,.18
ncicolors = [colors.CMYKColor(*ncipurple),colors.CMYKColor(*nciblue),colors.CMYKColor(*ncigreen),colors.CMYKColor(*ncitan),colors.CMYKColor(*ncigrey),colors.CMYKColor(*ncired)]
ncicolors.reverse()

ncilightpurple = .30,.22,.07,.0
ncilightblue = .43,.29,.09,0
ncilightgreen = .33,.1,.15,.0
ncilighttan = .23,.19,.26,.0
ncilightgrey = .26,.20,.16,.0
ncilightred = .24,.41,.29,.01
ncilightcolors = [colors.CMYKColor(*ncilightpurple),colors.CMYKColor(*ncilightblue),colors.CMYKColor(*ncilightgreen),colors.CMYKColor(*ncilighttan),colors.CMYKColor(*ncilightgrey),colors.CMYKColor(*ncilightred)]
ncilightcolors.reverse()


ncidarkblue = .70,.48,.15,.04
darkncigrey = .60,.51,.47,.04



PAGEHEIGHT=defaultPageSize[1]
PAGEWIDTH=defaultPageSize[0]
styles = getSampleStyleSheet()

pageinfo = "%s / %s / %s" % (Author, email, Title)
HeaderStyle1 = styles["Heading1"] # XXXX
HeaderStyle2 = styles["Heading2"]
HeaderStyle3 = styles["Heading3"]
HeaderStyle3.textColor=darkncigrey
ParaStyle = styles["Normal"]
PreStyle = styles["Code"]

SmallParaStyle = ParagraphStyle('SmallParaStyle',parent=ParaStyle,
                            fontSize=8, fontName='Times-Roman', leading=10,
                            textColor=colors.CMYKColor(*darkncigrey))

MedParaStyle = ParagraphStyle('MedParaStyle',parent=ParaStyle,
                            fontSize=8, fontName='Times-Roman', leading=10,
                            textColor=colors.CMYKColor(*darkncigrey))

titleStyle = ParagraphStyle('titleStyle',parent=HeaderStyle1,
                            fontSize=24, fontName='Helvetica-Bold',
                            textColor=colors.white)

subtitleStyle = ParagraphStyle('subtitleStyle',parent=ParaStyle,
                            fontSize=10, fontName='Helvetica-Bold', alignment=TA_RIGHT,
                            textColor=colors.white)

nciStyle = ParagraphStyle('nciStyle',parent=ParaStyle,
                            fontSize=10, fontName='Helvetica',
                            textColor=colors.white)


detailStyle = ParagraphStyle('detailStyle',parent=ParaStyle,fontSize=8)
defStyle = ParagraphStyle('defStyle',parent=ParaStyle,leftIndent=8)

# ratio = 1.618033989

def header(txt, style=HeaderStyle1, klass=Paragraph, sep=0.3):
    # s = Spacer(0.2*inch, sep*inch)
    # Elements.append(s)
    para = klass(txt, style)
    return para


def p(txt,style=ParaStyle):
    return Paragraph(txt, style=style)


def sub(txt,style=PreStyle):
  return  Paragraph(txt, style)


def pre(txt):
    s = Spacer(0.1*inch, 0.1*inch)
    Elements.append(s)
    p = Preformatted(txt, PreStyle)
    Elements.append(p)


headerTS = TableStyle([('ALIGN',(0,0),(1,1),'LEFT'),
                       ('ALIGN',(-1,-1),(-2,-2),'RIGHT'),
                          ('SIZE',(0,0),(0,0),18),
                          ('FACE',(0,0),(0,0),'Times-Bold'),
                          ('FACE',(-1,-1),(-2,-2),'Helvetica-Bold'),
                          # ('BACKGROUND', (0,0), (-1,-2), goldencolors[0][2]),
                          # ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                          ])


def make_ncibanner(x,y,width,height,SPACER,canvas):
  canvas.saveState()
  canvas.setFillColorCMYK(*ncired)
  canvas.rect(x*inch,y*inch,width*inch,height*inch, stroke=0, fill=1)
  canvas.setFillColor(colors.white)
  canvas.setFontSize(10)
  canvas.drawString((x+SPACER)*inch,(y+(SPACER/2))*inch,'National Cancer Institute')
  # canvas.drawRightString((x+width-(SPACER*3))*inch,(y+(SPACER/2))*inch,'U.S. National Institutes of Health | www.cancer.gov')
  canvas.drawRightString((x+width-SPACER)*inch,(y+(SPACER/2))*inch,'Core Genotyping Facility')
  canvas.restoreState()

def make_ncifooter(x,y,width,height,SPACER,canvas):
  canvas.saveState()
  canvas.setFillColorCMYK(*ncired)
  canvas.rect(x*inch,y*inch,width*inch,height*inch, stroke=0, fill=1)
  canvas.setFillColor(colors.white)
  # canvas.drawString((x+SPACER)*inch,(y+(SPACER/2))*inch,'Core Genotyping Facility')
  # canvas.drawRightString((x+width-(SPACER*3))*inch,(y+(SPACER/2))*inch,'Division of Cancer Epidemiology and Genetics')
  canvas.restoreState()

def backdrop(x,y,w,h,color,canvas):
  canvas.saveState()
  canvas.setFillColorCMYK(*color)
  canvas.setStrokeColor(colors.black)
  canvas.rect(x*inch,y*inch,w*inch,h*inch, stroke=0, fill=1)
  canvas.restoreState()

def frontdrop(x,y,w,h,color,canvas):
  canvas.saveState()
  canvas.setStrokeColor(color)
  canvas.rect(x*inch,y*inch,w*inch,h*inch, stroke=1, fill=0)
  canvas.restoreState()

SPACER=.125
BORDER=SPACER*2
class SummaryPage(object):
  def __init__(self, canvas, width, height):
    self.canvas = canvas
    self.canvas.saveState()
    self.buildpage(width, height)
    self.canvas.restoreState()

  def section(self, x,y,width,height, bgcolor=None):
    self.canvas.saveState()
    if bgcolor:
      self.canvas.setFillColor(bgcolor)
      self.canvas.rect(x*inch, y*inch, width*inch, height*inch, stroke=0, fill=1)

    f = Frame(x*inch, y*inch, width*inch, height*inch, showBoundary=0)
    self.canvas.restoreState()
    return f

  def buildpage(self,width, height):
    self.canvas.saveState()
    width,height = width/inch,height/inch

    xb,yb,wb,hb = BORDER,BORDER,width-(BORDER*2),height-(BORDER*2)
    backdrop(xb,yb,wb,hb,nciblue,self.canvas)

    nciheight = SPACER*1.5
    x,y,w,h = xb,height-yb-nciheight,wb,nciheight
    make_ncibanner(x,y,w,h,SPACER,self.canvas)

    titleheight = SPACER*4
    headwidth = 2
    x,y,w,h = x,y-titleheight,w-headwidth,titleheight
    self.title = Frame(x*inch, y*inch, w*inch, h*inch, showBoundary=0)
    self.title.addFromList([ p('CGF Project Summary', style=titleStyle)], self.canvas)

    import time
    self.header = Frame((x+w)*inch, (y-SPACER)*inch, headwidth*inch, (h+SPACER)*inch, showBoundary=0)
    self.header.addFromList([ p('CGF Reporting v0.2(BETA)', style=subtitleStyle),
                              p(time.asctime(), style=subtitleStyle)], self.canvas)

    projectheight,projectwidth = SPACER*10,width*.75
    x,y,w,h = x+SPACER,y-SPACER-projectheight,projectwidth-(BORDER*2),projectheight

    self.canvas.setFillColor(colors.white)
    self.canvas.rect(x*inch,y*inch,(width-(BORDER*2)-(SPACER*2))*inch,h*inch, stroke=0, fill=1)
    self.canvas.setFillColor(colors.black)

    self.project = self.section(x,y,w,h)

    managerheight,managerwidth = SPACER*10,width*.25
    x1,y1,w1,h1 = x+w, y, managerwidth-(SPACER*2), managerheight
    self.cgfproject = self.section(x1,y1,w1,h1)

    sampleheight = SPACER*9
    x,y,w,h = x,y-SPACER-sampleheight,width-(BORDER*4),sampleheight

    self.canvas.setFillColor(colors.white)
    self.canvas.rect(x*inch,y*inch,(width-(BORDER*2)-(SPACER*2))*inch,h*inch, stroke=0, fill=1)
    self.canvas.setFillColor(colors.black)

    self.sample = self.section(x,y,w/2,h)

    assayheight = SPACER*9
    x1,y1,w1,h1 = x+(w/2),y,(width-(SPACER*4))/2,assayheight
    self.assay = self.section(x1,y1,w1,h1)

    analysisheight = y-(SPACER*4)-nciheight
    x,y,w,h = x,y-SPACER-analysisheight,width-(BORDER*2)-(SPACER*2),analysisheight
    self.analysis = self.section(x,y,w,h,colors.white)

    footernci = os.path.join(os.path.dirname(imagedir),'footer_ncis.gif')
    footerhhs = os.path.join(os.path.dirname(imagedir),'footer_hhss.gif')
    footernih = os.path.join(os.path.dirname(imagedir),'footer_nihs.gif')
    footerfirstgov = os.path.join(os.path.dirname(imagedir),'footer_firstgovs.gif')

    ncigif = Image(footernci)
    hhsgif = Image(footerhhs)
    nihgif = Image(footernih)
    firstgovgif = Image(footerfirstgov)
    gifwidth = ncigif.drawWidth+hhsgif.drawWidth+nihgif.drawWidth+firstgovgif.drawWidth

    x2,y2 = ((width-(SPACER*4))*inch)-gifwidth,(y+SPACER)*inch
    self.canvas.drawImage(footernci, x2,y2)

    x2 = x2+ncigif.drawWidth
    self.canvas.drawImage(footerhhs, x2,y2)

    x2 = x2+hhsgif.drawWidth
    self.canvas.drawImage(footernih, x2,y2)

    x2 = x2+nihgif.drawWidth
    self.canvas.drawImage(footerfirstgov, x2,y2)

    make_ncifooter(x-SPACER,y-nciheight-SPACER,w+(SPACER*2),nciheight,SPACER,self.canvas)
    frontdrop(xb,yb,wb,hb,colors.black,self.canvas)

    self.canvas.restoreState()


class GoldenSummaryPage(SummaryPage):
  def buildpage(self,width, height):
    from reporttemplates import GoldenRect

    self.canvas.saveState()
    width,height = width/inch,height/inch

    # backdrop
    c,m,y,b = goldencolors[0][2]
    self.canvas.setFillColorCMYK(c,m,y,b)
    x,y,w,h = .25*inch, .25*inch, (width-.5)*inch, (height-.5)*inch
    self.canvas.rect(x,y,w,h, stroke=0, fill=1)
    self.canvas.setFillColor(colors.black)

    gwidth = ((height-.25)/GoldenRect.RATIO)*.95
    gr0 = GoldenRect.grect(gwidth)
    gr0.x=gr0.y= SPACER*2

    gr1,gs1 = gr0.divide('up')
    gr1.x,gr1.y=gr1.x+SPACER,gr1.y+SPACER
    gs1.x,gs1.y=gs1.x+SPACER,gs1.y+SPACER
    self.summaries = self.section(gs1.x,gs1.y,gs1.width,gs1.width,1)

    graphsize=1.9
    x,y,w,h = gs1.width+gs1.x-1,gs1.y+gs1.width-graphsize-SPACER,graphsize*1.2,graphsize
    self.graph1 = self.section(x,y,w,h)

    x,y,w,h = x,y-graphsize-SPACER,w,h
    self.graph2 = self.section(x,y,w,h)

    x,y,w,h = x,y-graphsize-SPACER,w,h
    self.graph3 = self.section(x,y,w,h)

    gr1.width=gr1.width-SPACER
    gr1.height=gr1.height-SPACER
    gr2,gs2 = gr1.divide('right')
    gr2.x,gr2.y=gr2.x+SPACER,gr2.y+SPACER
    gs2.y=gs2.y+SPACER
    self.project = self.section(gs2.x,gs2.y,gs2.width,gs2.width)

    gr2.width=gr2.width+SPACER
    gr2.height=gr2.height+SPACER
    gr3,gs3 = gr2.divide('down')
    self.graphlegend = self.section(gr3.x,gr3.y,gr3.width-SPACER,gr3.height-SPACER)
    self.assay = self.section(gs3.x,gs3.y,gs3.width-SPACER,gs3.width-SPACER)
    self.canvas.restoreState()

completionTS = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                          ('VALIGN',(0,0),(-1,-1),'TOP'),
                          ])

compcompTS = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                                    # ('ALIGN',(3,0),(-4,-1),'RIGHT'),
                                    # ('ALIGN',(5,0),(-3,-1),'LEFT'),
                                    ('LINEABOVE',(1,1),(-1,-2),1,colors.black),
                                    ('SIZE',(0,0),(-1,-1),8),
                                    ('FACE',(0,0),(0,2),'Helvetica-Bold'), # row head
                                    ('FACE',(0,0),(-1,-3),'Helvetica-Bold'), # column head
                                    ('FACE',(-1,-3),(-1,-1),'Helvetica-Bold'), # row tail
                                    ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                                    # ('LINEBEFORE',(3,0),(3,3),1,colors.black),
                                    # ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                                    ])

compcentileTS = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                            ('FACE',(0,0),(0,4),'Helvetica-Bold'), # row head
                            ('FACE',(-1,-4),(-1,-1),'Helvetica-Bold'), # row tail
                            ('LINEABOVE',(1,1),(-1,-4),1,colors.black),
                            ('LINEABOVE',(1,3),(-1,-2),.5,colors.black),
                            ('SIZE',(0,0),(-1,-1),8),
                            ('FACE',(0,0),(-1,-5),'Helvetica-Bold'), # column head
                            ('VALIGN',(0,0),(-1,-1),'TOP'),
                            ('TEXTCOLOR',(0,3),(-1,-1),colors.CMYKColor(*ncilightred)),
                            ('TEXTCOLOR',(0,1),(-1,-3),colors.CMYKColor(*ncilightblue)),
                            ])



grid = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),

                  ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                  ('BOX', (0,0), (-1,-1), 0.25, colors.black)])

concordance_tablestyle = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                            ('LINEABOVE',(1,1),(-1,-2),1,colors.black),
                            ('SIZE',(0,0),(-1,-1),8),
                            ('FACE',(-1,-2),(-1,-1),'Helvetica-Bold'),
                            ('FACE',(1,0),(-1,-3),'Helvetica-Bold'),
                            ('FACE',(0,0),(0,2),'Helvetica-Bold'),
                            ('VALIGN',(1,2),(-1,-1),'MIDDLE'),
                            # ('LINEBEFORE',(3,0),(-5,-1),1,colors.black),
                            ])
hwp_tablestyle = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                            ('LINEABOVE',(0,1),(-1,-3),1,colors.black),
                            ('SIZE',(0,0),(-1,-1),8),
                            ('FACE',(0,0),(-1,-4),'Helvetica-Bold'),
                            # ('FACE',(1,0),(-1,-4),'Helvetica-Bold'),
                            # ('FACE',(0,0),(0,3),'Helvetica-Bold'),
                            ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                            ])


def findcentiles(data, centiles):
  labels = [0 for i in centiles]
  for num,perc in data:
    for i,c in enumerate(centiles):
      if perc<=c:
        labels[i]=num
  return labels


def compute_percents(data, total):
  counts = []
  percents = []
  for d in data[:-1]:
    n=int(d[1])
    assert n<= total,'count (%d) >= total (%d)' % (n,total)
    counts.append(n)
    percents.append(percent(n,total))
  return counts,percents




def summarize_completion(samdata,locdata,frame,canvas):
  out = []
  out.append(header("Completion", style=HeaderStyle3))

  samples = list(islice(samdata,2,None))
  sampletotal = len(samples)
  loci = list(islice(locdata,2,None))
  locustotal = len(loci)



  samplecounts,samplepercents = compute_percents(samples,locustotal)
  locuscounts,locuspercents = compute_percents(loci,sampletotal)

  quantiles = [0,.05,.5,.95,1]

  samplecount_quantile  =  [quantile(sorted(samplecounts),  c, presorted=True) for c in quantiles]
  locuscount_quantile   =  [quantile(sorted(locuscounts),   c, presorted=True) for c in quantiles]

  samplepercent_quantile  = [percent(c,locustotal)  for c in samplecount_quantile]
  locuspercent_quantile   = [percent(c,sampletotal) for c in locuscount_quantile]

  tabledata = [['','min','5th','med','95th','max',''],
              ['Samples by Locus']  +  ['%d'     % c for c in locuscount_quantile]    + ['#'],
              ['']       +  ['%5.2f'  % p for p in locuspercent_quantile]  + ['%'],
              ['Loci by Sample'] +  ['%d'     % c for c in samplecount_quantile]   + ['#'],
              ['']       +  ['%5.2f'  % p for p in samplepercent_quantile] + ['%']]
  table = Table(tabledata, hAlign='LEFT', style=compcentileTS)

  height = 120
  width = 200

  graphcolors = [colors.CMYKColor(*ncilightred),colors.CMYKColor(*ncilightblue)]

  xmax = 1000

  quantiles = quantiles[1:4]

  labels = [[q*xmax for q in quantiles],['5th','med','95th']]

  plotcentiles = [i/(xmax*1.00) for i in xrange(xmax+1)]

  samplecount_centile =   [(i,quantile(sorted(samplecounts), c, presorted=True)) for i,c in enumerate(plotcentiles)]
  locuscount_centile  =   [(i,quantile(sorted(locuscounts),  c, presorted=True)) for i,c in enumerate(plotcentiles)]

  samplepercent_centile = [(i,percent(c,locustotal))  for i,c in samplecount_centile]
  locuspercent_centile  = [(i,percent(c,sampletotal)) for i,c in locuscount_centile]

  # locuslabels = [findcentiles(locuspercent_centile, locuspercent_quantile[1:4]),labels]

  x,y,w,h = frame.x1+frame.width-frame.height-(SPACER*2*inch),frame.y1+(inch*SPACER),frame.width,frame.height
  graph = lineplot([samplepercent_centile,locuspercent_centile],labels,graphcolors,x,y,width,height)

  t = Table([[table,graph]], style=completionTS, rowHeights=height)

  out.append(t)
  return out


def make_disconcordant_grid(name, data):
  concord,hethet,hethom,homhet,homhom = map(int, data[1])
  total = concord+hethet+hethom+homhet+homhom
  disconcord = hethet+hethom+homhet+homhom
  out = [name]
  out.append('%5.2f' % percent(hethet,total))
  out.append('%5.2f' % percent(hethom,total))
  out.append('%5.2f' % percent(homhet,total))
  out.append('%5.2f' % percent(homhom,total))
  out.append('%5.2f' % percent(disconcord,total))
  return out


def make_concordanttable(samsum, locsum):

  data =[['', 'HetHet','HetHom','HomHet','HomHom','Total']]
  loc =  make_disconcordant_grid('Loci',locsum)
  sam =  make_disconcordant_grid('Samples',samsum)
  data.append( loc )
  data.append( sam )

  t = Table(data)
  t.setStyle(concordance_tablestyle)
  return t

def make_concordant_centiletable(concordance, samsum, locsum):
  concord,hethet,hethom,homhet,homhom = map(int, samsum[1])
  total = concord+hethet+hethom+homhet+homhom

  samples = islice(read_section(concordance, section='sample_concordance'), 1, None)
  sam_conc = [int(s[2]) for s in samples]
  sam_tot = len(sam_conc)
  loci = islice(read_section(concordance, section='locus_concordance'), 1, None)
  loc_conc = [int(l[2]) for l in loci]
  loc_tot = len(loc_conc)

  sam_percents = [ (percent(c,loc_tot)) for c in sam_conc]
  loc_percents = [ (percent(c,sam_tot)) for c in loc_conc]

  centiles = [0,.05,.5,.95,1]
  samplecentiles = [quantile(sorted(sam_percents), c, presorted=True) for c in centiles]
  # samplecentiles = calculate_percentiles(sam_percents)
  samplecentiles = ['%5.2f' % f for f in samplecentiles]
  loccentiles = [quantile(sorted(loc_percents), c, presorted=True) for c in centiles]
  # loccentiles = calculate_percentiles(loc_percents)
  loccentiles = ['%5.2f' % f for f in loccentiles]

  out = [['min','5th','med','95th','max']]
  out.append(loccentiles)
  out.append(samplecentiles)
  t = Table(out)
  # t.setStyle(compcentileTS)
  return t


def sum_concordances(rows):
  conc=dc_hethet=dc_hethom=dc_homhet=dc_homhom=0
  for row in rows:
    conc += int(row[2])
    dc_hethet += int(row[3])
    dc_hethom += int(row[4])
    dc_homhet += int(row[5])
    dc_homhom += int(row[6])

  return  dict( (('concordant',conc),('hethet',dc_hethet),('hethom',dc_hethom),('homhet',dc_homhet),('homhom',dc_homhom)) )


def percent_concordances(data):
  total = 0
  for v in data.values():
    total+=v
  for k,v, in data.iteritems():
    yield percent(v,total)

def summarize_concordance(concordance):
  out = []
  out.append(header("Concordance", style=HeaderStyle3))

  if concordance:
    samcon = islice(read_section(concordance, section='sample_concordance'),1,None)
    samdata = sum_concordances(samcon)
    loccon = islice(read_section(concordance, section='locus_concordance'),1,None)
    locdata = sum_concordances(loccon)

    samplepercents = list(percent_concordances(samdata))
    locuspercents = list(percent_concordances(locdata))

    sampleconcordant = samdata.pop('concordant')
    locusconcordant = locdata.pop('concordant')

    counthead = ['']
    counthead.extend(samdata.keys())
    locusdata = ['Locus']
    locusdata.extend(locdata.values())
    sampledata = ['Sample']
    sampledata.extend(samdata.values())

    counts = [counthead,locusdata,sampledata]

    lightpiecolors = [colors.CMYKColor(*ncilightred),colors.CMYKColor(*darkncigrey),colors.CMYKColor(*lightnciblue),colors.CMYKColor(*darkncigrey),colors.CMYKColor(*ncidarkblue)]

    samplelabels = list([('%5.2f' % per) for per in samplepercents])
    locuslabels = list([('%5.2f' % per) for per in locuspercents])
    samplepie= pie(samplepercents,samplelabels,lightpiecolors,50)
    locuspie= pie(locuspercents,locuslabels,lightpiecolors,50)

    table =[[Table(counts), Table([['Sample','Locus'],[samplepie,locuspie]])]]

    out.append( Table(table, style=completionTS) )
  else:
    out.append(p(MISSING_INPUT))

  return out


def generate_pie(data,total,colors,frame,canvas):
  data = [num for num,status in data]
  labels = ['%2.1f' % percent(num,total) for num in data]

  pieheight = frame.height*.65
  p= pie(data,labels,colors,pieheight)

  x,y = frame.x1+frame.width-frame.height-(SPACER*inch),frame.y1+(inch*SPACER)
  drawingat(p,x,y,canvas)


def summarize_assay(assaydef,summary,frame,canvas):
  out = []

  if summary:
    summary = dict(summary)

    assays = list(islice(assaydef,1,None))
    total =  len(assays)
    valid = int(summary.get('total'))
    notvalid = total-valid

    ### Fix counting error when 'empty' is ''
    failed = int(summary.get('dropped'))
    if summary.get('empty') != '':
      failed += len(summary.get('empty').split(','))

    informative = total-failed-notvalid

    out.append(header("%d / %d Informative Loci" % (informative,total), style=HeaderStyle3))
    out.append(Spacer(.05*inch, .05*inch))

    data =[[failed,'loci failed in assay'],
            [informative,'informative loci']]
    if notvalid!=0:
      data.insert(0, [notvalid,'loci failed validation'])

    ts = [('TEXTCOLOR',(0,e),(-1,-1),color) for e,color in enumerate(ncicolors)]
    ts = ts+summaryTS
    t = Table(data, style=TableStyle(ts), rowHeights=9, hAlign='LEFT')
    out.append(t)

    generate_pie(data,total,ncilightcolors,frame,canvas)

  return out


summaryTS = [('ALIGN',(0,0),(-1,-1),'RIGHT'),('ALIGN',(1,0),(-1,-1),'LEFT'),('SIZE',(0,0),(-1,-1),7)]


def summarize_samples(sampledef,summary,frame,canvas):
  out = []

  summary = dict(summary)
  samples = list(islice(sampledef,1,None))
  total =  len(samples)

  samples = sorted(samples, key=itemgetter(1))
  groups = groupby(samples, itemgetter(1))
  statues = []
  notassayed = 0
  for status,content in groups:
    num = len(list(content))
    if status!='TRUE':
      notassayed+=num
      statues.append((status,num))

  statues = dict(statues)

  ### Fix counting error when 'empty' is ''
  failedassay = int(summary.get('dropped'))
  if summary.get('empty') != '':
    failedassay+=len(summary.get('empty').split(','))

  statues['failed genotyping']=failedassay

  informative = total-failedassay-notassayed
  out.append(header('%d / %d Informative Samples' % (informative,total), style=HeaderStyle3))
  out.append(Spacer(.05*inch, .05*inch))

  data = list([[num,status] for status,num in statues.iteritems()])
  # data.append([informative,'informative samples'])

  ts = [('TEXTCOLOR',(0,e),(-1,-1),color) for e,color in enumerate(ncicolors)]
  ts = list(chain(ts,summaryTS))
  t = Table(data, style=TableStyle(ts), rowHeights=9, hAlign='LEFT')

  out.append(t)

  data.append([informative,'informative'])

  generate_pie(data,total,ncilightcolors,frame,canvas)

  return out


def lineplot(data,labels,graphcolor,x,y,w,h):
  from reportlab.graphics.charts.lineplots import LinePlot

  lp = LinePlot()
  lp.x = 0
  lp.y = SPACER*3*inch
  lp.height = h*.95
  lp.width = w*.95
  lp.data = data
  lp.lines[0].strokeColor = graphcolor[0]
  lp.lines[1].strokeColor = graphcolor[1]

  lp.xValueAxis.visibleGrid = 1
  lp.xValueAxis.valueSteps = labels[0]
  lp.xValueAxis.labelTextFormat=labels[1]
  lp.xValueAxis.labels.fontSize = 8
  lp.xValueAxis.labels.textAnchor = 'end'

  lp.yValueAxis.valueMin = 0
  lp.yValueAxis.valueMax = 100
  lp.yValueAxis.valueStep = 20

  drawing = Drawing(w, h, hAlign='RIGHT')
  drawing.add(lp)
  return drawing


def drawingat(drawing,x,y,canvas):
  from reportlab.graphics import renderPDF
  renderPDF.draw(drawing, canvas, x,y)

def pie(data,labels,piecolors,diameter):
  from reportlab.graphics.charts.piecharts import Pie
  d = Drawing(diameter,diameter)
  pc = Pie()
  pc.startAngle = 240
  pc.x = pc.y = 0
  pc.width = pc.height = diameter
  pc.data = data
  pc.labels = labels
  pc.slices.strokeWidth=0.5
  pc.slices.labelRadius = 1.4
  pc.slices.fontSize = 8
  for i,c in enumerate(piecolors):
    pc.slices[i].fillColor = c
    pc.slices[i].fontColor = c

  d.add(pc)
  return d


def summarize_project(project):
  project = dict(project)

  out = []
  out.append(header(project.get('name'), style=HeaderStyle2))
  out.append(p(project.get('description'), MedParaStyle))

  comment = project.get('comment')
  if comment:
    comment = '<b>NOTE: </b>'+comment
    out.append(p(comment, SmallParaStyle))
  return out


def summarize_cgfproject(project):
  project = dict(project)

  out = []
  t = Table ([['project:',project.get('projectid'),''],
              ['study:',project.get('studyid'),''],
              ['cost center:',project.get('cas'),''],
              ['investigator:',project.get('investigator'),''],
              ['coordinator:',project.get('cgfcontact'),''],
              ['analyst:',project.get('cgfanalyst'),'']],
              style=[('ALIGN',(0,0),(-2,-1),'LEFT'),
                    ('ALIGN',(1,0),(-2,-1),'RIGHT'),
                    ('SIZE',(0,0),(-1,-1),10),
                    ('TEXTCOLOR',(0,0),(-1,-1),darkncigrey)],
              rowHeights=12)

  out.append(t)
  return out


def sort_pvalues(hwp):
  hwp = read_section(hwp, section='hwp')

  passed = []
  fail01 = []
  fail001 = []
  all = []
  for locus in hwp[1:]:
    loc,pval = locus[0],float(locus[5])
    if pval > 0.011:
      passed.append( (loc,pval) )
    if pval < 0.01:
      fail01.append( (loc,pval) )
    if pval < 0.001:
      fail001.append( (loc,pval) )
    all.append((loc,pval))

  return passed,fail01,fail001,all

def summarize_hwp(hwp):
  out = []
  out.append(header("Hardy–Weinberg Proportions", style=HeaderStyle3))

  if hwp:
    centiles = [0,.01,.05,.10,.50,.95,1]

    hwp = islice(read_section(hwp, section='hwp'),1,None)
    hwp = [(l[0],int(l[4]),float(l[5])) for l in hwp]

    #
    # TODO: NEED TO STORE THE TOTAL NUMBER OF SMAPLES IN hwp.dat
    total = max([c for l,c,per in hwp])

    counts = [quantile(sorted([c for l,c,per in hwp]), cen, presorted=True) for cen in centiles]
    # counts = list(calculate_percentiles([c for l,c,per in hwp], centiles))
    counts = ['%d' % c for c in counts]
    pvals = [quantile(sorted([per for l,c,per in hwp]), cen, presorted=True) for cen in centiles]
    # pvals = list(calculate_percentiles([per for l,c,per in hwp], centiles))
    pvals = ['%5.4f' % pval for pval in pvals]

    percents = [quantile(sorted([percent(c,total) for l,c,per in hwp]), cen, presorted=True) for cen in centiles]
    # percents = list(calculate_percentiles([percent(c,total) for l,c,per in hwp], centiles))
    percents = ['%5.2f' % per for per in percents]

    data = [['min','1st','5th','10th','med','95th','max'],
            counts,
            pvals,
            percents]
    # passed,fail01,fail001,all = sort_pvalues(hwp)
    # data =[['','Counts','%']]
    # cnts = '%d / %d ' % ( len(passed),len(all) )
    # data.append( ['passed (> 0.011)',cnts,'%5.2f' % percent(len(passed),len(all))] )
    # cnts = '%d / %d ' % ( len(fail01),len(all) )
    # data.append( ['failed (< 0.01)',cnts,'%5.2f' % percent(len(fail01),len(all))] )
    # cnts = '%d / %d ' % ( len(fail001),len(all) )
    # data.append( ['failed (< 0.001)',cnts,'%5.2f' % percent(len(fail001),len(all))] )

    t = Table(data, hAlign='LEFT', style=hwp_tablestyle)
    out.append(t)
  else:
    out.append( p(MISSING_INPUT) )

  return out


def summarize_dupcheck(dupcheck):
  out = []
  out.append(header("Duplicate Checking", style=HeaderStyle3))
  if dupcheck:
    headings = ['expected_duplicates','expected_duplicates_detected','unexpected_duplicates_detected','expected_duplicates_undetected']
    sections = filter_sections(read_sections(dupcheck),headings)
    sections = index_sections(sections)
    exp_dups = sections.get(headings[0])
    exp_dups_detect = sections.get(headings[1])
    unexp_dups_detect = sections.get(headings[2])
    exp_dups_undetect = sections.get(headings[3])

    if len(exp_dups)>1 or len(exp_dups_detect)>1 or len(unexp_dups_detect)>1 or len(exp_dups_undetect)>1:
      raise NoneUniqueSectionHeadingError, 'More than one section with the heading(s): %s' % ','.join(headings)
      return

    # if not exp_dups and not unexp_nondups and not exp_undetect_dups:
    #   out.append(p('No duplicates found'))

    if exp_dups[0]:
      out.append(p('%d / %d duplicates detected from genotypes' % (len(exp_dups_detect[0]),len(exp_dups[0]))))
    if unexp_dups_detect[0]:
      out.append(p('%s unexpected duplicate(s) detected from genotypes' % len(unexp_dups_detect[0])))
    if exp_dups_undetect[0]:
      out.append(p('%d expected duplicate(s) not detected from genotypes' % len(exp_dups_undetect[0])))
  else:
    out.append(p(MISSING_INPUT))

  return out


def summarize_future():
  out = []
  out.append(header("Contingency Association Analysis", style=HeaderStyle3))
  out.append(p(MISSING_INPUT))
  return out


def line(x1,y1,x2,y2):
  height = 15
  d = Drawing(1,height, hAlign='CENTER')
  l = Line(x1,y1+(height/2),x2,y2+(height/2), strokeColor=colors.black, strokeWidth=.5)
  d.add(l)
  # d.dumpProperties()
  return d


def combined_summary(projectdef,assaydef,sampledef, completion, dupcheck, concordance, hwp, filename, freetext):
  from reportlab.pdfgen import canvas
  from reportlab.lib.pagesizes import letter

  canvas = canvas.Canvas(filename+'.pdf', pagesize=letter)
  width, height = letter  #keep for later

  sp = SummaryPage(canvas, width, height)

  headings = ['project']
  sections = filter_sections(read_sections(projectdef),headings)
  sections = index_sections(sections)
  project = sections.get(headings[0])
  if len(project)>1:
    raise NoneUniqueSectionHeadingError, 'More than one section with the heading(s): %s' % ','.join(headings)
    return
  project = project[0]
  project = [x[:2] for x in project]
  loc_summ = loc_data = sam_summ = sam_data = None
  if completion:
    COMPLETION_SECTIONS=['summary', 'data']
    completion    =  index_sections(filter_sections(read_sections(completion),COMPLETION_SECTIONS))
    sam_summ = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Sample')
    loc_data = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Loci')
    loc_summ  = retrieve_section(completion, COMPLETION_SECTIONS[0], section_type='Loci')
    sam_data  = retrieve_section(completion, COMPLETION_SECTIONS[1], section_type='Sample')

  sp.project.addFromList(summarize_project(project),canvas)
  sp.cgfproject.addFromList(summarize_cgfproject(project),canvas)
  sp.sample.addFromList(summarize_samples(sampledef,sam_summ,sp.sample,canvas),canvas)
  sp.assay.addFromList(summarize_assay(assaydef,loc_summ,sp.assay,canvas),canvas)

  sp.analysis.addFromList(summarize_completion(sam_data,loc_data,sp.analysis,canvas),canvas)

  sp.analysis.addFromList([line(0,0,sp.analysis.width-(SPACER*inch),0)],canvas)

  if freetext:
    freetext = list(autofile(freetext))
    freetext = [p(t) for t in freetext]
    sp.analysis.addFromList(freetext,canvas)

  # sp.analysis.addFromList(summarize_concordance(concordance),canvas)
  #
  # sp.analysis.addFromList([line(0,0,sp.analysis.width-(SPACER*inch),0)],canvas)
  #
  # sp.analysis.addFromList(summarize_hwp(hwp),canvas)
  #
  # sp.analysis.addFromList([line(0,0,sp.analysis.width-(SPACER*inch),0)],canvas)
  #
  # sp.analysis.addFromList(summarize_dupcheck(dupcheck),canvas)
  #
  # sp.analysis.addFromList([line(0,0,sp.analysis.width-(SPACER*inch),0)],canvas)
  #
  # sp.analysis.addFromList(summarize_future(),canvas)

  canvas.showPage()
  canvas.save()


def gcombined_summary(completion, dupcheck, concordance, hwp, filename):
  from reportlab.pdfgen import canvas
  from reportlab.lib.pagesizes import letter

  canvas = canvas.Canvas(filename+'.pdf', pagesize=letter)
  width, height = letter  #keep for later

  summary = GoldenSummaryPage(canvas, width, height)

  summary.project.addFromList(summarize_project(completion),canvas)
  # summary.graphlegend.addFromList(graph_legend(),canvas)
  summary.assay.addFromList(summarize_assay(completion),canvas)

  if completion:
    comp = summarize_completion(completion)
    summary.summaries.addFromList(comp,canvas)
    summary.graph1.addFromList([p('completion graph')],canvas)
  if dupcheck:
    dups = summarize_dupcheck(dupcheck)
    summary.summaries.addFromList(dups,canvas)
  if concordance:
    concor = summarize_concordance(concordance)
    summary.summaries.addFromList(concor,canvas)
    summary.graph2.addFromList([p('concordance graph')],canvas)
  if hwp:
    hw = summarize_hwp(hwp)
    summary.summaries.addFromList(hw,canvas)
    summary.graph3.addFromList([p('hwp graph')],canvas)

  canvas.showPage()
  canvas.save()

def option_parser():
  import optparse

  usage = 'usage: %prog [options] [args]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__)

  parser.add_option('-s', '--projectsummary', dest='projectsummary', metavar='N', type='str',
                    help='option not defined')


  parser.add_option('-c', '--completion', dest='completion', metavar='FILE',
                    help='Tabular completion output using SIDs as sample identifiers')
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
                    help='Sample.prn file that uses SID as primary index')
  parser.add_option('-t', '--freetext', dest='freetext', metavar='FILE',
                    help='Freetext to place in the bottom portion')

  return parser


def main():
  parser = option_parser()
  options,args = parser.parse_args()

  assaydef = sampledef = completion = dupcheck = concordance = hwp = None
  if options.completion:
    completion = list(csv.reader(autofile(options.completion), dialect='excel-tab'))
  if options.dupcheck:
    dupcheck = list(csv.reader(autofile(options.dupcheck), dialect='excel-tab'))
  if options.concordance:
    concordance = list(csv.reader(autofile(options.concordance), dialect='excel-tab'))
  if options.hwp:
    hwp = list(csv.reader(autofile(options.hwp), dialect='excel-tab'))

  projectdef = assaydef = sampledef = None
  if options.projectdef:
    projectdef = list(csv.reader(autofile(options.projectdef), dialect='excel-tab'))
  if options.assaydef:
    assaydef = list(csv.reader(autofile(options.assaydef), dialect='excel-tab'))
  if options.sampledef:
    sampledef = list(csv.reader(autofile(options.sampledef), dialect='excel-tab'))

  if options.projectsummary:
    combined_summary(projectdef, assaydef, sampledef, completion, dupcheck, concordance, hwp, options.projectsummary, options.freetext)


if __name__ == '__main__':
  main()
