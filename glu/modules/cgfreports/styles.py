# -*- coding: utf-8 -*-
'''
File:          styles.py

Authors:

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

from    reportlab.lib.styles      import ParagraphStyle
from    reportlab.lib             import colors
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.lib.enums       import TA_RIGHT
from    reportlab.platypus        import TableStyle,Paragraph


styles = getSampleStyleSheet()

HeaderStyle1  = styles["Heading1"]
HeaderStyle2  = styles["Heading2"]
HeaderStyle3  = styles["Heading3"]
ParaStyle     = styles["Normal"]
PreStyle      = styles["Code"]

def p(txt,style=ParaStyle):
    return Paragraph(txt, style=style)

# COLORS
ncipurple       = .45,.36,.10,.01
nciblue         = .49,.32,.08,.01
ncigreen        = .52,.25,.29,.01
ncitan          = .49,.43,.49,.04
ncigrey         = .41,.34,.31,0
ncired          = .06,.85,.98,.18
ncicolors       = [colors.CMYKColor(*ncired),colors.CMYKColor(*ncigrey),colors.CMYKColor(*ncitan),colors.CMYKColor(*ncigreen),colors.CMYKColor(*nciblue),colors.CMYKColor(*ncipurple)]

ncilightpurple  = .30,.22,.07,.0
ncilightblue    = .43,.29,.09,0
ncilightgreen   = .33,.1,.15,.0
ncilighttan     = .23,.19,.26,.0
ncilightgrey    = .26,.20,.16,.0
ncilightred     = .24,.41,.29,.01
ncilightcolors  = [colors.CMYKColor(*ncilightred),colors.CMYKColor(*ncilightgrey),colors.CMYKColor(*ncilighttan),colors.CMYKColor(*ncilightgreen),colors.CMYKColor(*ncilightblue),colors.CMYKColor(*ncilightpurple)]

ncidarkblue     = .70,.48,.15,.04
darkncigrey     = .60,.51,.47,.04

HeaderStyle3.textColor = darkncigrey

# PARAGRAPH STYLES
SmallParaStyle  = ParagraphStyle('SmallParaStyle',parent=ParaStyle,
                                  fontSize=8, fontName='Times-Roman', leading=10,
                                  textColor=colors.CMYKColor(*darkncigrey))

MedParaStyle    = ParagraphStyle('MedParaStyle',parent=ParaStyle,
                                  fontSize=8, fontName='Times-Roman', leading=10,
                                  textColor=colors.CMYKColor(*darkncigrey))

titleStyle      = ParagraphStyle('titleStyle',parent=HeaderStyle1,
                                  fontSize=24, fontName='Helvetica-Bold',
                                  textColor=colors.white)

subtitleStyle   = ParagraphStyle('subtitleStyle',parent=ParaStyle,
                                  fontSize=10, fontName='Helvetica-Bold', alignment=TA_RIGHT,
                                  textColor=colors.white)

nciStyle        = ParagraphStyle('nciStyle',parent=ParaStyle,
                                  fontSize=10, fontName='Helvetica',
                                  textColor=colors.white)

detailStyle     = ParagraphStyle('detailStyle',parent=ParaStyle,fontSize=8)
defStyle        = ParagraphStyle('defStyle',parent=ParaStyle,leftIndent=8)


# TABLE STYLES
summaryTS         = [('ALIGN',(0,0),(-1,-1),'RIGHT'),('ALIGN',(1,0),(-1,-1),'LEFT'),('SIZE',(0,0),(-1,-1),8)]

headerTS          = TableStyle([('ALIGN',(0,0),(1,1),'LEFT'),
                                ('ALIGN',(-1,-1),(-2,-2),'RIGHT'),
                                ('SIZE',(0,0),(0,0),18),
                                ('FACE',(0,0),(0,0),'Times-Bold'),
                                ('FACE',(-1,-1),(-2,-2),'Helvetica-Bold'),
                                ])

completionTS      = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                                ('VALIGN',(0,0),(-1,-1),'TOP'),
                                ])

compcentileTS     = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
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

grid              = TableStyle([('ALIGN',(0,0),(-1,-1),'CENTER'),
                                ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                                ('BOX', (0,0), (-1,-1), 0.25, colors.black)])