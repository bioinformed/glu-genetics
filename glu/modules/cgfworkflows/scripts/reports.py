# -*- coding: utf-8 -*-
'''
File:          reportlib.py

Authors:       Brian Staats (staatsb@mail.nih.gov)

Created:       2006-07-17

Abstract:      Supporting report functions

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

import  os
import  sys
import  time

from    reportlab.platypus        import BaseDocTemplate,SimpleDocTemplate,Paragraph,Frame,Image,Table,Spacer,PageTemplate
from    reportlab.lib.styles      import getSampleStyleSheet
from    reportlab.rl_config       import defaultPageSize
from    reportlab.lib.units       import inch
from    reportlab.graphics.shapes import Line,Rect
from    reportlab.pdfgen          import canvas

from    glu.lib.utils             import percent



IMAGEDIR = '/home/staatsb/images/'
MISSING_INPUT = 'Analysis not performed. Requires data from investigator to complete.'

PAGE_HEIGHT=defaultPageSize[1]
PAGE_WIDTH=defaultPageSize[0]
SPACER=.125*inch
BORDER=SPACER*2



def retrieve_section(sections, heading, section_type=None):
  section = sections.get(heading)

  if not section_type and len(section) >1:
    raise NoneUniqueSectionHeadingError, 'More than one section with the heading: ' + heading

  if not section_type:
    return section[0]

  for items in section:
    for item in items:
      if section_type in set(item):
        return items



def main():
  pass

if __name__ == '__main__':
  main()
