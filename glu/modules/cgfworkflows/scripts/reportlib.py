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


def retrieve_section(sections, heading, section_type=None):
  section = sections.get(heading)

  if not section_type and len(section) >1:
    raise NoneUniqueSectionHeadingError, 'More than one section with the heading: ' + heading

  if not section_type:
    return section

  for items in section:
    for item in items:
      if section_type in set(item):
        return items

def main():
  pass

if __name__ == '__main__':
  main()