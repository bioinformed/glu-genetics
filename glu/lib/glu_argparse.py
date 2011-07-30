# -*- coding: utf-8 -*-

__abstract__  = 'GLU argument parser'
__copyright__ = 'Copyright (c) 2011, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'

import sys
import argparse

from gettext import gettext as _


class GLUArgumentParser(argparse.ArgumentParser):
  def error(self, message):
    '''error(message: string)

    Prints help text and the message to stderr and exits.

    If you override this in a subclass, it should not return -- it
    should either exit or raise an exception.
    '''
    if message=='too few arguments':
      self.print_help()
    else:
      self.print_usage()

    self.exit(2, _('\n%s: command error: %s\n') % (self.prog, message))
