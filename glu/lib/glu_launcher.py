# -*- coding: utf-8 -*-
'''
File:          glu_lancher.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-06-10

Abstract:      Application runner

Requires:      Python 2.5, glu

Revision:      $Id$
'''

from __future__ import absolute_import

__version__   = '0.5'
__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'


import sys
import time
import optparse
import traceback

import glu

class GLUError(RuntimeError): pass


def format_elapsed_time(t):
  units = [('s',60),('m',60),('h',24),('d',365),('y',None)]

  elapsed = []
  for symbol,divisor in units:
    if divisor:
      t,e = divmod(t,divisor)
    else:
      t,e = 0,t

    if e:
      if symbol == 's':
        elapsed.append('%.2f%s' % (e,symbol))
      else:
        elapsed.append('%d%s' % (e,symbol))

    if not t:
      break

  elapsed.reverse()
  return ''.join(elapsed) or '0.00s'


def run_profile(options,progmain):
  if options.profiler == 'python':
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      return prof.runcall(progmain)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)

  elif options.profiler == 'hotshot':
    import hotshot, hotshot.stats

    if sys.argv:
      statfile = '%s.prof' % sys.argv[0]
    else:
      statfile = 'tmp.prof'

    prof = hotshot.Profile(statfile)

    try:
      return prof.runcall(progmain)
    finally:
      prof.close()
      stats = hotshot.stats.load(statfile)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)

  else:
    raise GLUError('ERROR: Unknown profiling option provided "%s"' % options.profiler)


def module_info(name,module,out=sys.stderr):
  program   = 'GLU module: %s' % name
  version   = getattr(module,'__version__',   None)
  authors   = getattr(module,'__authors__',   [])
  copyright = getattr(module,'__copyright__', None)

  if version:
    sys.stderr.write('%s version %s\n' % (program,__version__))
  else:
    sys.stderr.write(program + '\n')

  if authors:
    sys.stderr.write('Written by %s\n' % authors[0])
    for author in authors[1:]:
      sys.stderr.write('      and %s\n' % author)
    sys.stderr.write('\n')

  if copyright:
    sys.stderr.write(copyright)
    sys.stderr.write('\n\n')



def main():
  usage = 'usage: %prog [options] [module] [args...]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__, add_help_option=False)
  parser.disable_interspersed_args()

  parser.add_option('-h', '--help', dest='help', action='store_true',
                    help='show this help message, then exit')
  parser.add_option('-s', '--stats', dest='stats', action='store_true',
                    help='display program runtime statistics')
  parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
                    help='verbose error output')
  parser.add_option('--path', dest='path', action='store_true',
                    help='Display GLU path')
  parser.add_option('-p', '--profile', dest='profile', action='store_true', help=optparse.SUPPRESS_HELP)
  parser.add_option('--profiler', dest='profiler', metavar='P', default='python', help=optparse.SUPPRESS_HELP)

  # options:
  #   set module path
  #   list modules
  #   profile
  #   verbose

  glu_options,args = parser.parse_args()

  if glu_options.path:
    sys.stderr.write('GLU import path: %s\n\n' % glu.__file__)

  if glu_options.help or not args:
    parser.print_help(sys.stderr)
    sys.stderr.write('\nFor information on how to get started run the "intro" module,\n'
                     'usually as "glu intro".  For a list of available modules run\n'
                     '"glu list".\n\n')
    return

  modulename = args[0]
  module_options = args[1:]

  try:
    module = __import__('glu.modules.' + modulename, fromlist=['main'])
  except ImportError:
    sys.stderr.write('Unable to import module %s.  Please verify module name and try again.\n' % modulename)
    return

  try:
    progmain = module.main
  except AttributeError:
    sys.stderr.write('Unable to execute module %s.\n' % modulename)
    return

  sys.argv = ['glu %s' % modulename] + module_options

  module_info(modulename,module)

  if glu_options.stats:
    cstart = time.clock()
    tstart = time.time()
    sys.stderr.write('[%s] Execution start\n' % time.asctime())

  try:
    if glu_options.profile:
      run_profile(glu_options,progmain)
    else:
      progmain()

  except KeyboardInterrupt:
    sys.stderr.write('\n[%s] Execution aborted by user\n' % time.asctime())

    if glu_options.verbose:
      sys.stderr.write('\nTraceback:  %s\n' % (traceback.format_exc().replace('\n','\n  ')))

  except GLUError, e:
    sys.stderr.write('\n%s\n\n[%s] Execution aborted due to reported error\n' % (e,time.asctime()))

    if glu_options.verbose:
      sys.stderr.write('\nTraceback:  %s\n' % (traceback.format_exc().replace('\n','\n  ')))

  except IOError, e:
    sys.stderr.write('\n%s\n\n[%s] Execution aborted due to input/output error\n' % (e,time.asctime()))

    if glu_options.verbose:
      sys.stderr.write('\nTraceback:  %s\n' % (traceback.format_exc().replace('\n','\n  ')))

  except SystemExit:
    pass
  except:
    sys.stderr.write('''
Execution aborted due to a problem with the program input, parameters
supplied, an error in the program.  Please examine the following failure
trace for clues as to what may have gone wrong.  When in doubt, please send
this message and a complete description of the analysis you are attempting
to perform to the software developers.

Traceback:
  %s

[%s] Execution aborted due to unhandled error
''' % (traceback.format_exc().replace('\n','\n  '),time.asctime()))

  else:
    if glu_options.stats:
      sys.stderr.write('[%s] Execution completed successfully\n' % time.asctime())

  if glu_options.stats:
    sys.stderr.write('[%s] Clock time: %s, CPU time: %s\n' % (time.asctime(),
                                                              format_elapsed_time(time.time()  - tstart),
                                                              format_elapsed_time(time.clock() - cstart)))


if __name__ == '__main__':
  main()
