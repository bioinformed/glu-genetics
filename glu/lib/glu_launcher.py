# -*- coding: utf-8 -*-

from __future__ import absolute_import

__abstract__  = 'GLU application shell and launcher'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import gc
import sys
import time
import argparse
import traceback

import pkgutil
import pkg_resources

import glu

from   glu.lib.glu_argparse import GLUArgumentParser


try:
  __version__ = pkg_resources.resource_string('glu','VERSION').strip()
except (IOError,ValueError):
  __version__ = '(unknown)'


class GLUError(RuntimeError): pass
class ModuleMissingError(ImportError): pass


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
        format = '%05.2f%s' if t else '%.2f%s'
      else:
        format = '%d%s'
      elapsed.append(format % (e,symbol))

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
      stats = pstats.Stats(prof,stream=sys.stderr)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
      stats.print_callers(25)

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
      stats.stream = sys.stderr
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)
      stats.print_callers(25)

  else:
    raise GLUError('ERROR: Unknown profiling option provided "%s"' % options.profiler)


def module_info(name,module,out=sys.stderr):
  program   = 'GLU %s module: %s' % (__version__,name)
  version   = getattr(module,'__version__',   None)
  authors   = getattr(module,'__authors__',   [])
  copyright = getattr(module,'__copyright__', None)

  if version:
    out.write('%s version %s\n' % (program,version))
  else:
    out.write(program + '\n')

  if authors:
    out.write('Written by %s\n' % authors[0])
    for author in authors[1:]:
      out.write('      and %s\n' % author)
    out.write('\n')

  if copyright:
    out.write(copyright)
    out.write('\n\n')


def option_parser():
  descr  = 'Driver program used to launch GLU modules\n'

  epilog = '''For information on how to get started run the "intro" module,'
              usually as "glu intro".  For a list of available modules run'
              "glu list".'''

  parser = GLUArgumentParser(description=descr, epilog=epilog, add_help=False)

  parser.add_argument('module', metavar='module', type=str, nargs='?',
                      help='GLU module to launch')
  parser.add_argument('module_args', nargs=argparse.REMAINDER,
                      help='module options (see specific module help for a list)')

  parser.add_argument('-h', '--help', action='store_true',
                      help='show this help message, then exit')
  parser.add_argument('-s', '--stats', action='store_true',
                      help='display program runtime statistics')
  parser.add_argument('-v', '--verbose', action='store_true',
                      help='verbose error output')
  parser.add_argument('--version', action='version', version=__version__)

  devopts = parser.add_argument_group('Options for software developers and power users')

  devopts.add_argument('--path', action='store_true',
                       help='Display GLU package installation path')
  devopts.add_argument('-p', '--profile', action='store_true',
                       help='Profile GLU code to find performance bottlenecks')
  devopts.add_argument('--profiler', metavar='P', default='python',
                       help='Set the profiler to use when -p is specified')
  devopts.add_argument('--gcstats', action='store_true',
                       help='Generate statistics from the runtime object garbage collector')
  devopts.add_argument('--gcthreshold', metavar='N', type=int, default=1000000,
                       help='Set the threshold for triggering collection of generation-0 objects')

  return parser


def write_traceback():
  sys.stderr.write('\nWell, this is embarrassing.\n')
  sys.stderr.write('\nTraceback:  %s\n' % (traceback.format_exc().replace('\n','\n  ')))


def main():
  parser  = option_parser()
  options = parser.parse_args()

  if options.path:
    sys.stderr.write('GLU interpreter version: %s\n' % sys.version)
    sys.stderr.write('GLU import path: %s\n\n' % glu.__file__)

  if options.help or not options.module:
    parser.print_help(sys.stderr)
    return 2

  if options.gcstats:
    gc.set_debug(gc.DEBUG_STATS)

  if options.gcthreshold >= 0:
    gc.set_threshold(options.gcthreshold,100,10)

  module_name     = options.module
  module_fullname = 'glu.modules.' + module_name
  module_options  = options.module_args

  if options.stats:
    cstart = time.clock()
    tstart = time.time()

  ret = 0

  try:
    loader = pkgutil.get_loader(module_fullname)

    if loader is None:
      raise ModuleMissingError()

    module = loader.load_module(module_fullname)

    module_info(module_name,module)

    try:
      progmain = module.main
    except AttributeError:
      sys.stdout.flush()
      sys.stderr.write('Unable to execute module %s.\n' % module_name)
      return

    if options.stats:
      sys.stdout.flush()
      sys.stderr.write('[%s] Execution start\n' % time.asctime())

    sys.argv = ['glu %s' % module_name] + module_options

    if options.profile:
      run_profile(options,progmain)
    else:
      progmain()

  except ModuleMissingError:
    sys.stdout.flush()
    sys.stderr.write('Unable to find module %s.  Please verify the module name and try again.\n' % module_name)
    return 1

  except KeyboardInterrupt:
    sys.stdout.flush()
    sys.stderr.write('\n[%s] Execution aborted by user\n' % time.asctime())

    if options.verbose:
      write_traceback()

    ret = 1

  except GLUError, e:
    sys.stdout.flush()
    sys.stderr.write('\n%s\n\n[%s] Execution aborted due to reported error\n' % (e,time.asctime()))

    if options.verbose:
      write_traceback()

    ret = 1

  except IOError, e:
    sys.stdout.flush()
    sys.stderr.write('\n%s\n\n[%s] Execution aborted due to input/output error\n' % (e,time.asctime()))

    if options.verbose:
      write_traceback()

    ret = 1

  except SystemExit, e:
    ret = e.code

  except BaseException, e:
    ret = 1

    sys.stdout.flush()
    sys.stderr.write('''\nWell, that could have gone better.\n''')

    if not options.verbose:
      sys.stderr.write('''
Execution aborted due to a problem with the program input, parameters
supplied, an error in the program.  Please examine the error message below.
If the cause of the error is not clear, re-run your command with 'glu -v'
for more information.

Error: %s
''' % e)

    else:
      sys.stderr.write('''
Execution aborted due to a problem with the program input, parameters
supplied, an error in the program.  Please examine the following failure
trace for clues as to what may have gone wrong.  When in doubt, please send
this message and a complete description of the analysis you are attempting
to perform to the software developers.

Error: %s

Command line:
  %s
''' % (e, ' '.join(sys.argv)))

      write_traceback()

    sys.stderr.write('\n[%s] Execution aborted due to a fatal error\n' % time.asctime())

  else:
    if options.stats:
      sys.stdout.flush()
      sys.stderr.write('[%s] Execution completed successfully\n' % time.asctime())

  if options.stats:
    sys.stdout.flush()
    sys.stderr.write('[%s] Clock time: %s, CPU time: %s\n' % (time.asctime(),
                                                              format_elapsed_time(time.time()  - tstart),
                                                              format_elapsed_time(time.clock() - cstart)))

  return ret


if __name__ == '__main__':
  sys.exit(main())
