# -*- coding: utf-8 -*-

from __future__ import absolute_import

__abstract__  = 'GLU application shell and launcher'
__copyright__ = 'Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'


import gc
import sys
import time
import optparse
import traceback

import pkgutil
import pkg_resources

import glu


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
  usage = 'usage: %prog [options] [module] [args...]'
  parser = optparse.OptionParser(usage=usage, version='%%prog %s' % __version__, add_help_option=False)
  parser.disable_interspersed_args()

  parser.add_option('-h', '--help', dest='help', action='store_true',
                    help='show this help message, then exit')
  parser.add_option('-s', '--stats', dest='stats', action='store_true',
                    help='display program runtime statistics')
  parser.add_option('-v', '--verbose', dest='verbose', action='store_true',
                    help='verbose error output')

  devopts = optparse.OptionGroup(parser, 'Options for software developers & power users')

  devopts.add_option('--path', dest='path', action='store_true',
                     help='Display GLU package installation path')
  devopts.add_option('-p', '--profile', dest='profile', action='store_true',
                     help='Profile GLU code to find performance bottlenecks')
  devopts.add_option('--profiler', dest='profiler', metavar='P', default='python',
                     help='Set the profiler to use when -p is specified')
  devopts.add_option('--gcstats', dest='gcstats', action='store_true',
                     help='Generate statistics from the runtime object garbage collector')
  devopts.add_option('--gcthreshold', dest='gcthreshold', metavar='N', type='int', default=1000000,
                     help='Set the threshold for triggering collection of generation-0 objects')

  parser.add_option_group(devopts)

  return parser


def main():
  parser = option_parser()
  glu_options,args = parser.parse_args()

  if glu_options.path:
    sys.stderr.write('GLU import path: %s\n\n' % glu.__file__)

  if glu_options.help or not args:
    parser.print_help(sys.stderr)
    sys.stderr.write('\nFor information on how to get started run the "intro" module,\n'
                     'usually as "glu intro".  For a list of available modules run\n'
                     '"glu list".\n\n')
    return

  if glu_options.gcstats:
    gc.set_debug(gc.DEBUG_STATS)

  if glu_options.gcthreshold >= 0:
    gc.set_threshold(glu_options.gcthreshold,100,10)

  module_name     = args[0]
  module_fullname = 'glu.modules.' + module_name
  module_options  = args[1:]

  try:
    loader = pkgutil.get_loader(module_fullname)

    if loader is None:
      raise ModuleMissingError()

    module = loader.load_module(module_fullname)

    module_info(module_name,module)

    try:
      progmain = module.main
    except AttributeError:
      sys.stderr.write('Unable to execute module %s.\n' % module_name)
      return

    if glu_options.stats:
      cstart = time.clock()
      tstart = time.time()
      sys.stderr.write('[%s] Execution start\n' % time.asctime())

    sys.argv = ['glu %s' % module_name] + module_options

    if glu_options.profile:
      run_profile(glu_options,progmain)
    else:
      progmain()

  except ModuleMissingError:
    sys.stderr.write('Unable to find module %s.  Please verify the module name and try again.\n' % module_name)
    return

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
