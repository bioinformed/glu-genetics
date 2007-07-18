# -*- coding: utf-8 -*-
'''
File:          launcher.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:

Abstract:

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

def run_profile(progmain,options,args):
  if not getattr(options,'profile',None):
    return progmain(options,args)

  if options.profile == 'python':
    try:
      import cProfile as profile
    except ImportError:
      import profile
    import pstats

    prof = profile.Profile()
    try:
      return prof.runcall(progmain, options, args)
    finally:
      stats = pstats.Stats(prof)
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)

  elif options.profile == 'hotshot':
    import hotshot, hotshot.stats
    prof = hotshot.Profile('tmp.prof')
    try:
      return prof.runcall(progmain, options, args)
    finally:
      prof.close()
      stats = hotshot.stats.load('tmp.prof')
      stats.strip_dirs()
      stats.sort_stats('time', 'calls')
      stats.print_stats(25)

  else:
    raise TagZillaError, 'ERROR: Unknown profiling option provided "%s"' % options.profile


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
  return ''.join(elapsed) or '0s'


def check_accelerators(accelerators, quiet=False):
    failed = []
    for a in accelerators:
      try:
        __import__(a)
      except ImportError:
        failed.append(a)

    if failed and not quiet:
      print 'WARNING: Failed to import the following native code accelerators:'
      print '            ',', '.join(failed)
      print '''\
         This will result in significantly increased computation times.
         Please do not post comparitive timing or benchmarking data when
         running in this mode.
'''
    return not failed


def launcher(progmain, opt_parser,
                       __program__      = '',
                       __version__      = '0.0',
                       __authors__      = [],
                       __copyright__    = '',
                       __license__      = '',
                       __accelerators__ = [],
                       **kwargs):

  if __program__ and __version__:
    print '%s version %s\n' % (__program__,__version__)
  elif __program__:
    print __program__
  elif __version__:
    print 'Version %s\n' % __version__

  if __authors__:
    print 'Written by %s' % __authors__[0]
    for author in __authors__[1:]:
      print '      and %s' % author
    print

  if __copyright__:
    print __copyright__
    print

  parser = opt_parser()
  options,args = parser.parse_args()

  if __license__ and options.license:
    print __license__
    return

  check_accelerators(__accelerators__)

  if options.help:
    parser.print_help(sys.stdout)
    return

  if not args:
    parser.print_usage(sys.stdout)
    print '''basic options:
      -h, --help            show detailed help on program usage'''
    if __copyright__ or __license__:
      print "      --license             show program's copyright and license terms and exit"
    print
    return

  start = time.clock()
  sys.stderr.write('[%s] Analysis start\n' % time.asctime())

  try:
    run_profile(progmain,options,args)

  except KeyboardInterrupt:
    sys.stderr.write('\n[%s] Analysis aborted by user\n' % time.asctime())

  except TagZillaError, e:
    sys.stderr.write('\n%s\n\n[%s] Analysis aborted due to reported error\n' % (e,time.asctime()))

  except IOError, e:
    sys.stderr.write('\n%s\n\n[%s] Analysis aborted due to input/output error\n' % (e,time.asctime()))

  except:
    import traceback
    sys.stderr.write('''
Analysis aborted due to a problem with the program input, parameters
supplied, an error in the program.  Please examine the following failure
trace for clues as to what may have gone wrong.  When in doubt, please send
this message and a complete description of the analysis you are attempting
to perform to the software developers.

Traceback:
  %s

[%s] Analysis aborted due to unhandled error
''' % (traceback.format_exc().replace('\n','\n  '),time.asctime()))

  else:
    sys.stderr.write('[%s] Analysis completed successfully\n' % time.asctime())
    sys.stderr.write('[%s] CPU time: %s\n' % (time.asctime(),format_elapsed_time(time.clock()-start)))
