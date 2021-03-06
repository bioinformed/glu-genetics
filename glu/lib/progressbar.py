#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
#
# progressbar  - Text progressbar library for python.
# Copyright (c) 2005 Nilton Volpato
# Copyright (c) 2009, BioInformed LLC
#
# Modified extensively for inclusion in GLU by Kevin Jacobs
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

from __future__ import division

'''Text progressbar library for python.

This library provides a text mode progressbar. This is typically used
to display the progress of a long running operation, providing a
visual clue that processing is underway.

The ProgressBar class manages the progress, and the format of the line
is given by a number of widgets. A widget is an object that may
display differently depending on the state of the progress. There are
three types of widget:
- a string, which always shows itself;
- a ProgressBarWidget, which may return a different value every time
it's update method is called; and
- a ProgressBarWidgetHFill, which is like ProgressBarWidget, except it
expands to fill the remaining width of the line.

The progressbar module is very easy to use, yet very powerful. And
automatically supports features like auto-resizing when available.
'''

__version__ = '2.4-glu'

import os
import sys
import time
import struct


try:
    from fcntl import ioctl
    import termios
except ImportError:
    pass

import signal


class ProgressBarWidget(object):
    '''This is an element of ProgressBar formatting.

    The ProgressBar object will call it's update value when an update
    is needed. It's size may change between call, but the results will
    not be good if the size changes drastically and repeatedly.
    '''
    def update(self, pbar):
        '''Returns the string representing the widget.

        The parameter pbar is a reference to the calling ProgressBar,
        where one can access attributes of the class for knowing how
        the update must be made.

        At least this function must be overridden.'''
        pass


class ProgressBarWidgetHFill(object):
    '''This is a variable width element of ProgressBar formatting.

    The ProgressBar object will call it's update value, informing the
    width this object must the made. This is like TeX \\hfill, it will
    expand to fill the line. You can use more than one in the same
    line, and they will all have the same width, and together will
    fill the line.
    '''
    def update(self, pbar, width):
        '''Returns the string representing the widget.

        The parameter pbar is a reference to the calling ProgressBar,
        where one can access attributes of the class for knowing how
        the update must be made. The parameter width is the total
        horizontal width the widget must have.

        At least this function must be overridden.'''
        pass


class ETA(ProgressBarWidget):
    'Widget for the Estimated Time of Arrival'
    def format_time(self, seconds):
        return time.strftime('%H:%M:%S', time.gmtime(seconds))
    def update(self, pbar):
        if pbar.currval == 0:
            return 'ETA:  --:--:--'
        elif pbar.finished:
            return 'Time: %s' % self.format_time(pbar.seconds_elapsed)
        else:
            elapsed = pbar.seconds_elapsed
            eta = elapsed * pbar.maxval / pbar.currval - elapsed
            return 'ETA:  %s' % self.format_time(eta)


class ProgressRate(ProgressBarWidget):
    'Widget for showing the transfer speed (useful for file transfers).'
    def __init__(self, unit='B'):
        self.unit = unit
        self.fmt = '%6.2f %s'
        self.prefixes = ['', 'K', 'M', 'G', 'T', 'P']
    def update(self, pbar):
        if pbar.seconds_elapsed < 2e-6:#== 0:
            bps = 0.0
        else:
            bps = pbar.currval / pbar.seconds_elapsed
        spd = bps
        for u in self.prefixes:
            if spd < 1000:
                break
            spd /= 1000
        return self.fmt % (spd, u + self.unit + '/s')


class RotatingMarker(ProgressBarWidget):
    'A rotating marker for filling the bar of progress.'
    def __init__(self, markers='|/-\\'):
        self.markers = markers
        self.curmark = -1
    def update(self, pbar):
        if pbar.finished:
            return self.markers[0]
        self.curmark = (self.curmark + 1)%len(self.markers)
        return self.markers[self.curmark]


class Percentage(ProgressBarWidget):
    'Just the percentage done.'
    def update(self, pbar):
        return '%3d%%' % pbar.percentage()


class SimpleProgress(ProgressBarWidget):
    "Returns what is already done and the total, e.g.: '5 of 47'"
    def update(self, pbar):
        return '%d of %d' % (pbar.currval, pbar.maxval)


class Completed(ProgressBarWidget):
    "Returns what is already done"
    def __init__(self, unit=None):
        self.unit = unit
    def update(self, pbar):
        if self.unit:
          return '%d %s' % (pbar.currval,self.unit)
        else:
          return '%d' % pbar.currval


class Bar(ProgressBarWidgetHFill):
    'The bar of progress. It will stretch to fill the line.'
    def __init__(self, marker='#', left='|', right='|'):
        self.marker = marker
        self.left = left
        self.right = right
    def _format_marker(self, pbar, width):
        if isinstance(self.marker, (str, unicode)):
            m = self.marker
        else:
            m = self.marker.update(pbar)
        return (m*(width//len(m)+1))[:width]

    def update(self, pbar, width):
        percent = pbar.percentage()
        cwidth = width - len(self.left) - len(self.right)
        marked_width = percent * cwidth // 100
        m = self._format_marker(pbar,marked_width)
        bar = (self.left + m.ljust(cwidth) + self.right)
        return bar


class ReverseBar(Bar):
    'The reverse bar of progress, or bar of regress. :)'
    def update(self, pbar, width):
        percent = pbar.percentage()
        cwidth = width - len(self.left) - len(self.right)
        marked_width = percent * cwidth // 100
        m = self._format_marker(pbar)
        bar = (self.left + (m*marked_width).rjust(cwidth) + self.right)
        return bar


default_widgets = [Percentage(), ' ', Bar()]
class ProgressBar(object):
    '''This is the ProgressBar class, it updates and prints the bar.

    The term_width parameter must be an integer or None. In the latter case
    it will try to guess it, if it fails it will default to 80 columns.

    The simple use is like this:
    >> pbar = ProgressBar().start()
    >> for i in xrange(100):
    ..    # do something
    ..    pbar.update(i+1)
    ..
    >> pbar.finish()

    But anything you want to do is possible (well, almost anything).
    You can supply different widgets of any type in any order. And you
    can even write your own widgets! There are many widgets already
    shipped and you should experiment with them.

    When implementing a widget update method you may access any
    attribute or function of the ProgressBar object calling the
    widget's update method. The most important attributes you would
    like to access are:
    - currval: current value of the progress, 0 <= currval <= maxval
    - maxval: maximum (and final) value of the progress
    - finished: True if the bar has finished (reached 100%), False o/w
    - start_time: the time when start() method of ProgressBar was called
    - seconds_elapsed: seconds elapsed since start_time
    - percentage(): percentage of the progress [0..100]. This is a method.

    The attributes above are unlikely to change between different versions,
    the other ones may change or cease to exist without notice, so try to rely
    only on the ones documented above if you are extending the progress bar.
    '''

    __slots__ = ('maxval', 'currval', 'term_width', 'start_time',
                 'last_update_time', 'seconds_elapsed', 'finished', 'fd',
                 'signal_set', 'widgets', 'update_interval', 'next_update',
                 'num_intervals')

    def __init__(self, maxval=100, widgets=default_widgets, term_width=None,
                 fd=sys.stderr):

        self.maxval     = maxval
        self.widgets    = widgets
        self.fd         = fd
        self.signal_set = False
        if term_width is not None:
            self.term_width = term_width
        else:
            try:
                self._handle_resize(None, None)
                signal.signal(signal.SIGWINCH, self._handle_resize)
                self.signal_set = True
            except (SystemExit, KeyboardInterrupt):
                raise
            except:
                raise
                self.term_width = int(os.environ.get('COLUMNS', 80)) - 1

        self.num_intervals = max(100, self.term_width)
        self.update_interval = self.maxval // self.num_intervals if self.maxval else 1
        self.next_update = 0

        self.currval = 0
        self.finished = False
        self.start_time = None
        self.last_update_time = None
        self.seconds_elapsed = 0

    def _handle_resize(self, signum, frame):
        s   = struct.pack('HHHH', 0, 0, 0, 0)
        h,w = struct.unpack('HHHH', ioctl(self.fd, termios.TIOCGWINSZ, s))[:2]
        self.term_width = w

    def percentage(self):
        'Returns the percentage of the progress.'
        return self.currval * 100 // self.maxval if self.maxval else None

    def _format_widgets(self):
        r = []
        hfill_inds = []
        num_hfill = 0
        currwidth = 0
        for i, w in enumerate(self.widgets):
            if isinstance(w, ProgressBarWidgetHFill):
                r.append(w)
                hfill_inds.append(i)
                num_hfill += 1
            elif isinstance(w, (str, unicode)):
                r.append(w)
                currwidth += len(w)
            else:
                weval = w.update(self)
                currwidth += len(weval)
                r.append(weval)
        for iw in hfill_inds:
            r[iw] = r[iw].update(self, (self.term_width - currwidth) // num_hfill)
        return r

    def _format_line(self):
        return ''.join(self._format_widgets()).ljust(self.term_width)

    def _next_update(self):
        if self.maxval:
          return int((int(self.num_intervals * self.currval / self.maxval) + 1) * self.update_interval)
        else:
          return self.currval

    def _need_update(self):
        '''Returns true when the progressbar should print an updated line.

        You can override this method if you want finer grained control over
        updates.

        The current implementation is optimized to be as fast as possible and
        as economical as possible in the number of updates. However, depending
        on your usage you may want to do more updates. For instance, if your
        progressbar stays in the same percentage for a long time, and you want
        to update other widgets, like ETA, then you could return True after
        some time has passed with no updates.

        Ideally you could call self._format_line() and see if it's different
        from the previous _format_line() call, but calling _format_line() takes
        around 20 times more time than calling this implementation of
        _need_update().
        '''
        return self.currval >= self.next_update

    def update(self, value, force=False, show=True):
        'Updates the progress bar to a new value.'
        #assert 0 <= value <= self.maxval
        self.currval = value
        if not force and not self._need_update():
            return
        now = time.time()
        self.seconds_elapsed = now - self.start_time
        self.next_update = self._next_update()
        if show:
          self.fd.write(self._format_line() + '\r')
        self.last_update_time = now

    def start(self,show=True):
        '''Start measuring time, and prints the bar at 0%.

        It returns self so you can use it like this:
        >> pbar = ProgressBar().start()
        >> for i in xrange(100):
        ..    # do something
        ..    pbar.update(i+1)
        ..
        >> pbar.finish()
        '''
        self.start_time = self.last_update_time = time.time()
        self.update(0,show=show)
        return self

    def finish(self, show=True):
        '''Used to tell the progress is finished.'''
        self.finished = True
        self.update(self.maxval or self.currval,show=show,force=True)
        if show:
          self.fd.write('\n')
        if self.signal_set:
          signal.signal(signal.SIGWINCH, signal.SIG_DFL)


class DummyProgressBar(object):
    '''This is the ProgressBar class, it updates and prints the bar.

    The term_width parameter must be an integer or None. In the latter case
    it will try to guess it, if it fails it will default to 80 columns.

    The simple use is like this:
    >> pbar = ProgressBar().start()
    >> for i in xrange(100):
    ..    # do something
    ..    pbar.update(i+1)
    ..
    >> pbar.finish()

    But anything you want to do is possible (well, almost anything).
    You can supply different widgets of any type in any order. And you
    can even write your own widgets! There are many widgets already
    shipped and you should experiment with them.

    When implementing a widget update method you may access any
    attribute or function of the ProgressBar object calling the
    widget's update method. The most important attributes you would
    like to access are:
    - currval: current value of the progress, 0 <= currval <= maxval
    - maxval: maximum (and final) value of the progress
    - finished: True if the bar has finished (reached 100%), False o/w
    - start_time: the time when start() method of ProgressBar was called
    - seconds_elapsed: seconds elapsed since start_time
    - percentage(): percentage of the progress [0..100]. This is a method.

    The attributes above are unlikely to change between different versions,
    the other ones may change or cease to exist without notice, so try to rely
    only on the ones documented above if you are extending the progress bar.
    '''

    __slots__ = ()

    def percentage(self):
        'Returns the percentage of the progress.'
        return 100.

    def update(self, value, force=False):
        'Updates the progress bar to a new value.'

    def start(self):
        '''Start measuring time, and prints the bar at 0%.

        It returns self so you can use it like this:
        >> pbar = ProgressBar().start()
        >> for i in xrange(100):
        ..    # do something
        ..    pbar.update(i+1)
        .
        >> pbar.finish()
        '''
        return self

    def finish(self):
        '''Used to tell the progress is finished.'''


def _isatty(f):
  if hasattr(f,'isatty'):
    return f.isatty()
  elif hasattr(f,'fileno'):
    return os.isatty(f.fileno())
  return False


def progress_bar(output=None, widgets=None, length=None, label='Progress: ', bar_marker='#', units=None):
  if output is None:
    output = sys.stderr

  if not _isatty(output):
    return None

  if widgets is None:
    widgets = []

    if label:
      widgets += [label, ' ']

    if length is not None:
      widgets += [Percentage(), ' ', Bar(marker=bar_marker), ' ', ETA()]
    else:
      widgets += ['  ',Completed(unit=units)]

    if units:
      widgets += [' ',ProgressRate(unit=units)]

  return ProgressBar(widgets=widgets, maxval=length)


try:
  import posix

  def _find_tty():
    for f in (sys.stdin,sys.stderr,sys.stdout):
      if _isatty(f):
        return f
    return open('/dev/tty')

  # NB: The tty is kept open and never closed for performance reasons.
  #     Hopefully this won't cause assplosions.
  posix_tty    = _find_tty()
  posix_tty_fd = posix_tty.fileno()

  def is_foreground():
    tpgrp = posix.tcgetpgrp(posix_tty_fd)
    pgrp  = posix.getpgrp()
    return tpgrp == pgrp

except (ImportError,IOError):
  def is_foreground():
    return True


try:
  class IntervalTimer(object):
    def __init__(self, interval=1):
      self.interval = interval
      self.jobs = {}
      self.seq = 0
      self.old_handler = signal.signal(signal.SIGALRM, self._handler)
      self.old_alarm = 0

    def add(self, handler):
      jobid = self.seq
      self.seq += 1
      hasjobs = bool(self.jobs)

      self.jobs[jobid] = handler

      if not hasjobs:
        self.old_alarm = signal.alarm(self.interval)

      return jobid

    def remove(self, jobid):
      del self.jobs[jobid]
      if not self.jobs:
        signal.alarm(self.old_alarm)

    def _handler(self, signum, frame):
      for job in self.jobs.itervalues():
        job(signum, frame)

      if self.old_handler and self.old_alarm>0:
        self.old_alarm -= self.interval
        if self.old_alarm<=0:
          self.old_handler(signum, frame)
          self.old_alarm = 0

      if self.jobs:
        signal.alarm(self.interval)


  interval_timer = IntervalTimer(1)


  def progress_loop(items, update_interval=None, length=None, **kwargs):
    if length is None:
      try:
        length = len(items)
        update_interval = max(1,min(length//100,250))
      except TypeError:
        pass

    bar = progress_bar(length=length,**kwargs)

    if bar is None:
      return items

    def _progress():
      i=0

      bar.start(show=is_foreground())

      def progress_handler(signum, frame):
        # Do not print status if we're in the background
        # Uses the current value of the loop counter
        if is_foreground():
          bar.update(i, force=True)

      jobid = interval_timer.add(progress_handler)

      try:
        for i,item in enumerate(items):
          yield item

      finally:
        interval_timer.remove(jobid)

      bar.update(i)
      bar.finish(show=is_foreground())

    return _progress()

except ImportError:

  # Fall back to a more portable dumb implementation
  def progress_loop(items, length=None, update_interval=None, **kwargs):
    if length is None:
      try:
        length = len(items)
        update_interval = max(1,min(length//100,250))
      except TypeError:
        pass

    bar = progress_bar(length=length,**kwargs)

    if bar is None:
      return items

    def _progress():
      bar.start()

      for i,item in enumerate(items):
        # Do not print status if we're in the background
        if i%update_interval == 0 and is_foreground():
          bar.update(i, force=True)

        yield item

      bar.update(i)
      bar.finish()

    return _progress()
