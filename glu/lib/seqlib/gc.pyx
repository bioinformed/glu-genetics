# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Helper functions for fast GC/CpG estimation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id: $'

cimport cython

from libc.stdio  cimport sprintf
from libc.stdlib cimport malloc, realloc, free
from cpython     cimport PyString_Size, PyErr_NoMemory


cdef class gc_window:
  cdef str        py_seq
  cdef Py_ssize_t slen, start, end, half, gc, valid, win
  cdef int        state
  cdef char*      seq

  def __init__(self, seq, win):
    cdef char base

    self.py_seq = seq
    self.seq    = seq
    self.slen   = len(seq)
    self.start  = 0
    self.end    = 0
    self.half   = (win+1)//2
    self.win    = win
    self.gc     = 0
    self.valid  = 0
    self.state  = 0

    # Initialize window
    while self.end < self.half and self.end < self.slen:
      base = self.seq[self.end]
      #print 'Adding base:',chr(base)
      if base=='G' or base=='C':
        self.gc    += 1
        self.valid += 1
      elif base=='A' or base=='T':
        self.valid += 1

      self.end += 1

    if self.win>=self.slen:
      self.state = 2

  def __iter__(self):
    return self

  def __next__(self):
    cdef char base

    if self.state>2:
      raise StopIteration

    #print self.state,self.start,self.end,self.seq[self.start:self.end],self.gc
    pgc = self.gc/self.valid if self.valid else 0.

    if self.state==0:
      # Remove trailing base if window shifts
      if self.end-self.start==self.win:
        assert self.start<self.slen
        base = self.seq[self.start]
        self.start += 1
        #print 'Removing base:',chr(base)
        if base=='G' or base=='C':
          self.gc    -= 1
          self.valid -= 1
        elif base=='A' or base=='T':
          self.valid -= 1

      # Add next base
      if self.end<self.slen:
        base = self.seq[self.end]
        #print 'Adding base:',chr(base)
        if base=='G' or base=='C':
          self.gc    += 1
          self.valid += 1
        elif base=='A' or base=='T':
          self.valid += 1

        self.end += 1

      if self.end==self.slen:
        self.state = 1

    elif self.state==1:
      assert self.start<self.slen
      base = self.seq[self.start]
      #print 'Removing base:',chr(base)
      if base=='G' or base=='C':
        self.gc    -= 1
        self.valid -= 1
      elif base=='A' or base=='T':
        self.valid -= 1

      self.start += 1

      if self.start>=self.slen-self.win+self.half:
        self.state = 3

    elif self.state==2:
      self.start += 1
      if self.start==self.slen:
        self.state = 3

    return pgc


cdef class cpg_window:
  cdef str        py_seq
  cdef Py_ssize_t slen, start, end, half, gc, cpg, valid, win
  cdef int        state
  cdef char*      seq

  def __init__(self, seq, win):
    cdef char base, last

    self.py_seq = seq
    self.seq    = seq
    self.slen   = len(seq)
    self.start  = 0
    self.end    = 0
    self.half   = (win+1)//2
    self.win    = win
    self.cpg    = 0
    self.gc     = 0
    self.valid  = 0
    self.state  = 0

    # Initialize window
    last = 0
    while self.end < self.half and self.end < self.slen:
      base = self.seq[self.end]
      if base=='G' and last=='C':
        self.cpg   += 1
      if base=='G' or base=='C':
        self.gc    += 1
        self.valid += 1
      elif base=='A' or base=='T':
        self.valid += 1

      last      = base
      self.end += 1

    if self.win>=self.slen:
      self.state = 2

  def __iter__(self):
    return self

  def __next__(self):
    cdef char  base, last

    if self.state>2:
      raise StopIteration

    #print self.state,self.start,self.end,self.seq[self.start:self.end],self.gc
    pgc  = self.gc / self.valid    if self.valid else 0.
    pcpg = self.cpg/(self.valid/2) if self.valid else 0.

    if self.state==0:
      # Remove trailing base if window shifts
      if self.end-self.start==self.win:
        assert self.start<self.slen
        base = self.seq[self.start]
        last = self.seq[self.start-1] if self.start else 0

        self.start += 1
        #print 'Removing base:',chr(base)
        if base=='G' and last=='C':
          self.cpg   -= 1
        if base=='G' or base=='C':
          self.gc    -= 1
          self.valid -= 1
        elif base=='A' or base=='T':
          self.valid -= 1

      # Add next base
      if self.end<self.slen:
        base = self.seq[self.end]
        last = self.seq[self.end-1] if self.end else 0
        #print 'Adding base:',chr(base)
        if base=='G' and last=='C':
          self.cpg   += 1
        if base=='G' or base=='C':
          self.gc    += 1
          self.valid += 1
        elif base=='A' or base=='T':
          self.valid += 1

        self.end += 1

      if self.end==self.slen:
        self.state = 1

    elif self.state==1:
      assert self.start<self.slen
      base = self.seq[self.start]
      last = self.seq[self.start-1] if self.start else 0
      #print 'Removing base:',chr(base)
      if base=='G' and last=='C':
        self.cpg   -= 1
      if base=='G' or base=='C':
        self.gc    -= 1
        self.valid -= 1
      elif base=='A' or base=='T':
        self.valid -= 1

      self.start += 1

      if self.start>=self.slen-self.win+self.half:
        self.state = 3

    elif self.state==2:
      self.start += 1
      if self.start==self.slen:
        self.state = 3

    return pgc,pcpg
