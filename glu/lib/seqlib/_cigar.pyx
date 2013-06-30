# -*- coding: utf-8 -*-

from __future__ import division

__abstract__  = 'Helper functions for SAM/BAM format encoding and decoding'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'


from libc.stdio    cimport sprintf
from libc.stdlib   cimport malloc, realloc, free
from cpython       cimport PyString_Size, PyErr_NoMemory


cdef inline size_t imax(int a, int b):
  if a < b:
    return b
  else:
    return a


cdef inline void reserve_str(char **buffer, Py_ssize_t *size, Py_ssize_t extra, Py_ssize_t reserve):
  if extra>=size[0]:
    buffer[0] = <char *>realloc(buffer[0], reserve*sizeof(char))
    size[0]   = reserve


cdef inline void append_int(char **buffer, Py_ssize_t *size, Py_ssize_t *pos, Py_ssize_t value, Py_ssize_t reserve):
  reserve_str(buffer, size, pos[0]+10, reserve)
  pos[0] += sprintf(buffer[0]+pos[0], "%d", value)


def make_cigar(query, ref):
  cdef char  symbol = 0, new_symbol = 0
  cdef Py_ssize_t i, pos = 0, count = 0

  cdef char      *q = query
  cdef char      *r = ref
  cdef Py_ssize_t n = PyString_Size(query)
  cdef Py_ssize_t m = PyString_Size(ref)

  if n!=m:
    raise ValueError('query and reference must be of equal length')

  # Allocate max(n,20) characters for the CIGAR result string
  # More space can be added to accommodate very long results (not likely for good alignments)
  cdef Py_ssize_t size = imax(n,20)
  cdef char    *result = <char *>malloc(size*sizeof(char))

  for i in range(n):
    if q[i]=='-' and r[i]=='-':
      new_symbol = 'P'
    elif q[i]=='-':
      new_symbol = 'D'
    elif r[i]=='-':
      new_symbol = 'I'
    else:
      new_symbol = 'M'

    # When a new symbol is seen, write the previous RLE count and symbol, if any
    if count and symbol!=new_symbol:
      # Write the ASCII count and symbol to the result buffer and reset count to zero
      # (NULL termination of result will occur after the final loop iteration)
      append_int(&result, &size, &pos, count, size*2)

      result[pos] = symbol
      pos += 1
      count = 0

    count += 1
    symbol = new_symbol

  # Write remaining RLE count and symbol
  if count:
    # Write the ASCII count and symbol to the result buffer
    # (NULL termination of result will occur unconditionally next)
    append_int(&result, &size, &pos, count, size+10)
    result[pos] = symbol
    pos += 1

  # NULL terminate result
  result[pos] = 0

  # Convert result to a Python string, free memory, and return
  cdef bytes retval = result
  free(result)

  return retval


def make_ndiff(query,ref):
  cdef Py_ssize_t i, pos=0, end=0, eq=0, nm=0
  cdef char      *q = query
  cdef char      *r = ref
  cdef Py_ssize_t n = PyString_Size(query)
  cdef Py_ssize_t m = PyString_Size(ref)

  if n!=m:
    raise ValueError('query and reference must be of equal length')

  # Allocate 50 characters for the result MD string
  # More space can be added to accommodate very long results (not likely for good alignments)
  cdef Py_ssize_t size = 50
  cdef char *ndiff = <char *>malloc(size*sizeof(char))

  # Implement a finite state transducer
  while pos<n:
    # Find runs of gaps: ignore
    if q[pos]==r[pos]=='-':
      while pos<n and q[pos]==r[pos]=='-':
        pos += 1

    # Find runs of matches: count matches (eq)
    elif q[pos]==r[pos]!='N':
      while pos<n and 'N'!=q[pos]==r[pos]!='-':
        pos += 1
        eq  += 1

    # Find runs of insertions: count mismatches (nm)
    elif r[pos]=='-':
      while pos<n and q[pos]!='-' and r[pos]=='-':
        nm  += 1
        pos += 1

    # Find runs of deletions: write match count (eq), output each deleted base
    elif q[pos]=='-':
      append_int(&ndiff, &size, &end, eq, size*2)
      ndiff[end] = '^'
      end += 1
      eq   = 0

      while pos<n and q[pos]=='-' and r[pos]!='-':
        reserve_str(&ndiff, &size, end+1, size*2)
        ndiff[end] = r[pos]
        end += 1
        nm  += 1
        pos += 1

    # Find runs of mismatches: write match count (eq) for each mismatch base
    else:
      while pos<n and ('-'!=q[pos]!=r[pos]!='-' or 'N'==q[pos]==r[pos]):
        append_int(&ndiff, &size, &end, eq, size*2)

        ndiff[end] = r[pos]
        eq   = 0
        end += 1
        nm  += 1
        pos += 1

  # Write final number of matched bases (eq)
  append_int(&ndiff, &size, &end, eq, size+10)

  # Unconditionally null-terminate
  ndiff[end] = 0

  # Copy result ndiff to a Python string and free buffer
  cdef bytes md = ndiff
  free(ndiff)

  # Return count of mismatches and ndiff string (md)
  return nm,md
