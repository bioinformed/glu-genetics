from   operator  import getitem, itemgetter
from   itertools import izip, imap, groupby, repeat


CIGAR_map = { ('-','-'):'P' }
for a in 'NACGTacgt':
  CIGAR_map[a,'-'] = 'I'
  CIGAR_map['-',a] = 'D'
  for b in 'NACGTacgt':
    CIGAR_map[a,b] = 'M'


NDIFF_map = { ('-','-'): ('P',None) }
for a in 'NACGTacgt':
  NDIFF_map[a,'-'] = ('I',a)
  NDIFF_map['-',a] = ('D',a)
  # N,N is considered a mismatch(?)
  for b in 'NACGTacgt':
    NDIFF_map[a,b] = ('=' if a==b and a!='N' else 'X',b)


def make_cigar_py(query,ref):
  assert len(query)==len(ref)
  igar  = imap(getitem, repeat(CIGAR_map), izip(query,ref))
  cigar = ''.join('%d%s' % (len(list(run)),code) for code,run in groupby(igar))
  return cigar


def make_ndiff_py(query,ref):
  assert len(query)==len(ref)

  nd   = groupby(imap(getitem, repeat(NDIFF_map), izip(query,ref)),itemgetter(0))
  nm   = 0
  md   = ''
  eq   = 0

  for code,items in nd:
    if code=='P':
      continue
    elif code=='=':
      eq += len(list(items))
    elif code=='I':
      nm += len(list(items))
    else:
      md += '%d' % eq
      eq = 0

      bases = ''.join(b for c,b in items)
      nm   += len(bases)

      if code=='D':
        md   += '^'
      else:
        # For some silly reason mismatch bases must be 0 separated
        bases = '0'.join(bases)

      md += bases

  md += '%d' % eq

  return nm,md


# Try to import the optimized Cython version
# The Python version is pretty fast, but I wanted to play with Cython.
try:
  from glu.lib.seqlib._cigar import make_cigar, make_ndiff
except ImportError:
  make_cigar = make_cigar_py
  make_ndiff = make_ndiff_py


def cigar_add_trim(cigar, trim_char, left, right):
  if left and right:
    cigar = '%d%s%s%d%s' % (left,trim_char,cigar,right,trim_char)
  elif left:
    cigar = '%s%s%s' % (left,trim_char,cigar)
  elif right:
    cigar = '%s%d%s' % (cigar,right,trim_char)
  return cigar
