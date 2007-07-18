# -*- coding: utf-8 -*-
'''
File:          bitarray.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-04-10

Abstract:      Functions to get and set bits in a byte array

Requires:      Python 2.5

Revision:      $Id: $
'''

__version__ = '0.1'
__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC"). All rights reserved.'


__all__ = ['getbits','setbits']


def getbits_python(data, i, readlen):
  '''
  From data, read bits [i,i+readlen)

  @param    data: input data
  @type     data: byte array
  @param       i: bit index into data from which to begin reading
  @type        i: int
  @param readlen: number of bits to read in range 1..32
  @type  readlen: int
  @return       : integer comprised of bits [i,i+readlen) of data as the least significant bits
  @rtype        : int
  '''
  if i < 0:
    i += len(data)*8

  if i < 0:
    raise IndexError('read index out of range')

  if readlen <= 0 or readlen > 32:
    raise IndexError('read length must be between 1..32')

  if len(data) <= (i+readlen-1)>>3:
    raise IndexError('read length out of range')

  result = 0

  byte,offset = divmod(i,8)

  while readlen:
    b = data[byte]
    n = 8 - offset

    # mask high order bits
    if n != 8:
      b &= (1<<n)-1

    # shift off low order bits (the tail)
    tail = n - readlen
    if tail > 0:
      b >>= tail
      n  -= tail

    result   = (result<<n)|b
    readlen -= n
    byte    += 1
    offset   = 0

  return result


def setbits_python(data, i, value, writelen):
  '''
  From data, read bits [i,i+readlen)

  @param    data: input data
  @type     data: byte array
  @param       i: bit index into data from which to begin reading
  @type        i: int
  @param readlen: number of bits to read in range 1..32
  @type  readlen: int
  @return       : integer comprised of bits [i,i+readlen) of data as the least significant bits
  @rtype        : int
  '''
  if i < 0:
    i += len(data)*8

  if i < 0:
    raise IndexError('write index out of range')

  if writelen <= 0 or writelen > 32:
    raise IndexError('write length must be between 1..32')

  if len(data) <= (i+writelen-1)>>3:
    raise IndexError('write length out of range')

  byte,offset = divmod(i,8)

  while writelen:
    n          = 8 - offset
    tail       = max(0,n-writelen)
    n         -= tail
    mask       = ((1<<n)-1)<<tail
    writelen  -= n
    # Mask (xor) the old values, then set (or) the new bits
    data[byte] = (data[byte]&(~mask)) | ((value>>writelen<<tail)&mask)
    byte      += 1
    offset     = 0


try:
  from bitarrayc import getbits as getbits_c, setbits as setbits_c
  getbits,setbits = getbits_c,setbits_c
except ImportError:
  getbits,setbits = getbits_python,setbits_python


###################################################################################################

# Test functions
def bitstring(n, minwidth=1, units=8):
  result = []
  while n:
    result.append(str(n%2))
    n >>= 1

  m = len(result)%units

  if not result or m:
    result.append( '0'*(units-m) )

  r = ''.join(reversed(result))

  while len(r) < minwidth:
    r = '0'*units + r

  return r


def test_bits(get, set, data, pos, value, width):
  # Convert original data to bitstring
  n = ''.join(bitstring(d) for d in data)

  # Verify getbits on the original data to be over-written
  assert bitstring(get(data, pos, width),minwidth=width)[-width:] == n[pos:pos+width]

  # Turn the value into a bit string
  v = bitstring(value,minwidth=width)[-width:]

  # Replace insert the value into the bitstring
  n = n[:pos] + v + n[pos+width:]

  # Perform the same operation on the data array
  set(data,pos,value,width)

  # Convert updated data to bitstring
  m = ''.join(bitstring(d) for d in data)

  if n != m:
    print data, pos, bitstring(value)[-width:], width
    print 'exp',len(n),n
    print 'obs',len(m),m

  # Verify that value that was written correctly
  assert bitstring(get(data, pos, width),minwidth=width)[-width:] == v

  # Verify the resulting bitstrings math
  assert n==m


def main():
  import time
  import array
  import numpy

  def mknumpy(bytes):
    return numpy.array(bytes,dtype='uint8')

  def mkarray(bytes):
    return array.array('B',bytes)

  imps = ('Python',getbits_python,setbits_python),('C',getbits_c,setbits_c)

  # Walk 1s/0s across a byte array with no/all bits set
  # Test over the following parameters:
  #    bytes: bit arrays from 1..8 bytes
  #    width: read/write 1..min(bytes*8,32) bit values
  #      pos: read/write 0..bytes*8-width bit position
  # The test cases cover:
  #   1) Setting 1111.. in a 0000... array
  #   2) Setting 0000.. in a 1111... array
  #   3) Setting 0101.. in a 0000... array
  #   4) Setting 0000.. in a 0101... array
  #   5) Setting 1010.. in a 0000... array
  #   6) Setting 0000.. in a 1010... array
  for mkrep in mkarray,mknumpy:
    for impname, get,set in imps:
      for bytes in range(1,9):
        for width in range(1,min(33,bytes*8)):
          for pos in range(bytes*8-width+1):
            test_bits(get,set,mkrep([0x00]*bytes),pos,0xFFFFFFFF,width)
            test_bits(get,set,mkrep([0xFF]*bytes),pos,0x00000000,width)
            test_bits(get,set,mkrep([0x00]*bytes),pos,0x55555555,width)
            test_bits(get,set,mkrep([0x55]*bytes),pos,0x00000000,width)
            test_bits(get,set,mkrep([0x00]*bytes),pos,0xAAAAAAAA,width)
            test_bits(get,set,mkrep([0xAA]*bytes),pos,0x00000000,width)

  for repname,mkrep in ('array',mkarray),('numpy',mknumpy):
    for impname,get,set in imps:
      t0 = time.time()

      for reps in range(100):
        for bytes in range(1,9):
          d=mkrep([0x00]*bytes)
          for width in range(1,min(33,bytes*8)):
            for pos in range(bytes*8-width+1):
              set(d,pos,0xFFFFFFFF,width)
              get(d,pos,width)

      t1 = time.time()
      print '%-5s, %-6s: time=%7.4f' % (repname,impname,t1-t0)


if __name__ == '__main__':
  main()
