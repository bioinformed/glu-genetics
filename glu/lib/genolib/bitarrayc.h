/*
File:          bitarray.h

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:

Abstract:      manipulator functions for packed bit arrays

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.

See GLU license for terms by running: glu license
*/

#include <Python.h>

/*
def getbits(data, i, readlen):
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
*/


extern unsigned long
_bitarray_getbits_slowpath(const unsigned char *data, ssize_t byte, ssize_t offset, ssize_t readlen, char **status);


static inline unsigned long
bitarray_getbits(const unsigned char *data, ssize_t buflen, ssize_t i, ssize_t readlen, char **status)
{
	ssize_t byte, offset;

	*status = NULL;

	if(i < 0)
		i += 8*buflen;

	if(i < 0)
	{
		*status = "read index out of range";
		return -1;
	}

	if(buflen <= (i+readlen-1)>>3)
	{
		*status = "read length out of range";
		return -1;
	}

	byte   = i>>3;
	offset = i&7;

	/* Fast path for common cases */
	if(readlen == 1)
		return (data[byte]>>(7-offset)) & 1;
	else if(readlen == 2 && offset < 7)
		return (data[byte]>>(6-offset)) & 3;
	else if(readlen == 4 && offset < 5)
		return (data[byte]>>(4-offset)) & 15;
	else if(readlen == 8 && offset == 0)
		return data[byte];

	return _bitarray_getbits_slowpath(data, byte, offset, readlen, status);
}


/*
def setbits(data, i, value, writelen):
  if i < 0:
    i += len(data)*8

  if i < 0:
    raise IndexError('write index out of range')

  if readlen <= 0 or readlen > 32:
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

  return data
*/

void
bitarray_setbits(unsigned char *data, ssize_t buflen, ssize_t i, unsigned long value, ssize_t writelen, char **status);
