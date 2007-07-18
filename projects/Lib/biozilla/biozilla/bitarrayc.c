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

unsigned long
bitarray_getbits(const unsigned char *data, ssize_t buflen, ssize_t i, ssize_t readlen, char **status)
{
	ssize_t byte, offset, n, tail;
	unsigned long result;
	unsigned char b;

	*status = NULL;

	if(i < 0)
		i += 8*buflen;

	if(i < 0)
	{
		*status = "read index out of range";
		return -1;
	}

	if(readlen <= 0 || readlen > 32)
	{
		*status = "read length must be between 1..32";
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

	/* Slow path */
	result = 0;
	while(readlen)
	{
		b = data[byte];
		n = 8 - offset;

		/* mask high order bits */
		if(offset)
			b &= (1<<n)-1;

		/* shift off low order bits (the tail) */
		tail = n - readlen;
		if(tail > 0)
		{
			b >>= tail;
			n  -= tail;
		}

		result   = (result<<n)|b;
		readlen -= n;
		byte    += 1;
		offset   = 0;
	}

	return result;
}

static PyObject *
py_getbits(PyObject *self, PyObject *args)
{
	PyObject   *data;
	Py_ssize_t i, readlen, buflen;
	unsigned long result;
	unsigned char *dbuff;
	char *status;

	if(!PyArg_ParseTuple(args, "Onn:getbits", &data, &i, &readlen))
		return NULL;

	if(PyObject_AsReadBuffer(data, (const void**)&dbuff, &buflen) < 0)
	{
		PyErr_SetString(PyExc_TypeError, "data argument must support the buffer interface");
		return NULL;
	}

	result = bitarray_getbits(dbuff, buflen, i, readlen, &status);

	if(status)
	{
		PyErr_SetString(PyExc_IndexError, status);
		return NULL;
	}

	return PyInt_FromLong(result);
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
bitarray_setbits(unsigned char *data, ssize_t buflen, ssize_t i, unsigned long value, ssize_t writelen, char **status)
{
	ssize_t byte, offset, n, tail;
	unsigned long mask;

	*status = NULL;

	if(i < 0)
		i += 8*buflen;

	if(i < 0)
	{
		*status = "write index out of range";
		return;
	}

	if(writelen <= 0 || writelen > 32)
	{
		*status = "write length must be between 1..32";
		return;
	}

	if(buflen <= (i+writelen-1)>>3)
	{
		*status = "write length out of range";
		return;
	}

	byte   = i>>3;
	offset = i&7;

	/* Fast path for common cases */
	if(writelen == 1)
	{
		data[byte] = (data[byte]&~( 1<<(7-offset))) | ((value& 1)<<(7-offset));
		return;
	}
	else if(writelen == 2 && offset < 7)
	{
		data[byte] = (data[byte]&~( 3<<(6-offset))) | ((value& 3)<<(6-offset));
		return;
	}
	else if(writelen == 4 && offset < 5)
	{
		data[byte] = (data[byte]&~(15<<(4-offset))) | ((value&15)<<(4-offset));
		return;
	}
	else if(writelen == 8 && offset == 0)
	{
		data[byte] = value&0xFF;
		return;
	}

	/* Slow path */
	while(writelen)
	{
		n          = 8 - offset;
		tail       = n-writelen;
		if(tail<0)
			tail = 0;
		n         -= tail;
		mask       = ((1<<n)-1)<<tail;
		writelen  -= n;
		/* Mask (xor) the old values, then set (or) the new bits */
		data[byte] = (data[byte]&(~mask)) | ((value>>writelen<<tail)&mask);
		byte      += 1;
		offset     = 0;
	}
}

static PyObject *
py_setbits(PyObject *self, PyObject *args)
{
	PyObject   *data;
	Py_ssize_t i, writelen, buflen;
	unsigned long value;
	unsigned char *dbuff;
	char *status;

	if(!PyArg_ParseTuple(args, "Onnn:setbits", &data, &i, &value, &writelen))
		return NULL;

	if(PyObject_AsWriteBuffer(data, (void**)&dbuff, &buflen) < 0)
	{
		PyErr_SetString(PyExc_TypeError, "data argument must support the buffer interface");
		return NULL;
	}

	bitarray_setbits(dbuff, buflen, i, value, writelen, &status);

	if(status)
	{
		PyErr_SetString(PyExc_IndexError, status);
		return NULL;
	}

	Py_INCREF(Py_None);
	return Py_None;
}


static PyMethodDef bitarrayc_methods[] = {
	{"getbits",      py_getbits,      METH_VARARGS, "Get bits"},
	{"setbits",      py_setbits,      METH_VARARGS, "Set bits"},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initbitarrayc(void)
{
	(void) Py_InitModule("bitarrayc", bitarrayc_methods);
}
