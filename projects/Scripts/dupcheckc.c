#include <Python.h>

#if PY_VERSION_HEX < 0x02050000
typedef int Py_ssize_t;
#endif

/* Implements this function for various types of g1 and g2:

  def match(g1, g2):
    i = n = 0
    for a,b in izip(g1,g2):
      a = g1[x]
      b = g2[x]

      # If data are not missing
      if a and b:
        # Check for match
        if a is b:
          i += 1
        # Count total comparisons
        n += 1

    return i,n
*/

static inline long min(long a, long b)
{
  if(a<b)
    return a;
  return b;
}

static PyObject *
dupcheckc_match_object(PyObject *self, PyObject *args)
{
  PyObject   *g1, *g2, *a, *b;
  Py_ssize_t ilen, i, j, n;
  int        cmp;

  if(!PyArg_ParseTuple(args, "OO", &g1, &g2))
    return NULL;

  if(!PyList_Check(g1))
  {
    PyErr_SetString(PyExc_TypeError, "g1 argument must be a list");
    return NULL;
  }

  if(!PyList_Check(g2))
  {
    PyErr_SetString(PyExc_TypeError, "g2 argument must be a list");
    return NULL;
  }

  ilen = min(PyList_Size(g1),PyList_Size(g2));
  if(ilen < 0)
    return NULL;

  a = b = NULL;

  i = n = 0;
  for(j = 0; j < ilen; ++j)
  {
    a = PyList_GetItem(g1, j);
    b = PyList_GetItem(g2, j);

    if(!a || !b)
      return NULL;

    if(a != Py_None && b != Py_None)
    {
      cmp = PyObject_RichCompareBool(a, b, Py_EQ);

      if(cmp < 0)
        return NULL;
      else if(cmp == 1)
        i += 1;

      n += 1;
    }
  }

  return Py_BuildValue("(ii)", i, n);
}

static PyObject *
dupcheckc_match_byte_array(PyObject *self, PyObject *args)
{
  PyObject   *g1, *g2;
  Py_ssize_t blen1, blen2;
  Py_ssize_t ilen, i, j, n;
  const unsigned char *b1, *b2;
  unsigned char        a1, a2;

  if(!PyArg_ParseTuple(args, "OO", &g1, &g2))
    return NULL;

  if(PyObject_AsReadBuffer(g1, (const void**)&b1, &blen1) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g1 argument must support the buffer interface");
    return NULL;
  }

  if(PyObject_AsReadBuffer(g2, (const void**)&b2, &blen2) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g2 argument must support the buffer interface");
    return NULL;
  }

  ilen = min(blen1,blen2)/sizeof(char);

  if(ilen < 0)
    return NULL;

  i = n = 0;
  for(j = 0; j < ilen; ++j)
  {
    a1 = b1[j];
    a2 = b2[j];

    if(a1 && a2)
    {
      if(a1 == a2)
        i += 1;
      n += 1;
    }
  }

  return Py_BuildValue("(ii)", i, n);
}

static PyObject *
dupcheckc_match_short_array(PyObject *self, PyObject *args)
{
  PyObject   *g1, *g2;
  Py_ssize_t blen1, blen2;
  Py_ssize_t ilen, i, j, n;
  const unsigned short *b1, *b2;
  unsigned short        a1, a2;

  if(!PyArg_ParseTuple(args, "OO", &g1, &g2))
    return NULL;

  if(PyObject_AsReadBuffer(g1, (const void**)&b1, &blen1) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g1 argument must support the buffer interface");
    return NULL;
  }

  if(PyObject_AsReadBuffer(g2, (const void**)&b2, &blen2) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g2 argument must support the buffer interface");
    return NULL;
  }

  ilen = min(blen1,blen2)/sizeof(short);

  if(ilen < 0)
    return NULL;

  i = n = 0;
  for(j = 0; j < ilen; ++j)
  {
    a1 = b1[j];
    a2 = b2[j];

    if(a1 && a2)
    {
      if(a1 == a2)
        i += 1;
      n += 1;
    }
  }

  return Py_BuildValue("(ii)", i, n);
}


static PyObject *
dupcheckc_match_long_array(PyObject *self, PyObject *args)
{
  PyObject   *g1, *g2;
  Py_ssize_t blen1, blen2;
  Py_ssize_t ilen, i, j, n;
  const unsigned long *b1, *b2;
  unsigned long        a1, a2;

  if(!PyArg_ParseTuple(args, "OOd", &g1, &g2))
    return NULL;

  if(PyObject_AsReadBuffer(g1, (const void**)&b1, &blen1) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g1 argument must support the buffer interface");
    return NULL;
  }

  if(PyObject_AsReadBuffer(g2, (const void**)&b2, &blen2) < 0)
  {
    PyErr_SetString(PyExc_TypeError, "g2 argument must support the buffer interface");
    return NULL;
  }

  ilen = min(blen1,blen2)/sizeof(long);

  if(ilen < 0)
    return NULL;

  i = n = 0;
  for(j = 0; j < ilen; ++j)
  {
    a1 = b1[j];
    a2 = b2[j];

    if(a1 && a2)
    {
      if(a1 == a2)
        i += 1;
      n += 1;
    }
  }

  return Py_BuildValue("(ii)", i, n);
}


static PyMethodDef dupcheckc_methods[] = {
        {"match_object",      dupcheckc_match_object,      METH_VARARGS, "Match two lists of genotypes objects"},
        {"match_byte_array",  dupcheckc_match_byte_array,  METH_VARARGS, "Match two byte arrays of genotypes"},
        {"match_short_array", dupcheckc_match_short_array, METH_VARARGS, "Match two short int arrays of genotypes"},
        {"match_long_array",  dupcheckc_match_long_array,  METH_VARARGS, "Match two long int arrays of genotypes"},
        {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initdupcheckc(void)
{
  (void) Py_InitModule("dupcheckc", dupcheckc_methods);
}
