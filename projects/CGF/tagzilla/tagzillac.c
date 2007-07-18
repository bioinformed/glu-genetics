/*
File:          tagzillac.c

Authors:       Kevin Jacobs (jacobs@bioinformed.com)

Created:       November 8, 2005

Abstract:      A module that contains an accelerated C implementation of the
               functions in tagzilla.py.  This module is optional, since
               the tagzilla will fall-back to the slower, but pure Python,
               function if necessary.

Compatibility: CPython 2.4 and above

Requires:      No external dependencies, yet...

Version:       0.5

Revision:      $Id: $

Copyright (c) 2005 BioInformed Consulting Services.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.
*/

#include <Python.h>
#include <math.h>


static inline double dmax(double a, double b)
{
  if(a > b)
    return a;
  else
    return b;
}

static inline double dmin(double a, double b)
{
  if(a < b)
    return a;
  else
    return b;
}

static inline long max(long a, long b)
{
  if(a > b)
    return a;
  else
    return b;
}

static inline long min(long a, long b)
{
  if(a < b)
    return a;
  else
    return b;
}

static inline int heq(char *g, char a)
{
  return (g[0] == a || g[0] == ' ') && (g[1] == a || g[1] == ' ');
}

static inline int gmissing(char *g)
{
  return g[0] == ' ' && g[1] == ' ';
}

static inline int homozyg(char *g)
{
  return g[0] == g[1];
}

static inline void gswap(int *x, int *y)
{
  int tmp = *x;
  *x = *y;
  *y = tmp;
}


static int find_het_slowpath(char *a, char *b, char* g)
{
  /**** Slow path: Detect alleles ****/
  int g1, g2;

  /* Detect first allele */
  if(!*a && g[0] != ' ')
    *a = g[0];
  if(!*a && g[1] != ' ')
    *a = g[1];

  /* Detect second allele */
  if(*a && !*b && g[0] != ' ' && g[0] != *a)
    *b = g[0];
  if(*a && !*b && g[1] != ' ' && g[1] != *a)
    *b = g[1];

  /* Redetect more than 2 alleles */
  g1 = g[0] == ' ' || g[0] == *a || g[0] == *b;
  g2 = g[1] == ' ' || g[1] == *a || g[1] == *b;

  return g1 && g2;
}


static inline int find_het(char *a, char *b, char* g)
{
  /* Match genotype to missing or known alleles */
  int g1 = g[0] == ' ' || g[0] == *a || g[0] == *b;
  int g2 = g[1] == ' ' || g[1] == *a || g[1] == *b;

  /* Fast path -- 2 known alleles */
  if( g1 && g2 )
    return 1;

  /* Fast path -- >2 alleles */
  if( *a && *b )
    return 0;

  /**** Slow path: Detect alleles ****/
  return find_het_slowpath(a,b,g);
}


static PyObject *
count_haplotypes(PyObject *self, PyObject *args)
{
  PyObject *locus1, *locus2, *haplo_counts, *gc1, *gc2;
  char     *g1, *g2;
  char     a1,a2;
  char     b1,b2;
  long     i, a, b, n;
  int      x, mult;
  int      diplos[9];
  int      haplos[5];

  if(!PyArg_ParseTuple(args, "OO", &locus1, &locus2))
    return NULL;

  if(!PyList_Check(locus1))
  {
    PyErr_SetString(PyExc_TypeError, "locus1 argument must be a list");
    return NULL;
  }

  if(!PyList_Check(locus2))
  {
    PyErr_SetString(PyExc_TypeError, "locus2 argument must be a list");
    return NULL;
  }

  n = PyList_Size(locus1);

  if(n < 0)
    return NULL;

  if( n != PyList_Size(locus2) )
  {
    PyErr_SetString(PyExc_ValueError, "locus1 and locus2 must be the same length");
    return NULL;
  }

  a1 = a2 = '\0';
  b1 = b2 = '\0';

  for(i = 0; i<9; ++i)
    diplos[i] = 0;

  for(i = 0; i < n; ++i)
  {
    gc1 = PyList_GET_ITEM(locus1, i);
    gc2 = PyList_GET_ITEM(locus2, i);

    if( !PyString_CheckExact(gc1) || !PyString_CheckExact(gc2)
      || PyString_GET_SIZE(gc1) != 2 || PyString_GET_SIZE(gc2) != 2)
    {
      PyErr_SetString(PyExc_ValueError, "invalid genotype: genotypes must be strings of length 2");
      return NULL;
    }

    g1 = PyString_AS_STRING(gc1);
    g2 = PyString_AS_STRING(gc2);

    if(gmissing(g1) || gmissing(g2))
      continue;

    x = find_het(&a1, &a2, g1) + find_het(&b1, &b2, g2);

    if(x != 2)
    {
      PyErr_SetString(PyExc_ValueError, "invalid genotypes: loci may have no more than 2 alleles");
      return NULL;
    }

    if     (heq(g1,a1)) a = 0;
    else if(heq(g1,a2)) a = 2;
    else                a = 1;

    if     (heq(g2,b1)) b = 0;
    else if(heq(g2,b2)) b = 2;
    else                b = 1;

    mult = 1;
    if(homozyg(g1) && homozyg(g2))
      mult = 2;

    diplos[3*a + b] += mult;
  }

  /* Haplotype count matrix elements

                    haplotypes
                 0   1   2   3   4
                AC  AD  BC  BD  Dbl.Het. (AC,BD or AD,BC)
      0: AA CC [ 2,  0,  0,  0, 0 ]
   d  0: A* C* [ 1,  0,  0,  0, 0 ]
   i  1: AA CD [ 1,  1,  0,  0, 0 ]
   p  2: AA DD [ 0,  2,  0,  0, 0 ]
   l  2: A* D* [ 0,  1,  0,  0, 0 ]
   o  3: AB CC [ 1,  0,  1,  0, 0 ]
   t  4: AB CD [ 0,  0,  0,  0, 1 ]
   y  5: AB DD [ 0,  1,  0,  1, 0 ]
   p  6: BB CC [ 0,  0,  2,  0, 0 ]
   e  6: B* C* [ 0,  0,  1,  0, 0 ]
   s  7: BB CD [ 0,  0,  1,  1, 0 ]
      8: BB DD [ 0,  0,  0,  2, 0 ]
      8: B* D* [ 0,  0,  0,  1, 0 ]

      NOTE: Multiplicities (1 for hemizygous vs. 2 for double homozygous)
            are taken into account by the diplotype counting.
  */

  haplos[0] = diplos[0] + diplos[1] + diplos[3];
  haplos[1] = diplos[2] + diplos[1] + diplos[5];
  haplos[2] = diplos[6] + diplos[3] + diplos[7];
  haplos[3] = diplos[8] + diplos[5] + diplos[7];
  haplos[4] = diplos[4];

  if(a2 && a1 > a2)
  {
    gswap(&haplos[0],&haplos[2]);
    gswap(&haplos[1],&haplos[3]);
  }

  if(b2 && b1 > b2)
  {
    gswap(&haplos[0],&haplos[1]);
    gswap(&haplos[2],&haplos[3]);
  }

  haplo_counts = PyTuple_New(5);

  for(i = 0; i < 5; ++i)
    PyTuple_SET_ITEM(haplo_counts, i, PyInt_FromLong(haplos[i]) );

  return haplo_counts;
}

static PyObject *
estimate_ld(PyObject *self, PyObject *args)
{
  const double TOLERANCE = 10e-7;
  int    i;
  long   c11, c12, c21, c22, dh, n;
  double p, q, old_p11, p11, p12, p21, p22, a, nx1, nx2;
  double d, d_max, dprime, r2;
  PyObject *results;

  if(!PyArg_ParseTuple(args, "lllll", &c11, &c12, &c21, &c22, &dh))
    return NULL;

  if(!dh && (c11+c12 == 0 || c21+c22 == 0 || c11+c21 == 0 || c12+c22 == 0))
  {
    dprime = r2 = 0;
    goto bail;
  }

  /* Initial estimates */
  n = c11 + c12 + c21 + c22 + 2*dh;
  p = ((double)(c11 + c12 + dh))/n;
  q = ((double)(c11 + c21 + dh))/n;

  p11 = p*q;
  p12 = p*(1-q);
  p21 = (1-p)*q;
  p22 = (1-p)*(1-q);

  for(i = 0; i < 100; ++i)
  {
    old_p11 = p11;

    /* Force estimates away from boundaries */
    p11=dmax(10e-10, p11);
    p12=dmax(10e-10, p12);
    p21=dmax(10e-10, p21);
    p22=dmax(10e-10, p22);

    a = p11*p22 + p12*p21;

    nx1 = dh*p11*p22/a;
    nx2 = dh*p12*p21/a;

    p11 = (c11+nx1)/n;
    p12 = (c12+nx2)/n;
    p21 = (c21+nx2)/n;
    p22 = (c22+nx1)/n;

    if(fabs(old_p11-p11) < TOLERANCE)
      break;
  }

  d = p11*p22 - p12*p21;

  if(d > 0)
    d_max = dmin( p*(1-q), (1-p)*q );
  else
    d_max = -dmin( p*q, (1-p)*(1-q) );

  dprime = d/d_max;
  r2 = d*d/(p*(1-p)*q*(1-q));

bail:
  results = PyTuple_New(2);

  if(!results)
    return NULL;

  PyTuple_SET_ITEM(results, 0, PyFloat_FromDouble(r2));
  PyTuple_SET_ITEM(results, 1, PyFloat_FromDouble(dprime));

  return results;
}

static PyObject *
estimate_maf(PyObject *self, PyObject *args)
{
  PyObject *genos, *gc;
  char     *g;
  char     a1,a2;
  int      x;
  long     i, c1, c2;
  double   maf;

  if(!PyArg_ParseTuple(args, "O", &genos))
    return NULL;

  if(!PyList_Check(genos))
  {
    PyErr_SetString(PyExc_TypeError, "genos argument must be a list");
    return NULL;
  }

  long n = PyList_Size(genos);

  if(n < 0)
    return NULL;

  a1 = a2 = '\0';
  c1 = c2 = 0;

  for(i = 0; i < n; ++i)
  {
    gc = PyList_GET_ITEM(genos, i);

    if( !PyString_CheckExact(gc)|| PyString_GET_SIZE(gc) != 2 )
    {
      PyErr_SetString(PyExc_TypeError, "invalid genotype: genotypes must be strings of length 2");
      return NULL;
    }

    g = PyString_AS_STRING(gc);

    if(gmissing(g))
      continue;

    x = find_het(&a1, &a2, g);

    if(!x)
    {
      PyErr_SetString(PyExc_ValueError, "invalid genotypes: locus may have no more than 2 alleles");
      return NULL;
    }

    if(g[0] == a1) c1 += 1;
    if(g[0] == a2) c2 += 1;
    if(g[1] == a1) c1 += 1;
    if(g[1] == a2) c2 += 1;
  }

  if(c1 == 0 || c2 == 0)
    maf = 0;
  else
    maf = (dmin(c1,c2))/(c1+c2);

  return PyFloat_FromDouble(maf);
}


static PyMethodDef tagzillac_methods[] = {
        {"count_haplotypes", count_haplotypes, METH_VARARGS, "Count haplotypes at two loci"},
        {"estimate_ld",      estimate_ld,      METH_VARARGS, "Compute LD statics from haplotype counts"},
        {"estimate_maf",     estimate_maf,     METH_VARARGS, "Estimate minor allele frequency (MAF)"},
        {NULL, NULL, 0, NULL}   /* Sentinel */
};

PyMODINIT_FUNC
inittagzillac(void)
{
    (void) Py_InitModule("tagzillac", tagzillac_methods);
}
