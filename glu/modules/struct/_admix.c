/*
File:          _admix.c

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:       Mon Dec 14 09:40:27 EST 2009

Abstract:      GLU struct.admix accelerator module

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.
See GLU license for terms by running: glu license
*/

#include <Python.h>
#include <math.h>
#include <float.h>
#include "numpy/arrayobject.h"

#define Array_CheckType(op)    (PyArray_CheckExact(op)              && \
                                PyArray_CHKFLAGS(op,NPY_CARRAY_RO)  && \
                                PyArray_ISONESEGMENT(op)            && \
                                PyArray_TYPE(op)==PyArray_DOUBLE)


PyObject *
admixture_log_likelihood(PyObject *self, PyObject *args)
{
	PyObject *frequencies, *mixture;
	double *f, *x;
	Py_ssize_t n, k, m, i, j;
	double logl, lpart, lpart_try;
	double log_DBL_MAX = log(DBL_MAX);

	if(!PyArg_ParseTuple(args, "OO", &frequencies, &mixture))
		return NULL;

	if(!Array_CheckType(frequencies) || PyArray_NDIM(frequencies)!=2)
	{
		PyErr_SetString(PyExc_TypeError, "frequencies must be a 2d float ndarray");
		return NULL;
	}

	if(!Array_CheckType(mixture) || PyArray_NDIM(mixture)!=1)
	{
		PyErr_SetString(PyExc_TypeError, "mixture parameters must be a 1d float ndarray");
		return NULL;
	}

	n = PyArray_DIMS(frequencies)[0];
	k = PyArray_DIMS(frequencies)[1];
	m = PyArray_DIMS(mixture)[0];

	if(m != k)
	{
		PyErr_Format(PyExc_ValueError, "unexpected number of mixture parameters %zd."
		                               "  Expected %zd.", m, k);
		return NULL;
	}

	/* The naive log likelihood could be accumulated by taking the log of the likelihood of
	 * each locus.  Although straightforward, this approach is also extremely slow due to
	 * the cost of computing the log function.  Instead, we can accumulate log likelihood in
	 * parts, only taking logs of the likelihood when the current partial likelihood becomes
	 * too small.
	 *
	 * The likelihood starts at zero and the partial (non-log) likelihood begins with a
	 * positive bias of DBL_MAX (~1.79e308).  The likelihood is accumulated for each locus
	 * from this large value until it decreases to under DBL_MIN (~1.79e-308), where
	 * precision starts to be lost.  The previous value of the partial likelihood is then
	 * converted to the log-scale, the bias is subtracted (-log(DBL_MAX)), and it is
	 * accumulated.  The partial likelihood is re-initialized as DBL_MAX times the current
	 * locus.
	 *
	 * The following code is also manually unrolled for various common values of k.  This
	 * results in > 2x speedups by avoiding several loops, or else it wouldn't be worth the
	 * bother.
	 */

	logl  = 0;
	lpart = DBL_MAX;

	f = (double *)PyArray_DATA(frequencies);

	if(k==2)
	{
		x = (double *)PyArray_DATA(mixture);

		for(i=0; i<n; i++,f+=k)
		{
			double ll = f[0]*x[0] + f[1]*x[1];

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}
	else if(k==3)
	{
		x = (double *)PyArray_DATA(mixture);

		for(i=0; i<n; i++,f+=k)
		{
			double ll = f[0]*x[0] + f[1]*x[1] + f[2]*x[2];

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}
	else if(k==4)
	{
		x = (double *)PyArray_DATA(mixture);

		for(i=0; i<n; i++,f+=k)
		{
			double ll = f[0]*x[0] + f[1]*x[1] + f[2]*x[2] + f[3]*x[3];

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}
	else if(k==5)
	{
		x = (double *)PyArray_DATA(mixture);

		for(i=0; i<n; i++,f+=k)
		{
			double ll = f[0]*x[0] + f[1]*x[1] + f[2]*x[2] + f[3]*x[3] + f[4]*x[4];

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}
	else if(k==6)
	{
		x = (double *)PyArray_DATA(mixture);

		for(i=0; i<n; i++,f+=k)
		{
			double ll = f[0]*x[0] + f[1]*x[1] + f[2]*x[2]
			          + f[3]*x[3] + f[4]*x[4] + f[5]*x[5];

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			double ll=0;

			/* Compute the likelihood for the current locus */
			x = (double *)PyArray_DATA(mixture);
			for(j=0; j<k; ++j,++f,++x)
				ll += (*f) * (*x);

	                lpart_try = lpart*ll;

	                if(lpart_try>DBL_MIN)
			{
				/* If not, accept the value */
				lpart = lpart_try;
			}
			else
	                {
				logl += log(lpart)-log_DBL_MAX;
				lpart = DBL_MAX*ll;
			}
		}
	}

	/* Account for the final partial likelihood */
	logl += log(lpart) - log_DBL_MAX;

	return PyFloat_FromDouble(logl);
}


PyObject *
admixture_log_likelihood_derivative(PyObject *self, PyObject *args)
{
	PyObject *frequencies, *mixture, *derivative;
	double *f, *x, *d;
	Py_ssize_t n, k, m, i, j;

	if(!PyArg_ParseTuple(args, "OO", &frequencies, &mixture))
		return NULL;

	if(!Array_CheckType(frequencies) || PyArray_NDIM(frequencies)!=2)
	{
		PyErr_SetString(PyExc_TypeError, "frequencies must be a 2d float ndarray");
		return NULL;
	}

	if(!Array_CheckType(mixture) || PyArray_NDIM(mixture)!=1)
	{
		PyErr_SetString(PyExc_TypeError, "mixture parameters must be a 1d float ndarray");
		return NULL;
	}

	n = PyArray_DIMS(frequencies)[0];
	k = PyArray_DIMS(frequencies)[1];
	m = PyArray_DIMS(mixture)[0];

	if(m != k)
	{
		PyErr_Format(PyExc_ValueError, "unexpected number of mixture parameters %zd.  Expected %zd",
		                                m, k);
		return NULL;
	}

	derivative = PyArray_ZEROS(1, &k, NPY_DOUBLE, 0);
	if(!derivative) return NULL;

	f = (double *)PyArray_DATA(frequencies);

	/* The following code is also manually unrolled for various common values of k.  This
	 * results in > 2x speedups by avoiding several loops, or else it wouldn't be worth the
	 * bother.
	 */
	if(k==2)
	{
                x = (double *)PyArray_DATA(mixture);
                d = (double *)PyArray_DATA(derivative);

	        for(i=0; i<n; ++i,f+=k)
	        {
	                /* Pass 1: Compute denominator */
			double u = 1/(f[0]*x[0] + f[1]*x[1]);

	                /* Pass 2: Accumulate derivative */
			d[0] += f[0]*u;
			d[1] += f[1]*u;
	        }
	}
	else if(k==3)
	{
                x = (double *)PyArray_DATA(mixture);
                d = (double *)PyArray_DATA(derivative);

	        for(i=0; i<n; ++i,f+=k)
	        {
	                /* Pass 1: Compute denominator */
			double u = 1/(f[0]*x[0] + f[1]*x[1] + f[2]*x[2]);

	                /* Pass 2: Accumulate derivative */
			d[0] += f[0]*u;
			d[1] += f[1]*u;
			d[2] += f[2]*u;
	        }
	}
	else if(k==4)
	{
                x = (double *)PyArray_DATA(mixture);
                d = (double *)PyArray_DATA(derivative);

	        for(i=0; i<n; ++i,f+=k)
	        {
	                /* Pass 1: Compute denominator */
			double u = 1/(f[0]*x[0] + f[1]*x[1] + f[2]*x[2] + f[3]*x[3]);

	                /* Pass 2: Accumulate derivative */
			d[0] += f[0]*u;
			d[1] += f[1]*u;
			d[2] += f[2]*u;
			d[3] += f[3]*u;
	        }
	}
	else if(k==5)
	{
                x = (double *)PyArray_DATA(mixture);
                d = (double *)PyArray_DATA(derivative);

	        for(i=0; i<n; ++i,f+=k)
	        {
	                /* Pass 1: Compute denominator */
			double u = 1/(f[0]*x[0] + f[1]*x[1] + f[2]*x[2] + f[3]*x[3] + f[4]*x[4]);

	                /* Pass 2: Accumulate derivative */
			d[0] += f[0]*u;
			d[1] += f[1]*u;
			d[2] += f[2]*u;
			d[3] += f[3]*u;
			d[4] += f[4]*u;
	        }
	}
	else if(k==6)
	{
                x = (double *)PyArray_DATA(mixture);
                d = (double *)PyArray_DATA(derivative);

	        for(i=0; i<n; ++i,f+=k)
	        {
	                /* Pass 1: Compute denominator */
			double u = 1/(f[0]*x[0] + f[1]*x[1] + f[2]*x[2]
			            + f[3]*x[3] + f[4]*x[4] + f[5]*x[5]);

	                /* Pass 2: Accumulate derivative */
			d[0] += f[0]*u;
			d[1] += f[1]*u;
			d[2] += f[2]*u;
			d[3] += f[3]*u;
			d[4] += f[4]*u;
			d[5] += f[5]*u;
	        }
	}
	else
	{
	        for(i=0; i<n; ++i)
	        {
	                double  u  = 0;
	                double *ff = f;

	                /* Pass 1: Compute denominator */
	                x = (double *)PyArray_DATA(mixture);
	                for(j=0; j<k; ++j,++f,++x)
	                        u += (*f) * (*x);

	                /* Pass 2: Accumulate derivative */
	                d = (double *)PyArray_DATA(derivative);
                        u = 1/u;
	                for(j=0,f=ff; j<k; ++j,++f,++d)
	                        *d += (*f)*u;
	        }
	}
	return derivative;
}


PyMODINIT_FUNC
init_admix(void)
{
	PyObject *m;

	static PyMethodDef admixmodule_methods[] = {
	       {"admixture_log_likelihood", (PyCFunction)admixture_log_likelihood, METH_VARARGS,
	        "Compute admixture likelihood for given population frequencies and admixture estimates"},
	       {"admixture_log_likelihood_derivative", (PyCFunction)admixture_log_likelihood_derivative, METH_VARARGS,
	        "Compute the derivative of the admixture likelihood for given population frequencies and admixture estimates"},
	       {NULL}  /* Sentinel */
	};

	import_array();

	m = Py_InitModule3("_admix", admixmodule_methods, "GLU struct.admix accelerator module");
}
