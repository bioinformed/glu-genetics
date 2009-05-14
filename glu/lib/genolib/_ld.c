/*
File:          _ld.c

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:

Abstract:      fast LD calculation for biallelic loci

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.
See GLU license for terms by running: glu license
*/

#include <Python.h>
#include <math.h>
#include "_genoarray.h"


static inline double
dmax(double a, double b)
{
	if(a > b)
		return a;
	else
		return b;
}

static inline double
dmin(double a, double b)
{
	if(a < b)
		return a;
	else
		return b;
}

static inline int
homozygous_for(GenotypeObject *g, int a)
{
	return g->allele1_index == a && g->allele2_index == a;
}

static inline int
homozygous(GenotypeObject *g)
{
	return g->allele1_index == g->allele2_index;
}

static inline void
gswap(Py_ssize_t *x, Py_ssize_t *y)
{
	Py_ssize_t tmp = *x;
	*x = *y;
	*y = tmp;
}

static int
find_het_slowpath(Py_ssize_t *a, Py_ssize_t *b, GenotypeObject *g)
{
	/**** Slow path: Detect alleles ****/
	Py_ssize_t g1, g2;

	/* Detect first allele */
	if(*a==-1 && g->allele1_index)
		*a = g->allele1_index;
	if(*a==-1 && g->allele2_index)
		*a = g->allele2_index;

	/* Detect second allele */
	if(*a!=-1 && *b==-1 && g->allele1_index && g->allele1_index != *a)
		*b = g->allele1_index;
	if(*a!=-1 && *b==-1 && g->allele2_index && g->allele2_index != *a)
		*b = g->allele2_index;

	/* Redetect more than 2 alleles */
	g1 = !g->allele1_index || g->allele1_index == *a || g->allele1_index == *b;
	g2 = !g->allele2_index || g->allele2_index == *a || g->allele2_index == *b;

	return g1 && g2;
}

static inline int
find_het(Py_ssize_t *a, Py_ssize_t *b, GenotypeObject *g)
{
	/* Match genotype to missing or known alleles */
	int g1 = !g->allele1_index || g->allele1_index == *a || g->allele1_index == *b;
	int g2 = !g->allele2_index || g->allele2_index == *a || g->allele2_index == *b;

	/* Fast path -- 2 known alleles */
	if( g1 && g2 )
		return 1;

	/* Fast path -- >2 alleles */
	if( *a!=-1 && *b!=-1 )
		return 0;

	/**** Slow path: Detect alleles ****/
	return find_het_slowpath(a,b,g);
}

PyObject *diplos_to_haplos(Py_ssize_t diplos[9])
{
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

	PyObject *haplo_counts;
	Py_ssize_t i, haplos[5];

	haplos[0] = diplos[0] + diplos[1] + diplos[3];
	haplos[1] = diplos[2] + diplos[1] + diplos[5];
	haplos[2] = diplos[6] + diplos[3] + diplos[7];
	haplos[3] = diplos[8] + diplos[5] + diplos[7];
	haplos[4] = diplos[4];

	haplo_counts = PyTuple_New(5);
	if(!haplo_counts) return NULL;

	for(i = 0; i < 5; ++i)
		PyTuple_SET_ITEM(haplo_counts, i, PyInt_FromSsize_t(haplos[i]) );

	return haplo_counts;
}

static int*
genotype_categories(UnphasedMarkerModelObject *model)
{
	int *gcat=NULL;
	Py_ssize_t glen, homoz, i;

	if(!model || !UnphasedMarkerModel_Check(model) || !model->genotypes)
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return NULL;
	}

	glen = PyList_Size(model->genotypes);
	if(glen < 0) return NULL;

	gcat = malloc( sizeof(int)*glen );
	if(!gcat) goto error;

        homoz = 0;

	for(i=0; i<glen; ++i)
	{
		GenotypeObject *g = (GenotypeObject *)PyList_GetItem(model->genotypes, i);

		if(!g || !Genotype_CheckExact(g))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object at index %zd", i);
			goto error;
		}

		if(!g->allele1_index && !g->allele2_index)
		        gcat[i] = -1;
		else if(!g->allele1_index || !g->allele2_index)
		{
			PyErr_SetString(PyExc_ValueError, "Hemizygote LD estimation is not currently supported");
			goto error;
		}
		else if(g->allele1_index != g->allele2_index)
			gcat[i] = 1;
		else if(homoz > 2)
		{
			PyErr_SetString(PyExc_ValueError, "invalid genotypes: loci may have no more than 2 alleles");
			goto error;
		}
		else
		{
			gcat[i] = homoz;
			homoz += 2;
		}
	}
	return gcat;

error:
	if(gcat) free(gcat);
	return NULL;
}


static inline void
count_2bit(Py_ssize_t *diplos, int cat1, int cat2)
{
	if(cat1>=0 && cat2>=0)
		diplos[3*cat1+cat2] += (cat1!=1 && cat2!=1)? 2 : 1;
}


static PyObject *
count_haplotypes_2bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	UnphasedMarkerModelObject *model1=NULL, *model2=NULL;
	int *gcat1=NULL, *gcat2=NULL;
	const unsigned char *g1, *g2;
	const unsigned int *offsets;
	Py_ssize_t len1, len2, i, diplos[9];

	if(!PyArg_ParseTuple(args, "OO", &genos1, &genos2))
		return NULL;

	if(!GenotypeArray_CheckExact(genos1))
	{
		PyErr_SetString(PyExc_TypeError,"genos1 must be a GenotypeArray instance");
		return NULL;
	}

	if(!GenotypeArray_CheckExact(genos2))
	{
		PyErr_SetString(PyExc_TypeError,"genos2 must be a GenotypeArray instance");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 2 || genos2->descriptor->homogeneous != 2)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 2 bit width");
		return NULL;
	}

	len1 = PyList_Size(genos1->descriptor->models);
	if(len1==-1) goto error;
	len2 = PyList_Size(genos2->descriptor->models);
	if(len2==-1) goto error;

	if(len1 != len2)
	{
		PyErr_Format(GenotypeRepresentationError,"genotype array sizes do not match: %zd != %zd",
		             len1, len2);
		goto error;
	}

	for(i = 0; i<9; ++i)
		diplos[i] = 0;

	if(!len1) goto done;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	g1 = genos1->data;
	g2 = genos2->data;

	model1 = (UnphasedMarkerModelObject *)PyList_GetItem(genos1->descriptor->models, 0);
	if(!model1) goto error;
	model2 = (UnphasedMarkerModelObject *)PyList_GetItem(genos2->descriptor->models, 0);
	if(!model2) goto error;

	gcat1 = genotype_categories(model1);
	if(!gcat1) goto error;
	gcat2 = genotype_categories(model2);
	if(!gcat2) goto error;

	for(i=0; i<len1/4; ++i)
	{
		unsigned char a = g1[i];
		unsigned char b = g2[i];
		count_2bit(diplos, gcat1[a&0x03], gcat2[b&0x03]);
		a>>=2; b>>=2;
		count_2bit(diplos, gcat1[a&0x03], gcat2[b&0x03]);
		a>>=2; b>>=2;
		count_2bit(diplos, gcat1[a&0x03], gcat2[b&0x03]);
		a>>=2; b>>=2;
		count_2bit(diplos, gcat1[a&0x03], gcat2[b&0x03]);
	}

	len1 = len1&0x03;

	/* Handle len%4 remaining genotypes */
	if(len1)
	{
		unsigned char a = g1[i];
		unsigned char b = g2[i];
		           count_2bit(diplos, gcat1[(a>>6)&0x03], gcat2[(b>>6)&0x03]);
		if(len1>1) count_2bit(diplos, gcat1[(a>>4)&0x03], gcat2[(b>>4)&0x03]);
		if(len1>2) count_2bit(diplos, gcat1[(a>>2)&0x03], gcat2[(b>>2)&0x03]);
	}

done:
	if(gcat1) free(gcat1);
	if(gcat2) free(gcat2);
	return diplos_to_haplos(diplos);

error:
	if(gcat1) free(gcat1);
	if(gcat2) free(gcat2);
	return NULL;
}

PyObject *
count_haplotypes(PyObject *self, PyObject *args)
{
	PyObject *genos1, *genos2;
	GenotypeObject *g1=NULL, *g2=NULL;
	UnphasedMarkerModelObject *model1=NULL, *model2=NULL;
	Py_ssize_t len1, len2, i, a, b, diplos[9];
	Py_ssize_t a1,a2,b1,b2,x,mult;

	if(!PyArg_ParseTuple(args, "OO", &genos1, &genos2))
		return NULL;

	if(GenotypeArray_CheckExact(genos1) && GenotypeArray_CheckExact(genos2))
	{
		GenotypeArrayObject *g1 = (GenotypeArrayObject *)genos1;
		GenotypeArrayObject *g2 = (GenotypeArrayObject *)genos2;
		unsigned int *offsets = (unsigned int *)PyArray_DATA(g1->descriptor->offsets);
		if(offsets[0] == 0 && g1->descriptor->homogeneous == 2 && g2->descriptor->homogeneous == 2)
			return count_haplotypes_2bit(self, args);
	}

	len1 = PySequence_Size(genos1);
	if(len1==-1) goto error;
	len2 = PySequence_Size(genos2);
	if(len2==-1) goto error;

	if(len1 != len2)
	{
		PyErr_Format(GenotypeRepresentationError,"genotype array sizes do not match: %zd != %zd",
		             len1, len2);
		goto error;
	}

	a1 = a2 = -1;
	b1 = b2 = -1;

	for(i = 0; i<9; ++i)
		diplos[i] = 0;

	if(!len1) goto done;

	for(i = 0; i < len1; ++i)
	{
		g1 = (GenotypeObject *)PySequence_GetItem(genos1,i);
		if(!g1) goto error;
		g2 = (GenotypeObject *)PySequence_GetItem(genos2,i);
		if(!g2) goto error;

		if(!Genotype_CheckExact(g1))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos1 at index %zd", i);
			goto error;
		}

		if(!Genotype_CheckExact(g2))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos2 at index %zd", i);
			goto error;
		}

		if(!model1)
		{
			model1 = g1->model;
			if(!model1 || !UnphasedMarkerModel_Check(model1))
			{
				PyErr_SetString(PyExc_TypeError,"invalid genotype model");
				goto error;
			}
		}

		if(!model2)
		{
			model2 = g2->model;
			if(!model2 || !UnphasedMarkerModel_Check(model2))
			{
				PyErr_SetString(PyExc_TypeError,"invalid genotype model");
				goto error;
			}
		}

		if(model1 != g1->model || model2 != g2->model)
		{
			PyErr_Format(GenotypeRepresentationError,
			    "genotype models do not match at index %zd", i);
			goto error;
		}

		/* Skip missing and hemizygotes */
		if(!g1->index || !g2->index)
			continue;

		if(!g1->allele1_index || !g1->allele2_index || !g2->allele1_index || !g2->allele2_index)
		{
			PyErr_SetString(PyExc_ValueError, "Hemizygote LD estimation is not currently supported");
			goto error;
		}

		x = find_het(&a1, &a2, g1) + find_het(&b1, &b2, g2);

		if(x != 2)
		{
			PyErr_SetString(PyExc_ValueError, "invalid genotypes: loci may have no more than 2 alleles");
			goto error;
		}

		if     (homozygous_for(g1,a1)) a = 0;
		else if(homozygous_for(g1,a2)) a = 2;
		else                           a = 1;

		if     (homozygous_for(g2,b1)) b = 0;
		else if(homozygous_for(g2,b2)) b = 2;
		else                           b = 1;

		mult = 1;
		if(homozygous(g1) && homozygous(g2))
			mult = 2;

		diplos[3*a + b] += mult;

		Py_DECREF(g1);
		Py_DECREF(g2);
		g1 = g2 = NULL;
	}

done:
	return diplos_to_haplos(diplos);

error:
	Py_XDECREF(g1);
	Py_XDECREF(g2);
	return NULL;
}

PyObject *
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
		d_max =  dmin( p*(1-q), (1-p)*q );
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
