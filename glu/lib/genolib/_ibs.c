/*
File:          _ibs.c

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:

Abstract:      IBS statistics

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.
See GLU license for terms by running: glu license
*/

#include <Python.h>
#include <math.h>
#include "_genoarray.h"


PyObject *
genoarray_concordance_8bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	const unsigned char *g1, *g2;
	const unsigned int *offsets;
	Py_ssize_t len, concordant, comparisons, i;

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 8)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 8 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	g1 = genos1->data;
	g2 = genos2->data;

	concordant = comparisons = 0;
	for(i = 0; i < len; ++i)
	{
		const unsigned char a = g1[i];
		const unsigned char b = g2[i];

		/* If both genotypes are not missing */
		if(a && b)
		{
			if(a==b) concordant += 1;
			comparisons += 1;
		}
	}
	return Py_BuildValue("(ii)", concordant, comparisons);
}

PyObject *
genoarray_concordance_4bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	const unsigned char *g1, *g2;
	const unsigned int *offsets;
	Py_ssize_t len, concordant, comparisons, i;

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 4)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 4 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	g1 = genos1->data;
	g2 = genos2->data;

	concordant = comparisons = 0;
	for(i = 0; i < len/2; ++i)
	{
		const unsigned char a  = g1[i];
		const unsigned char b  = g2[i];
		const unsigned char a1 = a&0xF0;
		const unsigned char b1 = b&0xF0;
		const unsigned char a2 = a&0x0F;
		const unsigned char b2 = b&0x0F;

		/* If both genotypes are not missing */
		if( a1 && b1 )
		{
			if(a1==b1) concordant += 1;
			comparisons += 1;
		}

		if( a2 && b2 )
		{
			if(a2==b2) concordant += 1;
			comparisons += 1;
		}
	}

	if(len&1)
	{
		const unsigned char a  = g1[i];
		const unsigned char b  = g2[i];
		const unsigned char a1 = a&0xF0;
		const unsigned char b1 = b&0xF0;

		/* If both genotypes are not missing */
		if( a1 && b1 )
		{
			if(a1==b1) concordant += 1;
			comparisons += 1;
		}

	}
	return Py_BuildValue("(ii)", concordant, comparisons);
}

PyObject *
genoarray_concordance_2bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	const unsigned char *g1, *g2;
	const unsigned int *offsets;
	Py_ssize_t len, concordant, comparisons, i;

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 2)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 2 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	g1 = genos1->data;
	g2 = genos2->data;

	concordant = comparisons = 0;
	for(i = 0; i < len/4; ++i)
	{
		const unsigned char a  = g1[i];
		const unsigned char b  = g2[i];
		const unsigned char a1 = a&0xC0;
		const unsigned char a2 = a&0x30;
		const unsigned char a3 = a&0x0C;
		const unsigned char a4 = a&0x03;
		const unsigned char b1 = b&0xC0;
		const unsigned char b2 = b&0x30;
		const unsigned char b3 = b&0x0C;
		const unsigned char b4 = b&0x03;

		/* If both genotypes are not missing */
		if( a1 && b1 )
		{
			if(a1==b1) concordant += 1;
			comparisons += 1;
		}

		if( a2 && b2 )
		{
			if(a2==b2) concordant += 1;
			comparisons += 1;
		}
		if( a3 && b3 )
		{
			if(a3==b3) concordant += 1;
			comparisons += 1;
		}

		if( a4 && b4 )
		{
			if(a4==b4) concordant += 1;
			comparisons += 1;
		}
	}

	len = len&3;

	/* Handle len%4 remaining genotypes */
	if(len)
	{
		const unsigned char a  = g1[i];
		const unsigned char b  = g2[i];
		const unsigned char a1 = a&0xC0;
		const unsigned char a2 = a&0x30;
		const unsigned char a3 = a&0x0C;
		const unsigned char b1 = b&0xC0;
		const unsigned char b2 = b&0x30;
		const unsigned char b3 = b&0x0C;

		if( a1 && b1 )
		{
			if(a1==b1) concordant += 1;
			comparisons += 1;
		}

		if( len>1 && a2 && b2 )
		{
			if(a2==b2) concordant += 1;
			comparisons += 1;
		}
		if( len>2 && a3 && b3 )
		{
			if(a3==b3) concordant += 1;
			comparisons += 1;
		}
	}

	return Py_BuildValue("(ii)", concordant, comparisons);
}

PyObject *
genoarray_concordance(PyObject *self, PyObject *args)
{
	PyObject *genos1, *genos2;
	PyObject **items1, **items2;
	Py_ssize_t len1, len2, concordant, comparisons, i;

	if(!PyArg_ParseTuple(args, "OO", &genos1, &genos2))
		return NULL;

	if(GenotypeArray_CheckExact(genos1) && GenotypeArray_CheckExact(genos2))
	{
		GenotypeArrayObject *g1 = (GenotypeArrayObject *)genos1;
		GenotypeArrayObject *g2 = (GenotypeArrayObject *)genos2;
		unsigned int   *offsets = (unsigned int *)PyArray_DATA(g1->descriptor->offsets);
		if(g1->descriptor == g2->descriptor &&
		   descr_length(g1->descriptor) > 0 && offsets[0] == 0)
		{
			if(g1->descriptor->homogeneous == 8)
				return genoarray_concordance_8bit(self, args);
			if(g1->descriptor->homogeneous == 4)
				return genoarray_concordance_4bit(self, args);
			if(g1->descriptor->homogeneous == 2)
				return genoarray_concordance_2bit(self, args);
		}
	}

	genos1 = PySequence_Fast(genos1,"cannot convert genos1 into a sequence");
	if(!genos1) return NULL;

	genos2 = PySequence_Fast(genos2,"cannot convert genos2 into a sequence");
	if(!genos2) goto error;

	len1 = PySequence_Fast_GET_SIZE(genos1);
	len2 = PySequence_Fast_GET_SIZE(genos2);

	if(len1 != len2)
	{
		PyErr_Format(GenotypeRepresentationError,"genotype array sizes do not match: %zd != %zd",
		             len1, len2);
		goto error;
	}

	items1 = PySequence_Fast_ITEMS(genos1);
	items2 = PySequence_Fast_ITEMS(genos2);

	concordant = comparisons = 0;
	for(i = 0; i < len1; ++i)
	{
		GenotypeObject *a = (GenotypeObject *)items1[i];
		GenotypeObject *b = (GenotypeObject *)items2[i];

		if(!a || !Genotype_CheckExact(a))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos1 at index %zd", i);
			goto error;
		}

		if(!b || !Genotype_CheckExact(b))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos2 at index %zd", i);
			goto error;
		}

		/* If both genotypes are not missing */
		if(a->index && b->index)
		{
			if(a->index==b->index) concordant += 1;
			comparisons += 1;
		}
	}

	Py_DECREF(genos1);
	Py_DECREF(genos2);

	return Py_BuildValue("(ii)", concordant, comparisons);

error:
	Py_XDECREF(genos1);
	Py_XDECREF(genos2);
	return NULL;
}

/******************************************************************************************************/

int
genotype_ibs_from_indices(UnphasedMarkerModelObject *model, Py_ssize_t index1, Py_ssize_t index2)
{
	GenotypeObject *geno1, *geno2;

	assert(model && UnphasedMarkerModel_CheckExact(model));
	assert(model->genotypes && PyList_Check(model->genotypes));
	assert(index1 < PyList_Size(model->genotypes));
	assert(index2 < PyList_Size(model->genotypes));

	if(index1 == index2)
		return 2;

	geno1 = (GenotypeObject *)PyList_GET_ITEM(model->genotypes, index1);
	geno2 = (GenotypeObject *)PyList_GET_ITEM(model->genotypes, index2);

	assert(geno1 && Genotype_CheckExact(geno1));
	assert(geno2 && Genotype_CheckExact(geno2));
	assert(model == geno1->model);
	assert(model == geno2->model);

	if(geno1->allele1_index == geno2->allele1_index ||
	   geno1->allele1_index == geno2->allele2_index ||
	   geno1->allele2_index == geno2->allele1_index ||
	   geno1->allele2_index == geno2->allele2_index)
		return 1;

	return 0;
}

inline
UnphasedMarkerModelObject *
get_model(GenotypeArrayDescriptorObject *descr, Py_ssize_t i)
{
	UnphasedMarkerModelObject *model;

	model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(descr->models, i); /* borrowed ref */
	if(!model || !UnphasedMarkerModel_CheckExact(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return NULL;
	}
	return model;
}

PyObject *
genoarray_ibs_8bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	UnphasedMarkerModelObject *model;
	const unsigned char *geno1_data, *geno2_data;
	const unsigned int *offsets;
	Py_ssize_t len, i, ibs[3];

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 8)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 8 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	geno1_data = genos1->data;
	geno2_data = genos2->data;

	ibs[0] = ibs[1] = ibs[2] = 0;
	for(i = 0; i < len; ++i)
	{
		int a = geno1_data[i];
		int b = geno2_data[i];

		model = get_model(genos1->descriptor, i);
		if(!model) return NULL;

		/* If both genotypes are not missing */
		if(a && b) ibs[genotype_ibs_from_indices(model, a, b)] += 1;
	}
	return Py_BuildValue("(iii)", ibs[0], ibs[1], ibs[2]);
}

PyObject *
genoarray_ibs_4bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	UnphasedMarkerModelObject *model1, *model2;
	const unsigned char *geno1_data, *geno2_data;
	const unsigned int *offsets;
	Py_ssize_t len, i, ibs[3];

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 4)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 4 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	geno1_data = genos1->data;
	geno2_data = genos2->data;

	ibs[0] = ibs[1] = ibs[2] = 0;
	for(i=0; i < len/2; ++i)
	{
		const unsigned char a  = geno1_data[i];
		const unsigned char b  = geno2_data[i];
		const unsigned char a1 = (a>>4)&0x0F;
		const unsigned char b1 = (b>>4)&0x0F;
		const unsigned char a2 = (a   )&0x0F;
		const unsigned char b2 = (b   )&0x0F;

		model1 = get_model(genos1->descriptor, 2*i);
		model2 = get_model(genos1->descriptor, 2*i+1);
		if(!model1 || !model2) return NULL;

		/* If both genotypes are not missing */
		if(a1 && b1) ibs[genotype_ibs_from_indices(model1, a1, b1)] += 1;
		if(a2 && b2) ibs[genotype_ibs_from_indices(model2, a2, b2)] += 1;
	}

	if(len&1)
	{
		const unsigned char a  = geno1_data[i];
		const unsigned char b  = geno2_data[i];
		const unsigned char a1 = (a>>4)&0x0F;
		const unsigned char b1 = (b>>4)&0x0F;

		model1 = get_model(genos1->descriptor, 2*i);
		if(!model1) return NULL;

		/* If both genotypes are not missing */
		if(a1 && b1) ibs[genotype_ibs_from_indices(model1, a1, b1)] += 1;
	}
	return Py_BuildValue("(iii)", ibs[0], ibs[1], ibs[2]);
}

PyObject *
genoarray_ibs_2bit(PyObject *self, PyObject *args)
{
	GenotypeArrayObject *genos1, *genos2;
	UnphasedMarkerModelObject *model1, *model2, *model3, *model4;
	const unsigned char *geno1_data, *geno2_data;
	const unsigned int *offsets;
	Py_ssize_t len, i, ibs[3];

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

	if(genos1->descriptor != genos2->descriptor)
	{
		PyErr_SetString(GenotypeRepresentationError,"genos1 and genos2 must share the same descriptor");
		return NULL;
	}

	if(genos1->descriptor->homogeneous != 2)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype arrays must have homogeneous 2 bit width");
		return NULL;
	}

	len = descr_length(genos1->descriptor);
	if(len < 0) return NULL;

	offsets = (const unsigned int *)PyArray_DATA(genos1->descriptor->offsets);

	if(offsets[0] != 0)
	{
		PyErr_SetString(GenotypeRepresentationError, "genotype array data must begin with zero offset");
		return NULL;
	}

	geno1_data = genos1->data;
	geno2_data = genos2->data;

	ibs[0] = ibs[1] = ibs[2] = 0;
	for(i = 0; i < len/4; ++i)
	{
		const unsigned char a  = geno1_data[i];
		const unsigned char b  = geno2_data[i];
		const unsigned char a1 = (a>>6)&0x03;
		const unsigned char a2 = (a>>4)&0x03;
		const unsigned char a3 = (a>>2)&0x03;
		const unsigned char a4 = (a   )&0x03;
		const unsigned char b1 = (b>>6)&0x03;
		const unsigned char b2 = (b>>4)&0x03;
		const unsigned char b3 = (b>>2)&0x03;
		const unsigned char b4 = (b   )&0x03;

		model1 = get_model(genos1->descriptor, 4*i);
		if(!model1) return NULL;
		model2 = get_model(genos1->descriptor, 4*i+1);
		if(!model2) return NULL;
		model3 = get_model(genos1->descriptor, 4*i+2);
		if(!model3) return NULL;
		model4 = get_model(genos1->descriptor, 4*i+3);
		if(!model4) return NULL;

		/* If both genotypes are not missing */
		if(a1 && b1) ibs[genotype_ibs_from_indices(model1, a1, b1)] += 1;
		if(a2 && b2) ibs[genotype_ibs_from_indices(model2, a2, b2)] += 1;
		if(a3 && b3) ibs[genotype_ibs_from_indices(model3, a3, b3)] += 1;
		if(a4 && b4) ibs[genotype_ibs_from_indices(model4, a4, b4)] += 1;
	}

	len = len&3;

	/* Handle len%4 remaining genotypes */
	if(len)
	{
		const unsigned char a  = geno1_data[i];
		const unsigned char b  = geno2_data[i];
		const unsigned char a1 = (a>>6)&0x03;
		const unsigned char a2 = (a>>4)&0x03;
		const unsigned char a3 = (a>>2)&0x03;
		const unsigned char b1 = (b>>6)&0x03;
		const unsigned char b2 = (b>>4)&0x03;
		const unsigned char b3 = (b>>2)&0x03;

		if( a1 && b1 )
		{
			model1 = get_model(genos1->descriptor, 4*i);
			if(!model1) return NULL;
			ibs[genotype_ibs_from_indices(model1, a1, b1)] += 1;
		}

		if( len>1 && a2 && b2 )
		{
			model2 = get_model(genos1->descriptor, 4*i+1);
			if(!model2) return NULL;
			ibs[genotype_ibs_from_indices(model2, a2, b2)] += 1;
		}
		if( len>2 && a3 && b3 )
		{
			model3 = get_model(genos1->descriptor, 4*i+2);
			if(!model3) return NULL;
			ibs[genotype_ibs_from_indices(model3, a3, b3)] += 1;
		}
	}

	return Py_BuildValue("(iii)", ibs[0], ibs[1], ibs[2]);
}

PyObject *
genoarray_ibs(PyObject *self, PyObject *args)
{
	PyObject *genos1, *genos2;
	PyObject **items1, **items2;
	Py_ssize_t len1, len2, i, ibs[3];

	if(!PyArg_ParseTuple(args, "OO", &genos1, &genos2))
		return NULL;

	if(GenotypeArray_CheckExact(genos1) && GenotypeArray_CheckExact(genos2))
	{
		GenotypeArrayObject *g1 = (GenotypeArrayObject *)genos1;
		GenotypeArrayObject *g2 = (GenotypeArrayObject *)genos2;
		unsigned int   *offsets = (unsigned int *)PyArray_DATA(g1->descriptor->offsets);
		if(g1->descriptor == g2->descriptor &&
		   descr_length(g1->descriptor) > 0 && offsets[0] == 0)
		{
			if(g1->descriptor->homogeneous == 8)
				return genoarray_ibs_8bit(self, args);
			if(g1->descriptor->homogeneous == 4)
				return genoarray_ibs_4bit(self, args);
			if(g1->descriptor->homogeneous == 2)
				return genoarray_ibs_2bit(self, args);
		}
	}

	genos1 = PySequence_Fast(genos1,"cannot convert genos1 into a sequence");
	if(!genos1) return NULL;

	genos2 = PySequence_Fast(genos2,"cannot convert genos2 into a sequence");
	if(!genos2) goto error;

	len1 = PySequence_Fast_GET_SIZE(genos1);
	len2 = PySequence_Fast_GET_SIZE(genos2);

	if(len1 != len2)
	{
		PyErr_Format(GenotypeRepresentationError,"genotype array sizes do not match: %zd != %zd",
		             len1, len2);
		goto error;
	}

	items1 = PySequence_Fast_ITEMS(genos1);
	items2 = PySequence_Fast_ITEMS(genos2);

	ibs[0] = ibs[1] = ibs[2] = 0;
	for(i = 0; i < len1; ++i)
	{
		GenotypeObject *a = (GenotypeObject *)items1[i];
		GenotypeObject *b = (GenotypeObject *)items2[i];
		int n;

		if(!a || !Genotype_CheckExact(a))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos1 at index %zd", i);
			goto error;
		}

		if(!b || !Genotype_CheckExact(b))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos2 at index %zd", i);
			goto error;
		}

		/* If both genotypes are not missing */
		if(a->index && b->index)
		{
			if(a==b)
				n = 2;
			else if( a->allele1_index == b->allele1_index ||
			         a->allele1_index == b->allele2_index ||
			         a->allele2_index == b->allele1_index ||
			         a->allele2_index == b->allele2_index )
				n = 1;
			else
				n = 0;

			ibs[n] += 1;
		}
	}

	Py_DECREF(genos1);
	Py_DECREF(genos2);

	return Py_BuildValue("(iii)", ibs[0], ibs[1], ibs[2]);

error:
	Py_XDECREF(genos1);
	Py_XDECREF(genos2);
	return NULL;
}
