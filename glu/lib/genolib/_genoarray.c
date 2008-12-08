/*
File:          _genoarray.c

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:

Abstract:      fast implementation of a bit-packed genotype array

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2008, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.
See GLU license for terms by running: glu license
*/

#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"
#include "bitarrayc.h"
#include "_genoarray.h"

#ifdef STDC_HEADERS
#include <stddef.h>
#else
#include <sys/types.h>		/* For size_t */
#endif

#include <math.h>

#define MERGE_UNAMBIGUOUS 0
#define MERGE_CONCORDANT  1
#define MERGE_CONSENSUS   2
#define MERGE_DISCORDANT  3
#define MERGE_MISSING     4

#define CountArray_Check(op,len)               (PyArray_Check((op)) && PyArray_NDIM((op))== 1 && PyArray_DIMS((op))[0]==(len) \
                                                && PyArray_TYPE((op))==PyArray_LONG)
#define CountArray_Check2(op,len1,len2)        (PyArray_Check((op)) && PyArray_NDIM((op))== 2 && PyArray_DIMS((op))[0]==(len1) \
                                                && PyArray_DIMS((op))[1]==(len2) && PyArray_TYPE((op))==PyArray_LONG)
#define FastSequenceCheck(op)                  (PyTuple_Check(op) || PyList_Check(op))

/* Exceptions */
PyObject *GenotypeLookupError;
PyObject *GenotypeRepresentationError;

/******************************************************************************************************/

static int
genotype_init(GenotypeObject *self, PyObject *args, PyObject *kwds)
{
	UnphasedMarkerModelObject *model;
	Py_ssize_t index, allele1_index, allele2_index, allele_len, genotype_len;

	static char *kwlist[] = {"model", "allele1_index", "allele2_index", "index", NULL};

	self->allele1_index = self->allele2_index = self->index = 0;
	Py_CLEAR(self->model);

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "Onnn", kwlist,
		&model, &allele1_index, &allele2_index, &index))
		return -1;

	if(!model || !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return -1;
	}

	allele_len = PyList_Size(model->alleles);
	if(allele_len == -1) return -1;

	genotype_len = PyList_Size(model->genotypes);
	if(genotype_len == -1) return -1;

	if(allele1_index < 0 || allele1_index >= allele_len ||
	   allele2_index < 0 || allele2_index >= allele_len)
	{
		PyErr_SetString(PyExc_IndexError,"invalid allele index");
		return -1;
	}

	/* We allow the genotype to be the next one in the list */
	if(index < 0 || index > genotype_len)
	{
		PyErr_SetString(PyExc_IndexError,"invalid genotype index");
		return -1;
	}

	Py_INCREF(model);
	self->model   = model;
	self->index   = index;
	self->allele1_index = allele1_index;
	self->allele2_index = allele2_index;

	return 0;
}

static int
genotype_clear(GenotypeObject *self)
{
	Py_CLEAR(self->model);
	return 0;
}

static void
genotype_dealloc(GenotypeObject *self)
{
	genotype_clear(self);
	self->ob_type->tp_free((PyObject *)self);
}

static int
genotype_traverse(GenotypeObject *self, visitproc visit, void *arg)
{
	Py_VISIT(self->model);
	return 0;
}

static PyObject *
genotype_allele1_get(GenotypeObject *self)
{
	PyObject *allele1 = PyList_GetItem(self->model->alleles, self->allele1_index); /* borrowed ref */
	Py_XINCREF(allele1);
	return allele1;

}

static PyObject *
genotype_allele2_get(GenotypeObject *self)
{
	PyObject *allele2 = PyList_GetItem(self->model->alleles, self->allele2_index); /* borrowed ref */
	Py_XINCREF(allele2);
	return allele2;
}

static PyObject *
genotype_alleles(GenotypeObject *self)
{
	PyObject *result = PyList_GetItem(self->model->genotuples, self->index); /* borrowed ref */
	Py_XINCREF(result);
	return result;
}

static PyObject *
genotype_repr(GenotypeObject *self)
{
	PyObject *alleles, *result;

	alleles = genotype_alleles(self);
	if(!alleles) return NULL;

	result = PyObject_Repr(alleles);
	Py_DECREF(alleles);
	return result;
}

static PyObject *
genotype_heterozygote(GenotypeObject *self)
{
	if( self->allele1_index && self->allele2_index
	 && self->allele1_index != self->allele2_index )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_homozygote(GenotypeObject *self)
{
	if( self->allele1_index && self->allele1_index == self->allele2_index )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_hemizygote(GenotypeObject *self)
{
	if( !!self->allele1_index ^ !!self->allele2_index )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_missing(GenotypeObject *self)
{
	if( !self->allele1_index && !self->allele2_index )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static Py_ssize_t
genotype_category(GenotypeObject *self)
{
	Py_ssize_t n = !!self->allele1_index + !!self->allele2_index;
	if(n==2 && self->allele1_index != self->allele2_index)
		n++;
	return n;
}

static PyObject *
genotype_category_get(GenotypeObject *self)
{
	return PyInt_FromSsize_t(genotype_category(self));
}

static Py_ssize_t
genotype_length(GenotypeObject *self)
{
	return !!self->allele1_index + !!self->allele2_index;
}

static PyObject *
genotype_item(GenotypeObject *self, Py_ssize_t item)
{
	PyObject *allele;
	if(item == 0)
		allele = PyList_GetItem(self->model->alleles, self->allele1_index); /* borrowed ref */
	else if(item == 1)
		allele = PyList_GetItem(self->model->alleles, self->allele2_index); /* borrowed ref */
	else
	{
		PyErr_SetString(PyExc_IndexError,"genotype index out of range");
		return NULL;
	}
	Py_XINCREF(allele);
	return allele;
}

static int
genotype_contains(GenotypeObject *self, PyObject *allele)
{
	int cmp;
	PyObject *allele1, *allele2;

	allele1 = PyList_GetItem(self->model->alleles, self->allele1_index); /* borrowed ref */
	cmp = PyObject_RichCompareBool(allele, allele1, Py_EQ);
	if(cmp == 0)
	{
		allele2 = PyList_GetItem(self->model->alleles, self->allele2_index); /* borrowed ref */
		cmp = PyObject_RichCompareBool(allele, allele2, Py_EQ);
	}
	return cmp;
}

static long
genotype_hash(GenotypeObject *self)
{
	PyObject *alleles = PyList_GetItem(self->model->genotuples, self->index); /* borrowed ref */
	if(!alleles) return -1;
	return PyObject_Hash(alleles);
}

static PyObject *
genotype_normalized_tuple(PyObject *geno)
{
	PyObject *allele1, *allele2;
	int res;

	if(Genotype_CheckExact(geno))
		return genotype_alleles( (GenotypeObject *)geno );

	if(!PyTuple_Check(geno)  || PyTuple_GET_SIZE(geno)  != 2)
	{
		PyErr_SetString(GenotypeRepresentationError,"genotype must be specified as a 2-tuple");
		return NULL;
	}

	allele1 = PyTuple_GET_ITEM(geno, 0); /* borrowed ref */
	allele2 = PyTuple_GET_ITEM(geno, 1); /* borrowed ref */
	res = PyObject_RichCompareBool(allele1, allele2, Py_GT);

	if(res == -1) return NULL;
	else if(res == 1)
		return PyTuple_Pack(2, allele2, allele1);
	else
	{
		Py_INCREF(geno);
		return geno;
	}
}

static PyObject *
genotype_richcompare(PyObject *self, PyObject *other, int op)
{
	PyObject *geno1=NULL, *geno2=NULL, *ret=NULL;

	geno1 = genotype_normalized_tuple(self);
	if(!geno1) goto done;
	geno2 = genotype_normalized_tuple(other);
	if(!geno2) goto done;

	ret=PyObject_RichCompare(geno1, geno2, op);

done:
	Py_XDECREF(geno1);
	Py_XDECREF(geno2);
	return ret;
}

static PyMemberDef genotype_members[] = {
	{"model",   T_OBJECT_EX, offsetof(GenotypeObject, model),   RO, "UnphasedMarkerModel"},
	{"index",   T_UINT,      offsetof(GenotypeObject, index),   RO, "index"},
	{"allele1_index",	T_UINT, offsetof(GenotypeObject, allele1_index), RO, "allele1_index"},
	{"allele2_index",	T_UINT, offsetof(GenotypeObject, allele2_index), RO, "allele2_index"},
	{NULL}  /* Sentinel */
};

static PyGetSetDef genotype_getsetlist[] = {
	{"category",	(getter)genotype_category_get,	NULL,	NULL},
	{"allele1",	(getter)genotype_allele1_get,	NULL,	NULL},
	{"allele2",	(getter)genotype_allele2_get,	NULL,	NULL},
	{NULL, NULL, NULL, NULL},  /* Sentinel */
};

static PyMethodDef genotype_methods[] = {
	{"alleles",      (PyCFunction)genotype_alleles,      METH_NOARGS,
		"Returns a tuple of alleles"},
	{"heterozygote", (PyCFunction)genotype_heterozygote, METH_NOARGS,
		"Predicate to test if the genotype is heterozygous"},
	{"homozygote",   (PyCFunction)genotype_homozygote,   METH_NOARGS,
		"Predicate to test if the genotype is homozygous"},
	{"hemizygote",   (PyCFunction)genotype_hemizygote,   METH_NOARGS,
		"Predicate to test if the genotype is hemizygous"},
	{"missing",      (PyCFunction)genotype_missing,      METH_NOARGS,
		"Predicate to test if the genotype is missing (uninformative)"},
{NULL}  /* Sentinel */
};

static PySequenceMethods genotype_as_sequence = {
	(lenfunc)genotype_length,		/* sq_length */
	0,					/* sq_concat */
	0,					/* sq_repeat */
	(ssizeargfunc)genotype_item,		/* sq_item */
	0,					/* sq_slice */
	0,					/* sq_ass_item */
	0,					/* sq_ass_slice */
	(objobjproc)genotype_contains,		/* sq_contains */
	0,					/* sq_inplace_concat */
	0,					/* sq_inplace_repeat */
};

PyTypeObject GenotypeType = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size           */
	"Genotype",				/* tp_name           */
	sizeof(GenotypeObject),			/* tp_basicsize      */
	0,					/* tp_itemsize       */
	(destructor)genotype_dealloc,		/* tp_dealloc        */
	0,					/* tp_print          */
	0,					/* tp_getattr        */
	0,					/* tp_setattr        */
	0,					/* tp_compare        */
	(reprfunc)genotype_repr,		/* tp_repr           */
	0,					/* tp_as_number      */
	&genotype_as_sequence,			/* tp_as_sequence    */
	0,					/* tp_as_mapping     */
	(hashfunc)genotype_hash,		/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT,			/* tp_flags          */
	"Genotype objects",			/* tp_doc            */
	(traverseproc)genotype_traverse,	/* tp_traverse       */
	(inquiry)genotype_clear,		/* tp_clear          */
	genotype_richcompare,			/* tp_richcompare    */
	0,					/* tp_weaklistoffset */
	0,					/* tp_iter           */
	0,					/* tp_iternext       */
	genotype_methods,			/* tp_methods        */
	genotype_members,			/* tp_members        */
	genotype_getsetlist,			/* tp_getset         */
	0,					/* tp_base           */
	0,					/* tp_dict           */
	0,					/* tp_descr_get      */
	0,					/* tp_descr_set      */
	0,					/* tp_dictoffset     */
	(initproc)genotype_init,		/* tp_init           */
	0,					/* tp_alloc          */
	0,					/* tp_new            */
	0,					/* tp_free           */
};

/******************************************************************************************************/

static int
descr_clear(GenotypeArrayDescriptorObject *self)
{
	Py_CLEAR(self->models);
	Py_CLEAR(self->offsets);
	return 0;
}

static int
descr_init(GenotypeArrayDescriptorObject *self, PyObject *args, PyObject *kwds)
{
	PyObject *models;
	PyArrayObject *offsets;
	UnphasedMarkerModelObject *model;
	Py_ssize_t n,i;
	unsigned int initial_offset = 0;
	unsigned int *offset_data;
	npy_intp dims;
	static char *kwlist[] = {"models", "initial_offset", NULL};

	descr_clear(self);
	self->byte_size = self->bit_size = self->max_bit_size = self->homogeneous = 0;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|n:GenotypeArrayDescriptor", kwlist,
	                                             &models, &initial_offset))
		return -1;

	if(!models || !PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		return -1;
	}

	n = PyList_GET_SIZE(models);
	if(n == -1) return -1;

	dims = n+1;
	offsets = (PyArrayObject *)PyArray_SimpleNew(1, &dims, NPY_UINT);
	if(!offsets) return -1;

	offset_data = (unsigned int *)PyArray_DATA(offsets);
	offset_data[0] = initial_offset;

	for(i=0; i<n; ++i)
	{
		model = (UnphasedMarkerModelObject *)PyList_GetItem(models, i); /* borrowed ref */
		if(!model) goto error;
		if(!UnphasedMarkerModel_Check(model) || model->bit_size > 32)
		{
			PyErr_SetString(PyExc_TypeError,"invalid genotype model");
			goto error;
		}
		offset_data[i+1] = offset_data[i] + model->bit_size;

		if(model->bit_size > self->max_bit_size)
			self->max_bit_size = model->bit_size;

		if(i==0)
			self->homogeneous = model->bit_size;
		else if(self->homogeneous != model->bit_size)
			self->homogeneous = 0;
	}

	Py_INCREF(models);
	self->models    = models;
	self->offsets   = offsets;
	self->bit_size  = offset_data[n];
	self->byte_size = (self->bit_size+7)/8;

	return 0;

error:
	Py_DECREF(offsets);
	return -1;
}

static void
descr_dealloc(GenotypeArrayDescriptorObject *self)
{
	descr_clear(self);
	self->ob_type->tp_free((PyObject *)self);
}

static int
descr_traverse(GenotypeArrayDescriptorObject *self, visitproc visit, void *arg)
{
	Py_VISIT(self->models);
	Py_VISIT(self->offsets);
	return 0;
}

static Py_ssize_t
descr_length(GenotypeArrayDescriptorObject *self)
{
	Py_ssize_t n;

	if(!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return -1;
	}

	n = PyList_Size(self->models);

	if(n == -1 && PyErr_Occurred())
		return -1;

	return n;
}

static PyObject *
descr_iter(GenotypeArrayDescriptorObject *self)
{
	if(!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return NULL;
	}

	return PyObject_GetIter(self->models);
}

static PyObject *
descr_item(GenotypeArrayDescriptorObject *self, Py_ssize_t i)
{
	if(!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return NULL;
	}
	return PySequence_GetItem( (PyObject *)self->models, i); /* new ref */
}

static PyObject *
descr_subscript(GenotypeArrayDescriptorObject *self, PyObject *item)
{
	if(!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return NULL;
	}
	return PyObject_GetItem( (PyObject *)self->models, item); /* new ref */
}

static int
descr_ass_item(GenotypeArrayDescriptorObject *self, Py_ssize_t item, UnphasedMarkerModelObject *new_model)
{
	UnphasedMarkerModelObject *old_model;
	Py_ssize_t old_len, new_len;
	int ret;

	if(!new_model || !UnphasedMarkerModel_Check(new_model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return -1;
	}

	old_model = (UnphasedMarkerModelObject *)PyList_GetItem(self->models, item); /* borrowed ref */

	if(!old_model) return -1;


	if(old_model->bit_size != new_model->bit_size)
	{
		PyErr_SetString(GenotypeRepresentationError,"Model may not be replaced due to incompatible size");
		return -1;
	}

	old_len = PyList_Size(old_model->alleles);
	if(old_len == -1) return -1;
	new_len = PyList_Size(new_model->alleles);
	if(new_len == -1) return -1;

	if(old_len > new_len)
	{
		PyErr_SetString(GenotypeRepresentationError,"Model may not be replaced since new model defines too few alleles");
		return -1;
	}

	old_len = PyList_Size(old_model->genotypes);
	if(old_len == -1) return -1;
	new_len = PyList_Size(new_model->genotypes);
	if(new_len == -1) return -1;

	if(old_len > new_len)
	{
		PyErr_SetString(GenotypeRepresentationError,"Model may not be replaced since new model defines too few genotypes");
		return -1;
	}

	ret = PyList_SetItem(self->models, item, (PyObject *)new_model);

	if(ret==0) Py_INCREF(new_model);

	return ret;
}

static int
descr_ass_slice(GenotypeArrayDescriptorObject *self, PySliceObject *slice, PyObject *value)
{
	Py_ssize_t i, j, n;
	Py_ssize_t start, stop, step, slicelength;
	UnphasedMarkerModelObject *new_model;
	PyObject *seq;
	int ret = -1;

	if(!value)
	{
		/* delete slice */
		PyErr_SetString(PyExc_TypeError, "GenotypeArrayDescriptor objects do not support model deletion");
		return -1;
	}

	n = PyList_Size(self->models);
	if(n==-1) return -1;

	if(PySlice_GetIndicesEx(slice, n, &start, &stop, &step, &slicelength) < 0)
		return -1;

	seq = PyObject_GetIter(value);
	if(seq == NULL) goto error;

	for (i=start, j=0; j<slicelength; i+=step, ++j)
	{
		new_model = (UnphasedMarkerModelObject *)PyIter_Next(seq);

		if(!new_model) break;

		if( descr_ass_item(self, i, new_model) == -1 )
		{
			Py_DECREF(new_model);
			goto error;
		}
		Py_DECREF(new_model);
	}

	if(PyErr_Occurred())
	{
		if(j<slicelength && PyErr_ExceptionMatches(PyExc_StopIteration))
		{
			PyErr_Clear();
			PyErr_SetString(PyExc_ValueError,"attempt to assign too short sequence");
		}
		goto error;
	}
	else if(j!=slicelength)
	{
		PyErr_SetString(PyExc_ValueError,"attempt to assign too long a sequence");
		goto error;
	}

	ret = 0;

error:
	Py_XDECREF(seq);
	return ret;
}

static int
descr_ass_subscript(GenotypeArrayDescriptorObject *self, PyObject *item, PyObject *value)
{
	if(!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return -1;
	}
	else if(PyIndex_Check(item))
	{
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if(i == -1 && PyErr_Occurred())
			return -1;

		return descr_ass_item(self, i, (UnphasedMarkerModelObject *)value);
	}
	else if(PySlice_Check(item))
		return descr_ass_slice(self, (PySliceObject *)item, value);
	else
	{
		PyErr_SetString(PyExc_TypeError, "indices must be integers or slices");
		return -1;
	}
}

static PyMemberDef descr_members[] = {
	{"_models",	T_OBJECT_EX, offsetof(GenotypeArrayDescriptorObject, models),      RO, "UnphasedMarkerModels"},
	{"offsets",	T_OBJECT_EX, offsetof(GenotypeArrayDescriptorObject, offsets),     RO, "offsets"},
	{"bit_size",	T_UINT,      offsetof(GenotypeArrayDescriptorObject, bit_size),    RO, "bit_size"},
	{"homogeneous",	T_UINT,      offsetof(GenotypeArrayDescriptorObject, homogeneous), RO, "homogeneous"},
	{"max_bit_size",T_UINT,      offsetof(GenotypeArrayDescriptorObject, max_bit_size),RO, "max_bit_size"},
	{"byte_size",	T_UINT,      offsetof(GenotypeArrayDescriptorObject, byte_size),   RO, "byte_size"},
	{NULL}  /* Sentinel */
};

static PyMethodDef descr_methods[] = {
	{NULL}  /* Sentinel */
};

static PyMappingMethods descr_as_mapping = {
	(lenfunc)descr_length,
	(binaryfunc)descr_subscript,
	(objobjargproc)descr_ass_subscript,
};

static PySequenceMethods descr_as_sequence = {
	(lenfunc)descr_length,			/* sq_length */
	0,					/* sq_concat */
	0,					/* sq_repeat */
	(ssizeargfunc)descr_item,		/* sq_item */
	0,					/* sq_slice */
	(ssizeobjargproc)descr_ass_item,	/* sq_ass_item */
	0,					/* sq_ass_slice */
	0,					/* sq_contains */
	0,					/* sq_inplace_concat */
	0,					/* sq_inplace_repeat */
};

PyTypeObject GenotypeArrayDescriptorType = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size           */
	"GenotypeArrayDescriptor",		/* tp_name           */
	sizeof(GenotypeArrayDescriptorObject),	/* tp_basicsize      */
	0,					/* tp_itemsize       */
	(destructor)descr_dealloc,		/* tp_dealloc        */
	0,					/* tp_print          */
	0,					/* tp_getattr        */
	0,					/* tp_setattr        */
	0,					/* tp_compare        */
	0,					/* tp_repr           */
	0,					/* tp_as_number      */
	&descr_as_sequence,			/* tp_as_sequence    */
	&descr_as_mapping,			/* tp_as_mapping     */
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT,			/* tp_flags          */
	"Genotype Array Descriptor objects",	/* tp_doc            */
	(traverseproc)descr_traverse,		/* tp_traverse       */
	(inquiry)descr_clear,			/* tp_clear          */
	0,					/* tp_richcompare    */
	0,					/* tp_weaklistoffset */
	0,					/* tp_iter           */
	0,					/* tp_iternext       */
	descr_methods,				/* tp_methods        */
	descr_members,				/* tp_members        */
	0,					/* tp_getset         */
	0,					/* tp_base           */
	0,					/* tp_dict           */
	0,					/* tp_descr_get      */
	0,					/* tp_descr_set      */
	0,					/* tp_dictoffset     */
	(initproc)descr_init,			/* tp_init           */
	0,					/* tp_alloc          */
	0,					/* tp_new            */
	0,					/* tp_free           */
};

/* UnphasedMarkerModelObject object implementation*/

static int
genomodel_clear(UnphasedMarkerModelObject *self)
{
	Py_CLEAR(self->alleles);
	Py_CLEAR(self->genomap);
	Py_CLEAR(self->genotypes);
	Py_CLEAR(self->genotuples);
	return 0;
}

static void
genomodel_dealloc(UnphasedMarkerModelObject *self)
{
	genomodel_clear(self);
	self->ob_type->tp_free((PyObject *)self);
}

static int
genomodel_traverse(UnphasedMarkerModelObject *self, visitproc visit, void *arg)
{
	/* NOTE: Support for cyclic garbage collection is *NOT* vestigal.
	   Genotypes and models are circularly referential and *will*
	   accumulate if not for the garbage collector. */
	Py_VISIT(self->alleles);
	Py_VISIT(self->genomap);
	Py_VISIT(self->genotypes);
	Py_VISIT(self->genotuples);
	return 0;
}

static PyObject *
genomodel_get_allele(UnphasedMarkerModelObject *self, PyObject *allele)
{
	Py_ssize_t index;

	if(allele != Py_None && !PyString_Check(allele))
	{
		PyErr_SetString(GenotypeRepresentationError,"alleles must be None or strings");
		return NULL;
	}

	index = PySequence_Index(self->alleles, allele);

	if(index==-1)
	{
		PyErr_SetObject(GenotypeLookupError, allele);
		return NULL;
	}
	return PyInt_FromSsize_t(index);
}

static Py_ssize_t
genomodel_add_allele_internal(UnphasedMarkerModelObject *self, PyObject *allele)
{
	Py_ssize_t index;

	assert(!PyErr_Occurred());

	if(allele != Py_None && !PyString_Check(allele))
	{
		PyErr_SetString(GenotypeRepresentationError,"alleles must be None or strings");
		return -1;
	}

	index = PySequence_Index(self->alleles, allele);

	if(index==-1)
	{
		/* Note: alleles[0] is always the missing allele */
		index = PyList_Size(self->alleles);
		if(index==-1) return -1;

		if(index > self->max_alleles)
		{
			PyErr_SetString(GenotypeRepresentationError,"genotype model cannot accomodate additional alleles");
			return -1;
		}

		/* clear index error and append the new allele */
		PyErr_Clear();
		if(PyList_Append(self->alleles, allele) == -1)
			return -1;
	}
	return index;
}

static PyObject *
genomodel_add_allele(UnphasedMarkerModelObject *self, PyObject *allele)
{
	Py_ssize_t result = genomodel_add_allele_internal(self, allele);
	if(result == -1) return NULL;

	return PyInt_FromSsize_t(result);
}

static PyObject *
genomodel_add_genotype(UnphasedMarkerModelObject *self, PyObject *geno)
{
	PyObject *allele1, *allele2, *args, *g;
	Py_ssize_t index1, index2, n;
	int res;

	args = NULL;

	/* If geno is a genotype object */
	if(Genotype_CheckExact(geno))
	{
		GenotypeObject *genoobj = (GenotypeObject *)geno;
		/* If it belongs to this model, just return it */
		if(genoobj->model == self)
			return geno;
		else
		{
			/* Handle foreign genotypes by looking them up by their alleles */
			geno = genotype_alleles(genoobj);
			if(!geno) return NULL;
			Py_DECREF(geno); /* safe to treat as a borrowed reference */
		}

	}

	if(!PyTuple_CheckExact(geno) || PyTuple_GET_SIZE(geno) != 2)
	{
		PyErr_SetString(GenotypeRepresentationError,"genotype must be specified as a 2-tuple");
		return NULL;
	}

	/* Fast path: Genotype already exists */
	g = PyDict_GetItem(self->genomap, geno); /* borrowed ref */
	if(g)
	{
		Py_INCREF(g);
		return g;
	}

	/* Slow path: Add alleles and index the genotype */
	allele1 = PyTuple_GET_ITEM(geno, 0); /* borrowed ref */
	allele2 = PyTuple_GET_ITEM(geno, 1); /* borrowed ref */

	if( !self->allow_hemizygote && ((allele1 == Py_None) ^ (allele2 == Py_None)) )
	{
		PyErr_SetString(GenotypeRepresentationError,"model does not allow hemizgote genotypes");
		return NULL;
	}

	res = PyObject_RichCompareBool(allele1, allele2, Py_GT);
	if(res == -1) return NULL;
	else if(res == 1)
	{
		PyObject *tmp = allele1;
		allele1 = allele2;
		allele2 = tmp;
	}

	index1 = genomodel_add_allele_internal(self,allele1);
	if(index1==-1) return NULL;
	index2 = genomodel_add_allele_internal(self,allele2);
	if(index2==-1) return NULL;

	/* Replace alleles with internal representations to avoid immortalizing new reps */
	allele1 = PyList_GET_ITEM(self->alleles, index1); /* borrowed ref */
	if(!allele1) return NULL;
	allele2 = PyList_GET_ITEM(self->alleles, index2); /* borrowed ref */
	if(!allele2) return NULL;

	/* Create new genotype */
	n = PyList_Size(self->genotypes);
	if(n==-1) return NULL;

	g = PyObject_CallFunction( (PyObject *)&GenotypeType, "(Onnn)", self, index1, index2, n); /* new ref */
	if(!g) return NULL;

	/* Add it to the genotype list */
	if(PyList_Append(self->genotypes, g) == -1)
		goto error;

	/* Create cannonical genotype tuple */
	args = PyTuple_Pack(2, allele1, allele2);
	if(!args) goto error;

	/* Add to to the genotuples list */
	if(PyList_Append(self->genotuples, args) == -1)
	{
		Py_DECREF(args);
		goto error;
	}

	/* Add to the genotype map */
	res = PyDict_SetItem(self->genomap, args, g);
	Py_DECREF(args);
	if(res == -1) goto error;

	/* Create the reversed tuple (since we are an unphased representation) */
	if(index1 != index2)
	{
		args = PyTuple_Pack(2, allele2, allele1);
		if(!args) goto error;
		res = PyDict_SetItem(self->genomap, args, g);
		Py_DECREF(args);
		if(res == -1) goto error;
	}

	/* Return new reference */
	return g;

error:
	Py_XDECREF(g);
	return NULL;
}

static PyObject *
genomodel_get_genotype(UnphasedMarkerModelObject *self, PyObject *geno)
{
	PyObject *g;

	/* If geno is a genotype object */
	if(Genotype_CheckExact(geno))
	{
		GenotypeObject *genoobj = (GenotypeObject *)geno;
		/* If it belongs to this model, just return it */
		if(genoobj->model == self)
		{
			Py_INCREF(geno);
			return geno;
		}
		else
		{
			/* Handle foreign genotypes by looking them up by their alleles */
			geno = genotype_alleles(genoobj);
			if(!geno) return NULL;
			Py_DECREF(geno); /* safe to treat as a borrowed reference */
		}

	}

	if(!PyTuple_CheckExact(geno) || PyTuple_GET_SIZE(geno) != 2)
	{
		PyErr_SetString(GenotypeRepresentationError,"genotype must be specified as a 2-tuple");
		return NULL;
	}

	g = PyDict_GetItem(self->genomap, geno); /* borrowed ref */

	/* FIXME: Implicit model updates need to be reviewed very carefully */

	/* If genotype contains both alleles, then create it implicitly.
	   This is the most obvious behavior and deals with a slew of tricky
	   and inefficient alternative workarounds. */
	if(!g)
	{
		int c1,c2;
		PyObject *allele1, *allele2;

		allele1 = PyTuple_GET_ITEM(geno, 0); /* borrowed ref */
		allele2 = PyTuple_GET_ITEM(geno, 1); /* borrowed ref */

		/* Do not implicitly try to add hemizygotes if they are not allowed */
		if( !self->allow_hemizygote && ((allele1 == Py_None) ^ (allele2 == Py_None)) )
		{
			PyErr_SetObject(GenotypeLookupError, geno);
			return NULL;
		}

		c1 = PySequence_Contains(self->alleles, allele1);
		if(c1 == -1) return NULL;

		c2 = PySequence_Contains(self->alleles, allele2);
		if(c2 == -1) return NULL;

		if(c1 && c2)
			return genomodel_add_genotype(self, geno);
	}

	if(!g)
	{
		PyErr_SetObject(GenotypeLookupError, geno);
		return NULL;
	}

	Py_INCREF(g);
	return g;
}

static int
genomodel_init(UnphasedMarkerModelObject *self, PyObject *args, PyObject *kw)
{
	PyObject *missing, *geno;
	PyObject *allow_hemizygote = Py_False;
	unsigned short n = 2;
	static char *kwlist[] = {"allow_hemizygote","max_alleles", 0};

	genomodel_clear(self);

	if(!PyArg_ParseTupleAndKeywords(args, kw, "|OH:UnphasedMarkerModel", kwlist,
		&allow_hemizygote, &n))
		return -1;

	/* Force models to have at least two alleles */
	if(n < 2)
		n = 2;

	self->allow_hemizygote = PyObject_IsTrue(allow_hemizygote);
	if((short)self->allow_hemizygote == -1) return -1;

	self->max_alleles = n;
	self->genomap     = PyDict_New();
	self->alleles     = PyList_New(0);
	self->genotypes   = PyList_New(0);
	self->genotuples  = PyList_New(0);

	if(!self->genomap || !self->alleles || !self->genotypes || !self->genotuples)
		goto error;

	if(self->allow_hemizygote)
		self->bit_size = ceil(log2((n+1)*(n+2)/2));
	else
		self->bit_size = ceil(log2(n*(n+1)/2 + 1));

	/* Build missing allele and genotype */
	missing = PyTuple_Pack(2, Py_None, Py_None);
	if(!missing) goto error;

	geno = genomodel_add_genotype(self, missing);
	Py_DECREF(missing);
	if(!geno) goto error;
	Py_DECREF(geno);

	return 0;

error:
	genomodel_clear(self);
	return -1;
}

static PyMappingMethods genomodel_as_mapping = {
	(lenfunc)0,
	(binaryfunc)genomodel_get_genotype,
	(objobjargproc)0,
};

static PyMethodDef genomodel_methods[] = {
	{"get_genotype", (PyCFunction)genomodel_get_genotype, METH_O, "get an existing genotype"},
	{"add_genotype", (PyCFunction)genomodel_add_genotype, METH_O, "get or add a genotype"},
	{"get_allele",   (PyCFunction)genomodel_get_allele,   METH_O, "get an allele"},
	{"add_allele",   (PyCFunction)genomodel_add_allele,   METH_O, "get or add an allele"},
	{NULL,		NULL}		/* sentinel */
};

static PyMemberDef genomodel_members[] = {
	{"genomap",          T_OBJECT_EX, offsetof(UnphasedMarkerModelObject, genomap),          RO, "dictionary of alleles to genotype objects"},
	{"genotypes",        T_OBJECT_EX, offsetof(UnphasedMarkerModelObject, genotypes),        RO, "list of genotype objects"},
	{"genotuples",       T_OBJECT_EX, offsetof(UnphasedMarkerModelObject, genotuples),       RO, "list of genotypes tuples"},
	{"alleles",          T_OBJECT_EX, offsetof(UnphasedMarkerModelObject, alleles),          RO, "list of alleles"},
	{"max_alleles",      T_USHORT,    offsetof(UnphasedMarkerModelObject, max_alleles),      RO, "maximum number of alleles"},
	{"allow_hemizygote", T_USHORT,    offsetof(UnphasedMarkerModelObject, allow_hemizygote), RO, "allow_hemizygote"},
	{"bit_size",         T_USHORT,    offsetof(UnphasedMarkerModelObject, bit_size),         RO, "bit size"},
	{NULL}  /* Sentinel */
};


PyDoc_STRVAR(genomodel_doc,
"UnphasedMarkerModelObject(allow_hemizygote=False, max_alleles=2)\n");

PyTypeObject UnphasedMarkerModelType = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size           */
	"UnphasedMarkerModel",                  /* tp_name           */
	sizeof(UnphasedMarkerModelObject),      /* tp_basicsize      */
	0,                                      /* tp_itemsize       */
	(destructor)genomodel_dealloc,		/* tp_dealloc        */
	0,					/* tp_print          */
	0,					/* tp_getattr        */
	0,					/* tp_setattr        */
	0,					/* tp_compare        */
	0,					/* tp_repr           */
	0,					/* tp_as_number      */
	0,					/* tp_as_sequence    */
	&genomodel_as_mapping,			/* tp_as_mapping     */
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,		/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE,	/* tp_flags          */
	genomodel_doc,				/* tp_doc            */
	(traverseproc)genomodel_traverse,	/* tp_traverse       */
	(inquiry)genomodel_clear,		/* tp_clear          */
	0,					/* tp_richcompare    */
	0,					/* tp_weaklistoffset */
	0,					/* tp_iter           */
	0,					/* tp_iternext       */
	genomodel_methods,			/* tp_methods        */
	genomodel_members,			/* tp_members        */
	0,					/* tp_getset         */
	0,					/* tp_base           */
	0,					/* tp_dict           */
	0,					/* tp_descr_get      */
	0,					/* tp_descr_set      */
	0,					/* tp_dictoffset     */
	(initproc)genomodel_init,		/* tp_init           */
	0,					/* tp_alloc          */
	0,					/* tp_new            */
	0,					/* tp_free           */
};

/***************************************/
/* GenotypeArray object implementation */
/***************************************/

static int
genoarray_clear(GenotypeArrayObject *self)
{
	Py_CLEAR(self->descriptor);
	return 0;
}

static void
genoarray_dealloc(GenotypeArrayObject *self)
{
	genoarray_clear(self);
	self->ob_type->tp_free((PyObject *)self);
}

static Py_ssize_t
genoarray_length(GenotypeArrayObject *self)
{
	Py_ssize_t n;
	PyObject *models;

	if(!self->descriptor || !GenotypeArrayDescriptor_CheckExact(self->descriptor))
	{
		PyErr_SetString(PyExc_TypeError,"invalid descriptor");
		return -1;
	}

	models = self->descriptor->models;
	if(!models) return -1;

	n = PyList_Size(models);
	if(n == -1 && PyErr_Occurred())
		return -1;

	return n;
}

static PyObject *
genoarray_repr(GenotypeArrayObject *self)
{
	PyObject *junk, *result;

	junk = PySequence_GetSlice( (PyObject *)self, 0, genoarray_length(self));
	if(!junk) return NULL;

	result = PyObject_Repr(junk);
	Py_DECREF(junk);
	return result;
}

static int
genoarray_traverse(GenotypeArrayObject *self, visitproc visit, void *arg)
{
	Py_VISIT(self->descriptor);
	return 0;
}

static PyObject *
genoarray_alloc(PyTypeObject *type, Py_ssize_t nitems)
{
	PyObject *obj;
	const size_t size = type->tp_basicsize + nitems;

	if(PyType_IS_GC(type))
		obj = _PyObject_GC_Malloc(size);
	else
		obj = (PyObject *)PyObject_MALLOC(size);

	if(obj == NULL)
		return PyErr_NoMemory();

	memset(obj, '\0', size);

	if(type->tp_flags & Py_TPFLAGS_HEAPTYPE)
		Py_INCREF(type);

	PyObject_INIT(obj, type);

	if(PyType_IS_GC(type))
		PyObject_GC_Track(obj);
	return obj;
}

static PyObject *
genoarray_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
	unsigned int n;
	GenotypeArrayObject *self;
	PyObject *descr;
	GenotypeArrayDescriptorObject *descriptor = NULL;;
	PyObject *genos = NULL;
	static char *kwlist[] = {"descriptor","genos", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|O:GenotypeArray", kwlist, &descr, &genos))
		return NULL;

	if(GenotypeArrayDescriptor_CheckExact(descr))
		descriptor = (GenotypeArrayDescriptorObject *)descr;
	else if(GenotypeArray_Check(descr))
		descriptor = ((GenotypeArrayObject *)descr)->descriptor;
	else
	{
		PyErr_SetString(PyExc_TypeError,"invalid descriptor object");
		return NULL;
	}

	n = descriptor->byte_size;
	self = (GenotypeArrayObject *)type->tp_alloc(&GenotypeArrayType, n);
	if(!self) return NULL;

	self->descriptor = descriptor;
	Py_INCREF(self->descriptor);

	if( genos && genos!=Py_None )
	{
		if( PySequence_SetSlice( (PyObject *)self, 0, genoarray_length(self), genos) == -1)
		{
			Py_DECREF(self);
			return NULL;
		}
	}
	return (PyObject *)self;
}

static long
genoarray_nohash(PyObject *self)
{
	PyErr_SetString(PyExc_TypeError, "genoarray objects are unhashable");
	return -1;
}

static inline PyObject *
genoarray_inner_get(PyObject *models, const unsigned char *data, Py_ssize_t datasize,
                    unsigned int *offsets, Py_ssize_t i)
{
	Py_ssize_t k, offset1, offset2;
	PyObject *geno;
	UnphasedMarkerModelObject *model;
	char *status;

	model = (UnphasedMarkerModelObject *)PyList_GetItem(models, i); /* borrowed ref */

	if(!model)
		return NULL;

	if(!UnphasedMarkerModel_Check(model) || !model->genotypes)
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
			return NULL;
	}

	offset1 = offsets[i];
	offset2 = offsets[i+1];

	k = bitarray_getbits(data, datasize, offset1, offset2-offset1, &status);

	if(status)
	{
		PyErr_SetString(PyExc_IndexError, status);
		return NULL;
	}

	geno = PyList_GetItem(model->genotypes, k); /* borrowed ref */
	Py_XINCREF(geno);

	return geno;
}

static inline int
genoarray_checkstate(GenotypeArrayObject *self)
{
	if(!self->data || !self->descriptor
	               || !GenotypeArrayDescriptor_CheckExact(self->descriptor))
	{
		PyErr_SetString(PyExc_ValueError,"invalid genoarray state");
		return -1;
	}
	return 0;
}

static PyObject *
genoarray_item(GenotypeArrayObject *self, Py_ssize_t item)
{
	Py_ssize_t n, datasize;
	PyObject *models;
	unsigned int *offsets;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	offsets  = (unsigned int *)PyArray_DATA(self->descriptor->offsets);
	models   = self->descriptor->models;
	/* FIXME: Set error state */

	if(!PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		return NULL;
	}

	n = PyList_GET_SIZE(models);

	if( n+1 != PyArray_SIZE(self->descriptor->offsets) )
	{
		PyErr_SetString(PyExc_ValueError,"offsets and models sizes must agree");
		return NULL;
	}

	if(item < 0 || item >= n)
	{
		PyErr_SetString(PyExc_IndexError,"genotype index out of range");
		return NULL;
	}

	datasize = self->descriptor->byte_size;
	return genoarray_inner_get(models, self->data, datasize, offsets, item);
}

static PyObject *
genoarray_slice(GenotypeArrayObject *self, PySliceObject *slice)
{
	Py_ssize_t i, j, n, datasize;
	Py_ssize_t start, stop, step, slicelength;
	PyObject *result, *models, *geno;
	unsigned int *offsets;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	result   = NULL;
	datasize = self->descriptor->byte_size;

	offsets  = (unsigned int *)PyArray_DATA(self->descriptor->offsets);
	models   = self->descriptor->models;
	if(!models || !offsets) goto done;

	if(!PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		goto done;
	}

	n = PyList_GET_SIZE(models);

	if( n+1 != PyArray_SIZE(self->descriptor->offsets) )
	{
		PyErr_SetString(PyExc_ValueError,"offsets and models sizes must agree");
		goto done;
	}

	if(PySlice_GetIndicesEx(slice, n, &start, &stop, &step, &slicelength) < 0)
		goto done;

	result = PyList_New(slicelength);
	if(!result) goto done;

	for (i=start, j=0; j<slicelength; i+=step, j++)
	{
		geno = genoarray_inner_get(models, self->data, datasize, offsets, i);
		if(!geno) goto error;
		PyList_SET_ITEM(result, j, geno);
	}
	goto done;

error:
	Py_XDECREF(result);
	result = NULL;

done:
	return result;
}

static PyObject *
genoarray_pick(GenotypeArrayObject *self, PyObject *indices)
{
	Py_ssize_t i, j, n, datasize, len;
	PyObject *result=NULL, *models, *geno, *indseq=NULL;
	unsigned int *offsets;
	PyObject **indexitems;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	datasize = self->descriptor->byte_size;
	offsets  = (unsigned int *)PyArray_DATA(self->descriptor->offsets);
	models   = self->descriptor->models;
	if(!models || !offsets) goto error;

	if(!PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		goto error;
	}

	n = PyList_GET_SIZE(models);

	if( n+1 != PyArray_SIZE(self->descriptor->offsets) )
	{
		PyErr_SetString(PyExc_ValueError,"offsets and models sizes must agree");
		goto error;
	}

	assert(indices);
	indseq = PySequence_Fast(indices,"cannot convert indices into a sequence");
	if(!indseq) goto error;

	len = PySequence_Fast_GET_SIZE(indseq);
	indexitems = PySequence_Fast_ITEMS(indseq);

	result = PyList_New(len);
	if(!result) goto error;

	for (i=0; i<len; i++)
	{
		j = PyNumber_AsSsize_t(indexitems[i], PyExc_IndexError);

		if(j==-1 && PyErr_Occurred()) goto error;
		if(j<0 || j>=n)
		{
			PyErr_SetString(PyExc_IndexError,"genotype index out of range");
			goto error;
		}

		geno = genoarray_inner_get(models, self->data, datasize, offsets, j);
		if(!geno) goto error;

		PyList_SET_ITEM(result, i, geno);
	}

	Py_DECREF(indseq);
	return result;

error:
	Py_XDECREF(indseq);
	Py_XDECREF(result);
	return NULL;
}

static PyObject *
genoarray_subscript(GenotypeArrayObject *self, PyObject *item)
{
	if(PyIndex_Check(item))
	{
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if(i == -1 && PyErr_Occurred())
			return NULL;
		return genoarray_item(self, i);
	}
	else if(PySlice_Check(item))
		return genoarray_slice(self, (PySliceObject *)item);
	else if(PySequence_Check(item))
		return genoarray_pick(self, item);
	else
	{
		PyErr_SetString(PyExc_TypeError, "indices must be integers, slices, or lists of integers");
		return NULL;
	}
}

static int
genoarray_inner_set(PyObject *models, PyObject *geno, unsigned char *data, Py_ssize_t datasize, unsigned int *offsets, Py_ssize_t i)
{
	UnphasedMarkerModelObject *model;
	GenotypeObject *g;
	Py_ssize_t offset1, offset2;
	char *status;

	model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(models, i); /* borrowed ref */

	if(!model || !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return -1;
	}

	if(!Genotype_CheckExact(geno) ||
	   (Genotype_CheckExact(geno) && ((GenotypeObject *)geno)->model != model))
	{
		/* Handle foreign genotypes by looking them up by their alleles */
		geno = genomodel_get_genotype( (UnphasedMarkerModelObject *)model, geno);
		if(!geno) return -1;
		Py_DECREF(geno); /* safe to treat as a borrowed reference */
	}

	g = (GenotypeObject *)geno;

	offset1 = offsets[i];
	offset2 = offsets[i+1];

	if(g->index >= (1<<(offset2-offset1+1)))
	{
		PyErr_SetString(PyExc_ValueError,"genotype index out of range");
		return -1;
	}

	bitarray_setbits(data, datasize, offset1, g->index, offset2-offset1, &status);

	if(status)
	{
		PyErr_SetString(PyExc_IndexError, status);
		return -1;
	}
	return 0;
}

static int
genoarray_ass_item(GenotypeArrayObject *self, Py_ssize_t item, PyObject *value)
{
	Py_ssize_t n, datasize;
	PyObject *models;
	unsigned int *offsets;

	if(!value)
	{
		/* delete item */
		PyErr_SetString(PyExc_TypeError, "genoarray objects do not support item deletion");
		return -1;
	}

	if( genoarray_checkstate(self) == -1 )
		return -1;

	datasize = self->descriptor->byte_size;

	offsets  = (unsigned int *)PyArray_DATA(self->descriptor->offsets);
	models   = self->descriptor->models;
	if(!models || !offsets) return -1;

	if(!PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		return -1;
	}

	n = PyList_GET_SIZE(models);

	if( n+1 != PyArray_SIZE(self->descriptor->offsets) )
	{
		PyErr_SetString(PyExc_ValueError,"offsets and models sizes must agree");
		return -1;
	}

	if(item < 0 || item >= n)
	{
		PyErr_SetString(PyExc_IndexError,"genotype index out of range");
		return -1;
	}

	if( genoarray_inner_set(models, value, self->data, datasize, offsets, item) == -1 )
		return -1;

	return 0;
}

static int
genoarray_ass_slice(GenotypeArrayObject *self, PySliceObject *slice, PyObject *value)
{
	Py_ssize_t i, j, n, m, datasize;
	Py_ssize_t start, stop, step, slicelength;
	PyObject *models, *seq;
	PyObject **seqitems;
	unsigned int *offsets;
	int ret = -1;

	if(!value)
	{
		/* delete slice */
		PyErr_SetString(PyExc_TypeError, "genoarray objects do not support item deletion");
		return -1;
	}

	if( genoarray_checkstate(self) == -1 )
		return -1;

	seq = NULL;
	datasize = self->descriptor->byte_size;

	offsets  = (unsigned int *)PyArray_DATA(self->descriptor->offsets);
	models   = self->descriptor->models;
	if(!models || !offsets) goto error;

	if(!PyList_Check(models))
	{
		PyErr_SetString(PyExc_TypeError,"models must be a list");
		goto error;
	}

	n = PyList_GET_SIZE(models);

	if( n+1 != PyArray_SIZE(self->descriptor->offsets) )
	{
		PyErr_SetString(PyExc_ValueError,"offsets and models sizes must agree");
		goto error;
	}

	if(PySlice_GetIndicesEx(slice, n, &start, &stop, &step, &slicelength) < 0)
		goto error;

	if( PyList_CheckExact(value) || PyTuple_CheckExact(value) || (PyObject *)self == value )
	{
		seq = PySequence_Fast(value,"cannot convert slice into a sequence");
		if(!seq) goto error;

		m = PySequence_Fast_GET_SIZE(seq);

		if(m != slicelength)
		{
			PyErr_Format(PyExc_ValueError, "attempt to assign sequence of size %zd to slice of size %zd",
			             m, slicelength);
			goto error;
		}

		if(!slicelength) goto done;
		seqitems = PySequence_Fast_ITEMS(seq);

		for (i=start, j=0; j<slicelength; i+=step, j++)
		{
			if( genoarray_inner_set(models, seqitems[j], self->data, datasize, offsets, i) == -1 )
				goto error;
		}
	}
	else
	{
		PyObject *(*iternext)(PyObject *);

		seq = PyObject_GetIter(value);
		if(seq == NULL)
			goto error;
		iternext = *seq->ob_type->tp_iternext;

		for (i=start, j=0; j<slicelength; i+=step, j++)
		{
			PyObject *geno = iternext(seq);
			if(!geno)
			{
				if(PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_StopIteration))
				{
					PyErr_Clear();
					PyErr_SetString(PyExc_ValueError,"attempt to assign too short sequence");
				}
				goto error;
			}
			if( genoarray_inner_set(models, geno, self->data, datasize, offsets, i) == -1 )
			{
				Py_DECREF(geno);
				goto error;
			}
			Py_DECREF(geno);
		}
		if(j!=slicelength)
		{
			PyErr_SetString(PyExc_ValueError,"attempt to assign too long a sequence");
			goto error;
		}
	}
done:	ret = 0;

error:
	Py_XDECREF(seq);
	return ret;
}

static int
genoarray_ass_subscript(GenotypeArrayObject *self, PyObject *item, PyObject *value)
{
	if(PyIndex_Check(item))
	{
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if(i == -1 && PyErr_Occurred())
			return -1;

		return genoarray_ass_item(self, i, value);
	}
	else if(PySlice_Check(item))
		return genoarray_ass_slice(self, (PySliceObject *)item, value);
	else
	{
		PyErr_SetString(PyExc_TypeError, "indices must be integers or slices");
		return -1;
	}
}

static PyObject *
genoarray_descriptor_get(GenotypeArrayObject *self, void *closure)
{
	Py_INCREF(self->descriptor);
	return (PyObject *)self->descriptor;
}

static PyObject *
genoarray_data_get(GenotypeArrayObject *self, void *closure)
{
	PyArrayObject *data;
	npy_intp dims;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	dims = self->descriptor->byte_size;
	data = (PyArrayObject *)PyArray_SimpleNewFromData(1,&dims,NPY_UBYTE,self->data);
	if(!data) return NULL;

	data->base = (PyObject *)self;
	Py_INCREF(self);

	return (PyObject *)data;
}

static int
genoarray_data_set(GenotypeArrayObject *self, PyObject *new_data, void *closure)
{
	if( genoarray_checkstate(self) == -1 )
		return -1;

	if(!new_data)
	{
		PyErr_SetString(PyExc_ValueError,"cannot delete genoarray data");
		return -1;
	}

	if(!PyArray_CheckExact(new_data))
	{
		PyErr_SetString(PyExc_TypeError,"new data must be ndarray type");
		return -1;
	}

	if(PyArray_NBYTES(new_data) != self->descriptor->byte_size)
	{
		PyErr_Format(PyExc_ValueError,"new data must be same size as current (%zd != %d)",
		                               PyArray_NBYTES(new_data), self->descriptor->byte_size);
		return -1;
	}

	memcpy(self->data, PyArray_DATA(new_data), self->descriptor->byte_size);
	return 0;
}

static Py_ssize_t
genoarray_getreadbuf(GenotypeArrayObject *self, Py_ssize_t index, const void **ptr)
{
	if( index != 0 )
	{
		PyErr_SetString(PyExc_SystemError,
				"accessing non-existent genoarray segment");
		return -1;
	}
	*ptr = (void *)self->data;
	return self->descriptor->byte_size;
}

static Py_ssize_t
genoarray_getwritebuf(GenotypeArrayObject *self, Py_ssize_t index, const void **ptr)
{
	PyErr_SetString(PyExc_TypeError,
			"Cannot use genoarray as modifiable buffer");
	return -1;
}

static Py_ssize_t
genoarray_getsegcount(GenotypeArrayObject *self, Py_ssize_t *lenp)
{
	/* Return one segment and the length in bytes if lenp is not NULL */
	if(lenp)
		*lenp = self->descriptor->byte_size;
	return 1;
}

static Py_ssize_t
genoarray_getcharbuf(GenotypeArrayObject *self, Py_ssize_t index, const char **ptr)
{
	PyErr_SetString(PyExc_TypeError,
			"Cannot use genoarray as a character buffer");
	return -1;
}

/******************************************************************************************************/

typedef struct {
	unsigned int *indices;
	UnphasedMarkerModelObject *model;
} genotype_indices_state;

static int indices_foreach(Py_ssize_t i, GenotypeObject *geno, genotype_indices_state *state)
{
	if(state->model && state->model != geno->model)
	{
		PyErr_SetString(PyExc_TypeError,"heterogeneous genotype model");
		return -1;
	}
	state->indices[i] = geno->index;
	return 0;
}

static PyObject *
genotype_indices(PyObject *genos, UnphasedMarkerModelObject *model)
{
	genotype_indices_state state;
	PyObject *index_array;
	npy_intp len;

	len = PyObject_Size(genos);
	if(len==-1) return NULL;

	if(model && !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return NULL;
	}

	index_array = PyArray_SimpleNew(1,&len,NPY_UINT);
	if(!index_array) return NULL;

	state.model = model;
	state.indices = (unsigned int *)PyArray_DATA(index_array);

	if(for_each_genotype(genos, (geno_foreach)indices_foreach, &state) < 0)
	{
		Py_DECREF(index_array);
		return NULL;
	}
	return index_array;
}

static PyObject *
genoarray_indices_method(GenotypeArrayObject *self, PyObject *args, PyObject *kw)
{
	UnphasedMarkerModelObject *model = NULL;
	static char *kwlist[] = {"model", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "|O:indices", kwlist, &model))
		return NULL;

	return genotype_indices( (PyObject *)self, model);
}

static PyObject *
genotype_indices_func(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *genos, *model=NULL;
	static char *kwlist[] = {"genos", "model", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|O:genotype_indices", kwlist, &genos, &model))
		return NULL;

	return genotype_indices(genos, (UnphasedMarkerModelObject *)model);
}

/******************************************************************************************************/

typedef struct {
	unsigned long *counts;
	UnphasedMarkerModelObject *model;
} genotype_counts_state;

static int counts_foreach(Py_ssize_t i, GenotypeObject *geno, genotype_counts_state *state)
{
	if(state->model && state->model != geno->model)
	{
		PyErr_SetString(PyExc_TypeError,"heterogeneous genotype model");
		return -1;
	}
	state->counts[geno->index] += 1;
	return 0;
}

static PyObject *
count_genotypes(PyObject *genos, PyObject *count_array, PyObject *model)
{
	genotype_counts_state state;
	npy_intp model_len;
	int len;

	if(count_array==Py_None) count_array=NULL;
	if(model==Py_None)       model=NULL;

	len = PyObject_Size(genos);
	if(len==-1) return NULL;

	if(!len && !model)
	{
		PyErr_SetString(GenotypeRepresentationError,
		    "empty genotype sequence with unknown model");
		return NULL;
	}

	if(len && !model && GenotypeArray_Check(genos))
	{
		GenotypeArrayObject *garray = (GenotypeArrayObject *)genos;
		model = PyList_GetItem(garray->descriptor->models, 0); /* borrowed ref */
	}

	if(len && !model)
	{
		GenotypeObject *geno = (GenotypeObject *)PySequence_GetItem(genos, 0);
		if(!geno || !Genotype_CheckExact(geno))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos at index %zd", 0L);
			Py_XDECREF(geno);
			return NULL;
		}
		model = (PyObject *)geno->model; /* borrowed ref */
		Py_DECREF(geno);
	}

	if(!model || !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return NULL;
	}

	model_len = PyList_Size( ((UnphasedMarkerModelObject *)model)->genotypes);

	if(count_array)
	{
		if(!CountArray_Check(count_array,model_len))
		{
			PyErr_SetString(PyExc_ValueError,"invalid count array");
			return NULL;
		}
		Py_INCREF(count_array);
	}
	else
	{
		count_array = PyArray_SimpleNew(1,&model_len,NPY_LONG);
		if(!count_array) return NULL;
		PyArray_FILLWBYTE(count_array, 0);
	}

	state.model  = (UnphasedMarkerModelObject *)model; /* borrowed ref */
	state.counts = (unsigned long *)PyArray_DATA(count_array); /* borrowed ref */

	if(for_each_genotype(genos, (geno_foreach)counts_foreach, &state) < 0)
	{
		Py_DECREF(count_array);
		return NULL;
	}
	return count_array;
}

static PyObject *
genoarray_counts_method(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *model=NULL, *counts=NULL;
	static char *kwlist[] = {"counts", "model", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "|OO:counts", kwlist, &counts, &model))
		return NULL;

	return count_genotypes(self, model, counts);
}

static PyObject *
count_genotypes_func(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *genos, *model=NULL, *counts=NULL;
	static char *kwlist[] = {"genos", "counts", "model", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|OO:count_genotypes", kwlist, &genos, &counts, &model))
		return NULL;

	return count_genotypes(genos, counts, model);
}

/******************************************************************************************************/

static int tolist_foreach(Py_ssize_t i, GenotypeObject *geno, void *state)
{
	Py_INCREF(geno);
	PyList_SET_ITEM( (PyObject *)state, i, (PyObject *)geno);
	return 0;
}

static PyObject *
genoarray_tolist(GenotypeArrayObject *self)
{
	PyObject *result;
	int len;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	len = genoarray_length(self);
	if(len==-1) return NULL;

	result = PyList_New(len);
	if(!result) return NULL;

	if(for_each_genotype( (PyObject *)self, tolist_foreach, result) < 0)
	{
		Py_DECREF(result);
		return NULL;
	}
	return result;
}

/******************************************************************************************************/

static int category_foreach(Py_ssize_t i, GenotypeObject *geno, void *state)
{
	unsigned long *counts = state;
	counts[ genotype_category(geno) ] += 1;
	return 0;
}

static PyObject *
genotype_categories(PyObject *genos, PyObject *count_array)
{
	npy_intp count_len=4;
	int len;

	if(count_array==Py_None) count_array=NULL;

	len = PyObject_Size(genos);
	if(len==-1) return NULL;

	if(count_array)
	{
		if(!CountArray_Check(count_array,count_len))
		{
			PyErr_SetString(PyExc_ValueError,"invalid count array");
			return NULL;
		}
		Py_INCREF(count_array);
	}
	else
	{
		count_array = PyArray_SimpleNew(1,&count_len,NPY_LONG);
		if(!count_array) return NULL;
		PyArray_FILLWBYTE(count_array, 0);
	}

	if(for_each_genotype(genos, category_foreach, PyArray_DATA(count_array)) < 0)
	{
		Py_DECREF(count_array);
		return NULL;
	}
	return count_array;
}

static PyObject *
genoarray_categories_method(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *counts=NULL;
	static char *kwlist[] = {"counts", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "|O:counts", kwlist, &counts))
		return NULL;

	return genotype_categories(self, counts);
}

static PyObject *
genotypes_categories_func(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *genos, *counts=NULL;
	static char *kwlist[] = {"genos", "counts", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|O:count_genotypes", kwlist, &genos, &counts))
		return NULL;

	return genotype_categories(genos, counts);
}

/******************************************************************************************************/

static PyMappingMethods genoarray_as_mapping = {
	(lenfunc)genoarray_length,
	(binaryfunc)genoarray_subscript,
	(objobjargproc)genoarray_ass_subscript,
};

static PySequenceMethods genoarray_as_sequence = {
	(lenfunc)genoarray_length,			/* sq_length */
	0,						/* sq_concat */
	0,						/* sq_repeat */
	(ssizeargfunc)genoarray_item,			/* sq_item */
	0,						/* sq_slice */
	(ssizeobjargproc)genoarray_ass_item,		/* sq_ass_item */
	0,						/* sq_ass_slice */
	0,						/* sq_contains */
	0,						/* sq_inplace_concat */
	0,						/* sq_inplace_repeat */
};

static PyBufferProcs genoarray_as_buffer = {
	(readbufferproc)genoarray_getreadbuf,
	(writebufferproc)genoarray_getwritebuf,
	(segcountproc)genoarray_getsegcount,
	(charbufferproc)genoarray_getcharbuf,
};

static PyGetSetDef genoarray_getsetlist[] = {
	{"descriptor",
	 (getter)genoarray_descriptor_get,
	 NULL,
	 NULL},
	{"data",
	 (getter)genoarray_data_get,
	 (setter)genoarray_data_set,
	 NULL},
	{NULL, NULL, NULL, NULL},  /* Sentinel */
};

PyDoc_STRVAR(getitem_doc,
"x.__getitem__(y) <==> x[y]");

PyDoc_STRVAR(genoarray_doc,
"genoarray(descriptor) -> new empty genotype array\n"
"genoarray(descriptor, genos) -> new genotype array initialized with genos");

static PyMethodDef genoarray_methods[] = {
	{"indices", (PyCFunction)genoarray_indices_method, METH_KEYWORDS,
		"Return an array of integer genotype indices"},
	{"counts", (PyCFunction)genoarray_counts_method, METH_KEYWORDS,
		"Return an array of genotype counts by index"},
	{"categories", (PyCFunction)genoarray_categories_method, METH_KEYWORDS,
		"Return an array of counts by genotype category"},
	{"tolist", (PyCFunction)genoarray_tolist, METH_NOARGS,
		"Return a list of genotypes"},
	{"__getitem__", (PyCFunction)genoarray_subscript, METH_O|METH_COEXIST, getitem_doc},
	{NULL,		NULL}		/* sentinel */
};

PyTypeObject GenotypeArrayType = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size           */
	"GenotypeArray",                        /* tp_name           */
	sizeof(GenotypeArrayObject),            /* tp_basicsize      */
	0,					/* tp_itemsize       */
	(destructor)genoarray_dealloc,		/* tp_dealloc        */
	0,					/* tp_print          */
	0,					/* tp_getattr        */
	0,					/* tp_setattr        */
	0,					/* tp_compare        */
	(reprfunc)genoarray_repr,		/* tp_repr           */
	0,					/* tp_as_number      */
	&genoarray_as_sequence,			/* tp_as_sequence    */
	&genoarray_as_mapping,			/* tp_as_mapping     */
	genoarray_nohash,			/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	&genoarray_as_buffer,			/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT,			/* tp_flags          */
	genoarray_doc,				/* tp_doc            */
	(traverseproc)genoarray_traverse,	/* tp_traverse       */
	(inquiry)genoarray_clear,		/* tp_clear          */
	0,					/* tp_richcompare    */
	0,					/* tp_weaklistoffset */
	0,					/* tp_iter           */
	0,					/* tp_iternext       */
	genoarray_methods,			/* tp_methods        */
	0,					/* tp_members        */
	genoarray_getsetlist,			/* tp_getset         */
	0,					/* tp_base           */
	0,					/* tp_dict           */
	0,					/* tp_descr_get      */
	0,					/* tp_descr_set      */
	0,					/* tp_dictoffset     */
	0,					/* tp_init           */
	genoarray_alloc,			/* tp_alloc          */
	genoarray_new,				/* tp_new            */
	0,					/* tp_free           */
};

int
for_each_genotype_genoarray(GenotypeArrayObject *genos, geno_foreach func, void *state)
{
	UnphasedMarkerModelObject *model;
	GenotypeObject *geno;
	Py_ssize_t i, j, datasize, offset1, offset2;
	unsigned int *offsets;
	char *status;
	int len;

	if( genoarray_checkstate(genos) == -1 )
		return -1;

	len = genoarray_length(genos);
	if(len==-1) return -1;

	offsets = (unsigned int *)PyArray_DATA(genos->descriptor->offsets);

	if(offsets[0] == 0 && genos->descriptor->homogeneous == 8)
	{
		for(i=0; i < len; ++i)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, i);
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, genos->data[i]);
			if(!geno) goto error;
			if( func(i,geno,state) < 0) goto error;
		}
	}
	else if(offsets[0] == 0 && genos->descriptor->homogeneous == 4)
	{
		for(i=0,j=0; i < len/2; ++i)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>4)&0x0F);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;

			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes,  genos->data[i]    &0x0F);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
		if(len&1)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>4)&0x0F);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
	}
	else if(offsets[0] == 0 && genos->descriptor->homogeneous == 2)
	{
		for(i=0,j=0; i < len/4; ++i)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>6)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;

			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>4)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;

			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>2)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;

			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes,  genos->data[i]    &0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
		len = len&3;
		if(len>0)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>6)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
		if(len>1)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>4)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
		if(len>2)
		{
			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, j);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, (genos->data[i]>>2)&0x03);
			if(!geno) goto error;
			if( func(j++,geno,state) < 0) goto error;
		}
	}
	else
	{
		datasize = genos->descriptor->byte_size;
		for(i=0; i < len; ++i)
		{
			offset1 = offsets[i];
			offset2 = offsets[i+1];

			j = bitarray_getbits(genos->data, datasize, offset1, offset2-offset1, &status);

			if(status)
			{
				PyErr_SetString(PyExc_IndexError, status);
				goto error;
			}

			model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(genos->descriptor->models, i);
			if(!model) goto error;
			geno  = (GenotypeObject *)PyList_GetItem(model->genotypes, j);
			if(!geno) goto error;
			if( func(i,geno,state) < 0) goto error;
		}
	}
	return 0;

error:
	return -1;
}

int
for_each_genotype(PyObject *genos, geno_foreach func, void *state)
{
	PyObject **items;
	int len;
	Py_ssize_t i;

	if( GenotypeArray_Check(genos) )
		return for_each_genotype_genoarray( (GenotypeArrayObject *)genos, func, state);

	genos = PySequence_Fast(genos,"genos is not a sequence");
	if(!genos) return -1;

	len   = PySequence_Fast_GET_SIZE(genos);
	items = PySequence_Fast_ITEMS(genos);

	for(i = 0; i < len; ++i)
	{
		GenotypeObject *geno = (GenotypeObject *)items[i];

		if(!geno || !Genotype_CheckExact(geno))
		{
			PyErr_Format(GenotypeRepresentationError,
		             "invalid genotype object in genos at index %zd", i);
			goto error;
		}
		if( func(i,geno,state) < 0) goto error;
	}

	Py_DECREF(genos);
	return 0;

error:
	Py_DECREF(genos);
	return -1;
}

/******************************************************************************************************/

typedef struct {
	PyObject *sample_counts;
	unsigned long *locus_counts;
	UnphasedMarkerModelObject *model;
} locus_summary_state;

static int locus_summary_foreach(Py_ssize_t i, GenotypeObject *geno, locus_summary_state *state)
{
	if(state->model && state->model != geno->model)
	{
		PyErr_SetString(PyExc_TypeError,"heterogeneous genotype model");
		return -1;
	}
	state->locus_counts[geno->index] += 1;
	*(unsigned long *)PyArray_GETPTR2(state->sample_counts, i, genotype_category(geno)) += 1;
	return 0;
}

static PyObject *
locus_summary(PyObject *genos, PyObject *sample_counts, PyObject *locus_counts, PyObject *model)
{
	locus_summary_state state;
	PyObject *ret = NULL;
	npy_intp model_len;
	Py_ssize_t geno_len;
	npy_intp category_len[2];
	int len;

	if(sample_counts==Py_None) sample_counts = NULL;
	if(locus_counts==Py_None)  locus_counts  = NULL;
	if(model==Py_None)         model         = NULL;

	Py_XINCREF(sample_counts);
	Py_XINCREF(locus_counts);

	len = PyObject_Size(genos);
	if(len==-1) goto error;

	if(!len && !model)
	{
		PyErr_SetString(GenotypeRepresentationError,
		             "empty genotype sequence with unknown model");
		goto error;
	}

	if(len && !model && GenotypeArray_Check(genos))
	{
		GenotypeArrayObject *garray = (GenotypeArrayObject *)genos;
		model = PyList_GetItem(garray->descriptor->models, 0); /* borrowed ref */
	}
	else if(len && !model)
	{
		GenotypeObject *geno = (GenotypeObject *)PySequence_GetItem(genos, 0); /* new ref */
		if(!geno || !Genotype_CheckExact(geno))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos at index %zd", 0L);
			Py_XDECREF(geno);
			goto error;
		}
		model = (PyObject *)geno->model; /* borrowed ref */
		Py_DECREF(geno);
	}

	if(!model || !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		goto error;
	}

	// Count array is conservative and defaults to 1<<model->bit_size.
	model_len = 1<<((UnphasedMarkerModelObject *)model)->bit_size;

	// However a pre-existing count array must only be longer than the genotypes
	geno_len  = PyList_Size( ((UnphasedMarkerModelObject *)model)->genotypes);
	if(geno_len==-1) goto error;

	if(locus_counts)
	{
		if(!PyArray_Check(locus_counts)     || PyArray_TYPE(locus_counts)!=PyArray_LONG ||
		    PyArray_NDIM(locus_counts) != 1 || PyArray_DIMS(locus_counts)[0]<geno_len)
		{
			PyErr_SetString(PyExc_ValueError,"invalid locus count array");
			goto error;
		}
	}
	else
	{
		locus_counts = PyArray_SimpleNew(1,&model_len,NPY_LONG);
		if(!locus_counts) goto error;
		PyArray_FILLWBYTE(locus_counts, 0);
	}

	if(sample_counts)
	{
		if(!CountArray_Check2(sample_counts,len,4))
		{
			PyErr_SetString(PyExc_ValueError,"invalid sample count array");
			goto error;
		}
	}
	else
	{
		category_len[0] = len;
		category_len[1] = 4;
		sample_counts = PyArray_SimpleNew(2,category_len,NPY_LONG);
		if(!sample_counts) goto error;
		PyArray_FILLWBYTE(sample_counts, 0);
	}

	state.model         = (UnphasedMarkerModelObject *)model;
	state.locus_counts  = (unsigned long *)PyArray_DATA(locus_counts);
	state.sample_counts = sample_counts;

	if(for_each_genotype(genos, (geno_foreach)locus_summary_foreach, &state) < 0)
		goto error;

	ret = PyTuple_Pack(2, locus_counts, sample_counts);

error:
	Py_XDECREF(locus_counts);
	Py_XDECREF(sample_counts);

	return ret;
}

static PyObject *
locus_summary_func(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *genos, *sample_counts=NULL, *locus_counts=NULL, *model=NULL;
	static char *kwlist[] = {"genos", "sample_counts", "locus_counts", "model", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|OOO:locus_summary", kwlist,
	                                &genos, &sample_counts, &locus_counts, &model))
		return NULL;

	return locus_summary(genos, sample_counts, locus_counts, model);
}

/******************************************************************************************************/

typedef struct {
	unsigned long max_genos;
	unsigned long *sample_counts;
	PyObject *locus_counts;
} sample_summary_state;

static int sample_summary_foreach(Py_ssize_t i, GenotypeObject *geno, sample_summary_state *state)
{
	if(geno->index >= state->max_genos)
	{
		PyErr_SetString(PyExc_ValueError,"invalid locus count array");
		return -1;
	}
	state->sample_counts[genotype_category(geno)] += 1;
	*(unsigned long *)PyArray_GETPTR2(state->locus_counts, i, geno->index) += 1;
	return 0;
}

static PyObject *
sample_summary(PyObject *genos, PyObject *locus_counts, PyObject *sample_counts)
{
	sample_summary_state state;
	PyObject *ret = NULL;
	npy_intp category_len=4;
	npy_intp locus_len[2];
	int len;

	if(locus_counts==Py_None)  locus_counts  = NULL;
	if(sample_counts==Py_None) sample_counts = NULL;

	Py_XINCREF(locus_counts);
	Py_XINCREF(sample_counts);

	len = PyObject_Size(genos);
	if(len==-1) goto error;

	if(sample_counts)
	{
		if(!CountArray_Check(sample_counts,4))
		{
			PyErr_SetString(PyExc_ValueError,"invalid sample count array");
			goto error;
		}
	}
	else
	{
		sample_counts = PyArray_SimpleNew(1,&category_len,NPY_LONG);
		if(!sample_counts) goto error;
		PyArray_FILLWBYTE(sample_counts, 0);
	}

	if(locus_counts)
	{
		if(!CountArray_Check2(locus_counts,len,PyArray_DIMS(locus_counts)[1]))
		{
			PyErr_SetString(PyExc_ValueError,"invalid locus count array");
			goto error;
		}
		state.max_genos = PyArray_DIMS(locus_counts)[1];
	}
	else
	{
		if(!len)
			state.max_genos = 0;
		if(GenotypeArray_Check(genos))
		{
			GenotypeArrayObject *garray = (GenotypeArrayObject *)genos;
			state.max_genos = 1<<garray->descriptor->max_bit_size;
		}
		else
		{
			Py_ssize_t i;
			int max_bit_size=0;

			for(i = 0; i<len; ++i)
			{
				GenotypeObject *geno = (GenotypeObject *)PySequence_GetItem(genos, i); /* new ref */

				if(!geno || !Genotype_CheckExact(geno))
				{
					PyErr_Format(GenotypeRepresentationError,
					    "invalid genotype object in genos at index %zd", i);
					Py_XDECREF(geno);
					goto error;
				}
				if(geno->model->bit_size > max_bit_size)
					max_bit_size = geno->model->bit_size;
				Py_DECREF(geno);
			}
			state.max_genos = 1<<max_bit_size;
		}
		locus_len[0] = len;
		locus_len[1] = state.max_genos;
		locus_counts = PyArray_SimpleNew(2,locus_len,NPY_LONG);
		if(!locus_counts) goto error;
		PyArray_FILLWBYTE(locus_counts, 0);
	}

	state.sample_counts = (unsigned long *)PyArray_DATA(sample_counts);
	state.locus_counts  = locus_counts;

	if(for_each_genotype(genos, (geno_foreach)sample_summary_foreach, &state) < 0)
		goto error;

	ret = PyTuple_Pack(2, sample_counts, locus_counts);

error:
	Py_XDECREF(sample_counts);
	Py_XDECREF(locus_counts);

	return ret;
}

static PyObject *
sample_summary_func(PyObject *self, PyObject *args, PyObject *kw)
{
	PyObject *genos, *sample_counts=NULL, *locus_counts=NULL;
	static char *kwlist[] = {"genos", "locus_counts", "sample_counts", 0};

	if(!PyArg_ParseTupleAndKeywords(args, kw, "O|OO:sample_summary", kwlist,
	                                &genos, &locus_counts, &sample_counts))
		return NULL;

	return sample_summary(genos, locus_counts, sample_counts);
}

/******************************************************************************************************/

static PyObject *
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

static PyObject *
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

static PyObject *
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

static PyObject *
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

		if(a->model != b->model)
		{
			PyErr_Format(GenotypeRepresentationError,
			    "genotype models do not match at index %zd", i);
			goto error;
		}

		/* If both genotypes are not missing */
		if(a->index && b->index)
		{
			if(a==b) concordant += 1;
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

static PyObject *
pick(PyObject *self, PyObject *args)
{
	Py_ssize_t i, j, n, len;
	PyObject *genos, *indices, *item;
	PyObject *seq=NULL, *indseq=NULL, *result=NULL;
	PyObject **seqitems, **indexitems;

	if(!PyArg_ParseTuple(args, "OO", &genos, &indices))
		return NULL;

	if( GenotypeArray_Check(genos) )
		return genoarray_pick( (GenotypeArrayObject *)genos, indices);

	indseq = PySequence_Fast(indices,"cannot convert indices into a sequence");
	if(!indseq) return NULL;

	seq = PySequence_Fast(genos,"cannot convert genotypes into a sequence");
	if(!seq) goto error;

	n          = PySequence_Fast_GET_SIZE(seq);
	seqitems   = PySequence_Fast_ITEMS(seq);
	len        = PySequence_Fast_GET_SIZE(indseq);
	indexitems = PySequence_Fast_ITEMS(indseq);

	result = PyList_New(len);
	if(!result) goto error;

	for(i=0; i<len; i++)
	{
		j = PyNumber_AsSsize_t(indexitems[i], PyExc_IndexError);

		if(j==-1 && PyErr_Occurred()) goto error;
		if(j<0 || j>=n)
		{
			PyErr_SetString(PyExc_IndexError,"genotype index out of range");
			goto error;
		}
		item = seqitems[j];

		Py_INCREF(item);
		PyList_SET_ITEM(result, i, item);
	}

	Py_DECREF(indseq);
	Py_DECREF(seq);
	return result;

error:
	Py_XDECREF(indseq);
	Py_XDECREF(seq);
	Py_XDECREF(result);
	return NULL;
}

static PyObject *
pick_columns_index(PyObject *rows, Py_ssize_t index)
{
	Py_ssize_t i, len;
	PyObject *item, *result=NULL;
	PyObject **items;

	rows = PySequence_Fast(rows,"cannot convert rows into a sequence");
	if(!rows) goto error;

	len   = PySequence_Fast_GET_SIZE(rows);
	items = PySequence_Fast_ITEMS(rows);

	result = PyList_New(len);
	if(!result) goto error;

	for(i=0; i<len; i++)
	{
		item = PySequence_GetItem(items[i], index); /* new ref */
		if(!item) goto error;
		PyList_SET_ITEM(result, i, item); /* steal ref */
	}

	Py_DECREF(rows);
	return result;

error:
	Py_XDECREF(rows);
	Py_XDECREF(result);
	return NULL;
}

static PyObject *
pick_columns_slice(PyObject *rows, PyObject *slice)
{
	Py_ssize_t i, j, k, rowlen, itemlen;
	Py_ssize_t start, stop, step, slicelen;
	PyObject *item, *result=NULL;
	PyObject **rowitems, **resultitems;

	rows = PySequence_Fast(rows,"cannot convert rows into a sequence");
	if(!rows) return NULL;

	rowlen     = PySequence_Fast_GET_SIZE(rows);
	rowitems   = PySequence_Fast_ITEMS(rows);

	if(!rowlen)
	{
		result = PyList_New(0);
		goto done;
	}

	itemlen = PySequence_Size(rowitems[0]);
	if(itemlen == -1) goto error;

	if(!slice || slice == Py_None)
	{
		start    = 0;
		stop     = itemlen;
		step     = 1;
		slicelen = itemlen;
	}
	else if(PySlice_Check(slice))
	{
		if(PySlice_GetIndicesEx( (PySliceObject *)slice, itemlen,
		                         &start, &stop, &step, &slicelen) < 0)
			goto error;
	}
	else
	{
		PyErr_SetString(PyExc_TypeError, "invalid index");
		goto error;
	}

	result = PyList_New(slicelen);
	if(!result) goto error;

	resultitems = PySequence_Fast_ITEMS(result);

	for(i=0; i<slicelen; ++i)
	{
		item = PyList_New(rowlen);
		if(!item) goto error;
		PyList_SET_ITEM(result, i, item);
	}

	for(j=0; j<rowlen; j++)
		for (k=start, i=0; i<slicelen; k+=step, i++)
		{
			item = PySequence_GetItem(rowitems[j], k); /* new ref */
			if(!item) goto error;
			PyList_SET_ITEM(resultitems[i], j, item); /* steal ref */
		}

done:
	Py_DECREF(rows);
	return result;

error:
	Py_XDECREF(rows);
	Py_XDECREF(result);
	return NULL;
}

static PyObject *
pick_columns_indices(PyObject *rows, PyObject *indices)
{
	Py_ssize_t i, j, index, rowlen, indexlen;
	Py_ssize_t *indexarray=NULL;
	PyObject *item, *result=NULL;
	PyObject **rowitems, **indexitems, **resultitems;

	rows = PySequence_Fast(rows,"cannot convert rows into a sequence");
	if(!rows) return NULL;

	indices = PySequence_Fast(indices,"cannot convert indices into a sequence");
	if(!indices) goto error;

	rowlen     = PySequence_Fast_GET_SIZE(rows);
	rowitems   = PySequence_Fast_ITEMS(rows);
	indexlen   = PySequence_Fast_GET_SIZE(indices);
	indexitems = PySequence_Fast_ITEMS(indices);

	result = PyList_New(indexlen);
	if(!result) goto error;

	resultitems = PySequence_Fast_ITEMS(result);

	indexarray = malloc( sizeof(Py_ssize_t)*indexlen );
	if(!indexarray) goto error;

	for(i=0; i<indexlen; ++i)
	{
		item = PyList_New(rowlen);
		if(!item) goto error;
		PyList_SET_ITEM(result, i, item);
		index = PyNumber_AsSsize_t(indexitems[i], PyExc_IndexError);
		if(index == -1 && PyErr_Occurred()) goto error;
		indexarray[i] = index;
	}

	for(j=0; j<rowlen; j++)
		for(i=0; i<indexlen; i++)
		{
			item = PySequence_GetItem(rowitems[j], indexarray[i]); /* new ref */
			if(!item) goto error;
			PyList_SET_ITEM(resultitems[i], j, item); /* steal ref */
		}

	free(indexarray);
	Py_DECREF(rows);
	Py_DECREF(indices);
	return result;

error:
	if(indexarray)
		free(indexarray);
	Py_XDECREF(rows);
	Py_XDECREF(indices);
	Py_XDECREF(result);
	return NULL;
}

static PyObject *
pick_columns(PyObject *self, PyObject *args)
{
	PyObject *rows, *index = Py_None;

	if(!PyArg_ParseTuple(args, "O|O", &rows, &index))
		return NULL;

	if(PyIndex_Check(index))
	{
		Py_ssize_t i = PyNumber_AsSsize_t(index, PyExc_IndexError);
		if(i == -1 && PyErr_Occurred())
			return NULL;
		return pick_columns_index(rows, i);
	}
	else if(index == Py_None || PySlice_Check(index))
		return pick_columns_slice(rows, index);
	else if(PySequence_Check(index))
		return pick_columns_indices(rows, index);
	else
	{
		PyErr_SetString(PyExc_TypeError, "indices must be None, integer, sequence, or slice");
		return NULL;
	}
}

/******************************************************************************************************/

static PyObject *
place(PyObject *self, PyObject *args)
{
	Py_ssize_t i, j, srclen, indexlen, destlen;
	PyObject *dest, *src, *indices;
	PyObject *srcseq=NULL, *indseq=NULL;
	PyObject **srcitems, **indexitems;

	if(!PyArg_ParseTuple(args, "OOO", &dest, &src, &indices))
		return NULL;

	destlen = PySequence_Size(dest);
	if(destlen == -1) goto error;

	srcseq = PySequence_Fast(src,"cannot convert source data into a sequence");
	if(!srcseq) goto error;

	indseq = PySequence_Fast(indices,"cannot convert indices into a sequence");
	if(!indseq) goto error;

	srclen   = PySequence_Fast_GET_SIZE(srcseq);
	indexlen = PySequence_Fast_GET_SIZE(indseq);

	if(srclen != indexlen)
	{
		PyErr_SetString(PyExc_ValueError,"source data and index length mismatch");
		goto error;
	}

	srcitems   = PySequence_Fast_ITEMS(srcseq);
	indexitems = PySequence_Fast_ITEMS(indseq);

	for(i=0; i<indexlen; i++)
	{
		j = PyNumber_AsSsize_t(indexitems[i], PyExc_IndexError);
		if(j==-1 && PyErr_Occurred()) goto error;

		if( PySequence_SetItem(dest, j, srcitems[i]) < 0 ) goto error;
	}

	Py_DECREF(srcseq);
	Py_DECREF(indseq);
	Py_INCREF(dest);
	return dest;

error:
	Py_XDECREF(srcseq);
	Py_XDECREF(indseq);
	return NULL;
}

/******************************************************************************************************/

static PyObject *
place_list(PyObject *self, PyObject *args)
{
	Py_ssize_t i, srclen, indexlen, destlen;
	PyObject *dest, *src, *indices, *srcitem, *destitem;
	PyObject *srcseq=NULL, *indseq=NULL;
	PyObject **srcitems, **indexitems;

	if(!PyArg_ParseTuple(args, "OOO", &dest, &src, &indices))
		return NULL;

	destlen = PySequence_Size(dest);
	if(destlen == -1) goto error;

	srcseq = PySequence_Fast(src,"cannot convert source data into a sequence");
	if(!srcseq) goto error;

	indseq = PySequence_Fast(indices,"cannot convert indices into a sequence");
	if(!indseq) goto error;

	srclen   = PySequence_Fast_GET_SIZE(srcseq);
	indexlen = PySequence_Fast_GET_SIZE(indseq);

	if(srclen != indexlen)
	{
		PyErr_SetString(PyExc_ValueError,"source data and index length mismatch");
		goto error;
	}

	srcitems   = PySequence_Fast_ITEMS(srcseq);
	indexitems = PySequence_Fast_ITEMS(indseq);

	for(i=0; i<indexlen; i++)
	{
		int ret;
		Py_ssize_t j = PyNumber_AsSsize_t(indexitems[i], PyExc_IndexError);

		if(j==-1 && PyErr_Occurred()) goto error;

		srcitem  = srcitems[i];
		destitem = PySequence_GetItem(dest, j); /* new ref */

		if(!destitem)
		{
			goto error;
		}
		else if(srcitem == Py_None)
		{
			/* pass */
			Py_DECREF(destitem);
		}
		else if(destitem == Py_None)
		{
			Py_DECREF(destitem);
			if(PyList_CheckExact(srcitem) && PyList_GET_SIZE(srcitem)==1)
				srcitem = PyList_GET_ITEM(srcitem,0); /* borrowed ref */
			if( PySequence_SetItem(dest, j, srcitem) < 0)
				goto error;
		}
		else if(PyList_CheckExact(destitem) && PyList_CheckExact(srcitem))
		{
			PyObject *items = PySequence_InPlaceConcat(destitem, srcitem);
			/* items should equal destitem */
			Py_DECREF(destitem);
			if(!items || PySequence_SetItem(dest, j, items) < 0)
			{
				Py_XDECREF(items);
				goto error;
			}
			Py_DECREF(items);
		}
		else if(PyList_CheckExact(destitem))
		{
			int ret = PyList_Append(destitem, srcitem);
			Py_DECREF(destitem);
			if(ret < 0) goto error;
		}
		else if(PyList_CheckExact(srcitem))
		{
			ret = PyList_Append(srcitem, destitem);
			Py_DECREF(destitem);
			if(ret < 0 || PySequence_SetItem(dest, j, srcitem) < 0)
				goto error;
		}
		else /* neither None or lists */
		{
			PyObject *items = PyList_New(2);
			if(!items)
			{
				Py_DECREF(destitem);
				goto error;
			}
			PyList_SET_ITEM(items, 0, destitem);
			PyList_SET_ITEM(items, 1, srcitem);
			Py_INCREF(srcitem);

			ret = PySequence_SetItem(dest, j, items);
			Py_DECREF(items);

			if(ret < 0) goto error;
		}
	}

	Py_DECREF(srcseq);
	Py_DECREF(indseq);
	Py_INCREF(dest);
	return dest;

error:
	Py_XDECREF(srcseq);
	Py_XDECREF(indseq);
	return NULL;
}

/******************************************************************************************************/

static GenotypeObject *
merge_unanimous(UnphasedMarkerModelObject *model, PyObject *genos, Py_ssize_t *status)
{
	Py_ssize_t i, len, found=0;
	PyObject *genoseq=NULL;
	GenotypeObject *genofound=NULL, *missing=NULL;
	PyObject **items;

	*status = MERGE_MISSING;

	if(!UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		goto error;
	}

	if(genos==Py_None)
	{
		genofound = (GenotypeObject *)PyList_GetItem(model->genotypes, 0); /* borrowed ref */
		Py_XINCREF(genofound);
		return genofound;
	}

	if(Genotype_CheckExact(genos))
	{
		GenotypeObject *genofound = (GenotypeObject *)genos;
		if(genofound->model != model)
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype model in genos at index %zd", 0L);
			goto error;
		}
		if(genofound->index != 0)
			*status = MERGE_UNAMBIGUOUS;
		Py_INCREF(genofound);
		return genofound;
	}

	/* Add check temporarily to find slow code. */
	if(!FastSequenceCheck(genos))
	{
		PyErr_SetString(PyExc_ValueError,"genos are not a tuple or list");
		goto error;
	}

	genoseq = PySequence_Fast(genos,"genos is not a genotype or a sequence");
	if(!genoseq) goto error;

	len   = PySequence_Fast_GET_SIZE(genoseq);
	items = PySequence_Fast_ITEMS(genoseq);

	for(i=0, found=0; i<len; ++i)
	{
		GenotypeObject *geno = (GenotypeObject *)items[i];

		if(!geno || !Genotype_CheckExact(geno))
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype object in genos at index %zd", i);
			goto error;
		}

		if(geno->model != model)
		{
			PyErr_Format(GenotypeRepresentationError,
			    "invalid genotype model in genos at index %zd", i);
			goto error;
		}

		if(geno->index == 0)
		{
			missing = geno;
			continue;
		}
		else if(!genofound)
			genofound = geno;
		else if(genofound != geno)
		{
			*status = MERGE_DISCORDANT;
			found = 0;
			break;
		}
		found++;
	}

	if(!found)
	{
		if(missing)
			genofound = missing;
		else
			genofound = (GenotypeObject *)PyList_GetItem(model->genotypes, 0); /* borrowed ref */
	}
	else if(found == 1)
		*status = MERGE_UNAMBIGUOUS;
	else /* found > 1 */
		*status = MERGE_CONCORDANT;

	Py_XDECREF(genoseq);
	Py_XINCREF(genofound);
	return genofound;

error:
	Py_XDECREF(genoseq);
	return NULL;
}

static PyObject *
merge_unanimous_func(PyObject *self, PyObject *args)
{
	Py_ssize_t status;
	PyObject *genos, *geno, *result=NULL;
	UnphasedMarkerModelObject *model;

	if(!PyArg_ParseTuple(args, "OO", &model, &genos))
		return NULL;

	geno = (PyObject *)merge_unanimous(model, genos, &status);
	if(!geno) return NULL;

	result = PyTuple_Pack(2, PyInt_FromSsize_t(status), geno);
	Py_DECREF(geno);

	return result;
}

/******************************************************************************************************/

static int
genomerger_init(GenotypeMergerObject *self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"mergefunc", "trackstats", NULL};

	self->cmergefunc  = NULL;
	self->pymergefunc = self->samplestats = self->locusstats = NULL;
	self->trackstats  = 1;

	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist,
		&self->pymergefunc, &self->trackstats))
		return -1;

	if(PyCFunction_Check(self->pymergefunc) &&
	   PyCFunction_GetFunction(self->pymergefunc) == merge_unanimous_func)
	{
		self->cmergefunc = merge_unanimous;
	}

	/* FIXME: Add error checking */
	self->locusstats  = PyDict_New();
	self->samplestats = PyDict_New();

	Py_INCREF(self->pymergefunc);
	return 0;
}

static int
genomerger_clear(GenotypeMergerObject *self)
{
	Py_CLEAR(self->pymergefunc);
	Py_CLEAR(self->samplestats);
	Py_CLEAR(self->locusstats);
	return 0;
}

static void
genomerger_dealloc(GenotypeMergerObject *self)
{
	genomerger_clear(self);
	self->ob_type->tp_free((PyObject *)self);
}

static int
genomerger_traverse(GenotypeMergerObject *self, visitproc visit, void *arg)
{
	Py_VISIT(self->pymergefunc);
	return 0;
}

static PyObject *
genomerger_merge_genotype_py(GenotypeMergerObject *self, PyObject *model, PyObject *genos, Py_ssize_t *status)
{
	PyObject *result;
	GenotypeObject *geno;

	result = PyObject_CallFunction(self->pymergefunc, "(OO)", model, genos); /* new ref */
	if(!result) return NULL;

	if(!PyTuple_Check(result) || PyTuple_GET_SIZE(result) != 2
	|| !Genotype_CheckExact(PyTuple_GET_ITEM(result,1)))
	{
		PyErr_SetString(PyExc_TypeError,"genotype merger must return a 2-tuple of status and genotype");
		Py_DECREF(result);
		return NULL;
	}
	*status = PyInt_AsSsize_t(PyTuple_GET_ITEM(result,0));
	if(*status==-1 && PyErr_Occurred())
	{
		Py_DECREF(result);
		return NULL;
	}
	geno = (GenotypeObject *)PyTuple_GET_ITEM(result,1); /* borrowed ref */
	Py_INCREF(geno);
	Py_DECREF(result);
	return (PyObject *)geno;
}

static inline PyObject *
genomerger_merge_genotype(GenotypeMergerObject *self, PyObject *model, PyObject *genos, Py_ssize_t *status)
{
	if(self->cmergefunc)
		return (PyObject *)self->cmergefunc( (UnphasedMarkerModelObject *)model, genos, status);
	else
		return genomerger_merge_genotype_py(self, model, genos, status);
}

static inline PyObject *
genomerger_get_stats(PyObject *sdict, PyObject *key)
{
	PyObject *stats;

	stats = PyDict_GetItem(sdict, key); /* borrwed ref */

	if(!stats && PyErr_Occurred()) return NULL;

	if(!stats)
	{
		npy_intp dims = 5;

		stats = PyArray_SimpleNew(1, &dims, NPY_LONG);
		if(!stats) return NULL;
		PyArray_FILLWBYTE(stats, 0);

		if(PyDict_SetItem(sdict, key, stats) < 0)
		{
			Py_DECREF(stats);
			return NULL;
		}
		Py_DECREF(stats);
	}
	/* borrowed ref */
	return stats;
}

static inline int
genomerger_update_stats(PyObject *stats, Py_ssize_t status)
{
	Py_ssize_t *statarray;

	if(!CountArray_Check(stats,5))
	{
		PyErr_SetString(PyExc_ValueError,"invalid statistics array");
		return -1;
	}

	statarray = PyArray_DATA(stats);

	if(status < 0 || status > 4)
	{
		PyErr_SetString(PyExc_IndexError,"invalid genotype merge status");
		return -1;
	}
	statarray[status] += 1;
	return 0;
}

static PyObject *
genomerger_merge_geno(GenotypeMergerObject *self, PyObject *args)
{
	PyObject *sample, *locus, *model, *genos, *geno;
	Py_ssize_t status;

	if(!PyArg_ParseTuple(args, "OOOO", &sample, &locus, &model, &genos))
		return NULL;

	geno = genomerger_merge_genotype(self, model, genos, &status);

	if(!geno) goto error;

	if(self->trackstats)
	{
		PyObject *stats;

		/* borrowed ref: stats */
		stats = genomerger_get_stats(self->samplestats, sample);
		if(!stats) goto error;
		if( genomerger_update_stats(stats, status)<0 ) goto error;

		stats = genomerger_get_stats(self->locusstats, locus);
		if(!stats) goto error;
		if( genomerger_update_stats(stats, status)<0 ) goto error;
	}

	return geno;

error:
	Py_XDECREF(geno);
	return NULL;
}

static PyObject *
genomerger_merge_sample(GenotypeMergerObject *self, PyObject *args)
{
	PyObject *sample, *loci, *models, *genos, *geno, *samplestats=NULL;
	PyObject *locusseq=NULL, *modelseq=NULL, *genoseq=NULL, *result=NULL;
	PyObject **locusitems, **modelitems, **genoitems;
	Py_ssize_t locuslen, modellen, genolen, status, i;

	if(!PyArg_ParseTuple(args, "OOOO", &sample, &loci, &models, &genos))
		return NULL;

	/* Add checks temporarily to find slow code. */
	if(!FastSequenceCheck(loci))
	{
		PyErr_SetString(PyExc_ValueError,"loci are not a tuple or list");
		goto error;
	}

	if(!FastSequenceCheck(models))
	{
		PyErr_SetString(PyExc_ValueError,"models are not a tuple or list");
		goto error;
	}

	if(!FastSequenceCheck(genos))
	{
		PyErr_SetString(PyExc_ValueError,"genos are not a tuple or list");
		goto error;
	}

	locusseq = PySequence_Fast(loci,"cannot convert loci into a sequence");
	if(!locusseq) goto error;

	modelseq = PySequence_Fast(models,"cannot convert models into a sequence");
	if(!modelseq) goto error;

	genoseq = PySequence_Fast(genos,"cannot convert genotypes into a sequence");
	if(!genoseq) goto error;

	locuslen = PySequence_Fast_GET_SIZE(locusseq);
	modellen = PySequence_Fast_GET_SIZE(modelseq);
	genolen  = PySequence_Fast_GET_SIZE(genoseq);

	if(locuslen != modellen || modellen != genolen)
	{
		PyErr_SetString(PyExc_ValueError,"length mismatch among loci, models, and genotypes");
		goto error;
	}

	locusitems = PySequence_Fast_ITEMS(locusseq);
	modelitems = PySequence_Fast_ITEMS(modelseq);
	genoitems  = PySequence_Fast_ITEMS(genoseq);

	result = PyList_New(genolen);
	if(!result) goto error;

	/* borrowed ref: samplestats */
	if(self->trackstats)
	{
		samplestats = genomerger_get_stats(self->samplestats, sample);
		if(!samplestats) goto error;
	}

	for(i=0; i<genolen; ++i)
	{
		geno = genomerger_merge_genotype(self, modelitems[i], genoitems[i], &status);
		if(!geno) goto error;
		PyList_SET_ITEM(result, i, geno);

		if(self->trackstats)
		{
			PyObject *locusstats;

			if( genomerger_update_stats(samplestats, status)<0 ) goto error;

			/* borrowed ref: locusstats */
			locusstats = genomerger_get_stats(self->locusstats, locusitems[i]);
			if(!locusstats) goto error;
			if( genomerger_update_stats(locusstats, status)<0 ) goto error;
		}
	}

	Py_DECREF(locusseq);
	Py_DECREF(modelseq);
	Py_DECREF(genoseq);

	return result;

error:
	Py_XDECREF(locusseq);
	Py_XDECREF(modelseq);
	Py_XDECREF(genoseq);
	Py_XDECREF(result);

	return NULL;
}

static PyObject *
genomerger_merge_locus(GenotypeMergerObject *self, PyObject *args)
{
	PyObject *samples, *locus, *model, *genos, *geno, *locusstats=NULL;
	PyObject *sampleseq=NULL, *genoseq=NULL, *result=NULL;
	PyObject **sampleitems, **genoitems;
	Py_ssize_t samplelen, genolen, status, i;

	if(!PyArg_ParseTuple(args, "OOOO", &samples, &locus, &model, &genos))
		return NULL;

	/* Add checks temporarily to find slow code. */
	if(!FastSequenceCheck(samples))
	{
		PyErr_SetString(PyExc_ValueError,"samples are not a tuple or list");
		goto error;
	}

	if(!FastSequenceCheck(genos))
	{
		PyErr_SetString(PyExc_ValueError,"genos are not a tuple or list");
		goto error;
	}

	sampleseq = PySequence_Fast(samples,"cannot convert samples into a sequence");
	if(!sampleseq) goto error;

	genoseq = PySequence_Fast(genos,"cannot convert genotypes into a sequence");
	if(!genoseq) goto error;

	samplelen = PySequence_Fast_GET_SIZE(sampleseq);
	genolen   = PySequence_Fast_GET_SIZE(genoseq);

	if(samplelen != genolen)
	{
		PyErr_SetString(PyExc_ValueError,"length mismatch between samples and genotypes");
		goto error;
	}

	sampleitems = PySequence_Fast_ITEMS(sampleseq);
	genoitems   = PySequence_Fast_ITEMS(genoseq);

	result = PyList_New(genolen);
	if(!result) goto error;

	if(self->trackstats)
	{
		/* borrowed ref: locusstats */
		locusstats = genomerger_get_stats(self->locusstats, locus);
		if(!locusstats) goto error;
	}

	for(i=0; i<genolen; ++i)
	{
		geno = genomerger_merge_genotype(self, model, genoitems[i], &status);
		if(!geno) goto error;
		PyList_SET_ITEM(result, i, geno);

		if(self->trackstats)
		{
			PyObject *samplestats;

			if( genomerger_update_stats(locusstats, status)<0 ) goto error;

			/* borrowed ref: samplestats */
			samplestats = genomerger_get_stats(self->samplestats, sampleitems[i]);
			if(!samplestats) goto error;
			if( genomerger_update_stats(samplestats, status)<0 ) goto error;
		}
	}

	Py_DECREF(sampleseq);
	Py_DECREF(genoseq);

	return result;

error:
	Py_XDECREF(sampleseq);
	Py_XDECREF(genoseq);
	Py_XDECREF(result);

	return NULL;
}

static PyMemberDef genomerger_members[] = {
	{"mergefunc",   T_OBJECT_EX, offsetof(GenotypeMergerObject, pymergefunc),   RO, "Merge function"},
	{"locusstats",	T_OBJECT_EX, offsetof(GenotypeMergerObject, locusstats),    RO, "Locus statistics"},
	{"samplestats",	T_OBJECT_EX, offsetof(GenotypeMergerObject, samplestats),   RO, "Sample statistics"},
	{"trackstats",  T_UINT,      offsetof(GenotypeMergerObject, trackstats),    RO, "Track statistics"},
	{NULL}  /* Sentinel */
};

static PyMethodDef genomerger_methods[] = {
	{"merge_geno",      (PyCFunction)genomerger_merge_geno, METH_VARARGS,
		"Merge a single genotype"},
	{"merge_sample", (PyCFunction)genomerger_merge_sample,  METH_VARARGS,
		"Merge genotypes for a given sample"},
	{"merge_locus",   (PyCFunction)genomerger_merge_locus,  METH_VARARGS,
		"Merge genotypes for a given locus"},
{NULL}  /* Sentinel */
};

PyTypeObject GenotypeMergerType = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size           */
	"GenotypeMerger",			/* tp_name           */
	sizeof(GenotypeMergerObject),		/* tp_basicsize      */
	0,					/* tp_itemsize       */
	(destructor)genomerger_dealloc,		/* tp_dealloc        */
	0,					/* tp_print          */
	0,					/* tp_getattr        */
	0,					/* tp_setattr        */
	0,					/* tp_compare        */
	0,					/* tp_repr           */
	0,					/* tp_as_number      */
	0,					/* tp_as_sequence    */
	0,					/* tp_as_mapping     */
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT,			/* tp_flags          */
	"Genotype Merger",			/* tp_doc            */
	(traverseproc)genomerger_traverse,	/* tp_traverse       */
	(inquiry)genomerger_clear,		/* tp_clear          */
	0,					/* tp_richcompare    */
	0,					/* tp_weaklistoffset */
	0,					/* tp_iter           */
	0,					/* tp_iternext       */
	genomerger_methods,			/* tp_methods        */
	genomerger_members,			/* tp_members        */
	0,					/* tp_getset         */
	0,					/* tp_base           */
	0,					/* tp_dict           */
	0,					/* tp_descr_get      */
	0,					/* tp_descr_set      */
	0,					/* tp_dictoffset     */
	(initproc)genomerger_init,		/* tp_init           */
	0,					/* tp_alloc          */
	0,					/* tp_new            */
	0,					/* tp_free           */
};

/******************************************************************************************************/

PyMODINIT_FUNC
init_genoarray(void)
{
	PyObject *m;

	static PyMethodDef genoarraymodule_methods[] = {
		{"genotype_indices",            (PyCFunction)genotype_indices_func,	METH_KEYWORDS,
		 "Return an array of integer genotype indices"},
		{"count_genotypes",	(PyCFunction)count_genotypes_func,	METH_KEYWORDS,
		 "Return an array of genotype counts by index"},
		{"genotype_categories",	(PyCFunction)genotypes_categories_func,	METH_KEYWORDS,
		 "Count the number of occurances of each genotypes category"},
		{"locus_summary",	(PyCFunction)locus_summary_func,	METH_KEYWORDS,
		 "Count genotypes and categories for genotypes at a single locus"},
		{"sample_summary",	(PyCFunction)sample_summary_func,	METH_KEYWORDS,
		 "Count genotypes and categories for genotypes for a single sample"},
		{"genoarray_concordance",	genoarray_concordance,	METH_VARARGS,
		 "Generate simple concordance statistics from two genotype arrays"},
		{"genoarray_concordance_8bit",	genoarray_concordance_8bit,	METH_VARARGS,
		 "Generate simple concordance statistics from two genotype arrays stored in 8-bit format"},
		{"genoarray_concordance_4bit",genoarray_concordance_4bit,	METH_VARARGS,
		 "Generate simple concordance statistics from two genotype arrays stored in 4-bit format"},
		{"genoarray_concordance_2bit",	genoarray_concordance_2bit,	METH_VARARGS,
		 "Generate simple concordance statistics from two genotype arrays stored in 2-bit format"},
		{"pick",	pick,	METH_VARARGS,
		 "Pick items from a sequence given indices"},
		{"pick_columns",	pick_columns,	METH_VARARGS,
		 "Pick one or more columns from a two dimensional sequence"},
		{"place",	place,	METH_VARARGS,
		 "Place items from a sequence at indices into a destination sequence"},
		{"place_list",	place_list,	METH_VARARGS,
		 "Place items from a sequence at indices into a destination sequence concatenating into lists if items are not None"},
		{"merge_unanimous",	merge_unanimous_func,	METH_VARARGS,
		 "Merge a list of genotypes requiring non-missing genotypes to be unanimous or they are setting to missing"},
		{"count_haplotypes",	count_haplotypes,	METH_VARARGS, "Count haplotypes at two loci"},
		{"estimate_ld",	estimate_ld,	METH_VARARGS, "Compute LD statistics from haplotype counts"},
		{NULL}  /* Sentinel */
	};

	import_array();

	GenotypeType.tp_new                = PyType_GenericNew;
	UnphasedMarkerModelType.tp_new     = PyType_GenericNew;
	GenotypeArrayDescriptorType.tp_new = PyType_GenericNew;
	GenotypeMergerType.tp_new          = PyType_GenericNew;

	if(PyType_Ready(&GenotypeType)                < 0 ||
	   PyType_Ready(&GenotypeArrayType)           < 0 ||
	   PyType_Ready(&UnphasedMarkerModelType)     < 0 ||
	   PyType_Ready(&GenotypeArrayDescriptorType) < 0 ||
	   PyType_Ready(&GenotypeMergerType)          < 0)
		return;

	m = Py_InitModule3("_genoarray", genoarraymodule_methods,
		"Genotype array module");

	GenotypeLookupError = PyErr_NewException("_genoarray.GenotypeLookupError",
		PyExc_KeyError, NULL);
	if(GenotypeLookupError == NULL)
		return;

	GenotypeRepresentationError = PyErr_NewException("_genoarray.GenotypeRepresentationError",
		PyExc_ValueError, NULL);
	if(GenotypeRepresentationError == NULL)
		return;

	/* Increment reference counts on all static types to ensure they are not freed */
	Py_INCREF(&GenotypeArrayDescriptorType);
	Py_INCREF(&UnphasedMarkerModelType);
	Py_INCREF(&GenotypeArrayType);
	Py_INCREF(&GenotypeType);

	PyModule_AddObject(m, "GenotypeLookupError",        GenotypeLookupError);
	PyModule_AddObject(m, "GenotypeRepresentationError",GenotypeRepresentationError);
	PyModule_AddObject(m, "GenotypeArrayDescriptor",    (PyObject *)&GenotypeArrayDescriptorType);
	PyModule_AddObject(m, "UnphasedMarkerModel",        (PyObject *)&UnphasedMarkerModelType);
	PyModule_AddObject(m, "GenotypeArray",              (PyObject *)&GenotypeArrayType);
	PyModule_AddObject(m, "Genotype",                   (PyObject *)&GenotypeType);
	PyModule_AddObject(m, "GenotypeMerger",             (PyObject *)&GenotypeMergerType);
}
