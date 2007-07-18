#include "Python.h"
#include "structmember.h"
#include "numpy/arrayobject.h"

#ifdef STDC_HEADERS
#include <stddef.h>
#else
#include <sys/types.h>		/* For size_t */
#endif

#include <math.h>

/* Declarations from bitarray */
unsigned long
bitarray_getbits(const unsigned char *data, ssize_t buflen, ssize_t i, ssize_t readlen, char **status);

void
bitarray_setbits(unsigned char *data, ssize_t buflen, ssize_t i, unsigned long value, ssize_t writelen, char **status);

/* Genotype Array objects declarations */

typedef struct {
	PyObject_HEAD
	PyObject      *models;
	PyArrayObject *offsets;
	unsigned long  bit_size;
	unsigned long  byte_size;
} GenotypeArrayDescriptorObject;

typedef struct {
	PyObject_HEAD
	GenotypeArrayDescriptorObject *descriptor;
	unsigned char data[];
} GenotypeArrayObject;

typedef struct {
	PyObject_HEAD
	PyObject      *genomap;
	PyObject      *genotypes;
	PyObject      *genotuples;
	PyObject      *alleles;
	unsigned short max_alleles;
	unsigned short allow_hemizygote;
	unsigned short bit_width;
} UnphasedMarkerModelObject;

typedef struct {
	PyObject_HEAD
	UnphasedMarkerModelObject *model;
	unsigned long  index;
	unsigned int   allele1;
	unsigned int   allele2;
} GenotypeObject;

/* Forward declaration */
PyTypeObject UnphasedMarkerModelType;
PyTypeObject GenotypeArrayType;

#define GenotypeArray_Check(op)                PyObject_TypeCheck(op, &GenotypeArrayType)
#define GenotypeArray_CheckExact(op)           ((op)->ob_type == &GenotypeArrayType)
#define UnphasedMarkerModel_Check(op)         (((op)->ob_type == &UnphasedMarkerModelType) || PyObject_TypeCheck(op, &UnphasedMarkerModelType))
#define UnphasedMarkerModel_CheckExact(op)     ((op)->ob_type == &UnphasedMarkerModelType)
#define Genotype_CheckExact(op)                ((op)->ob_type == &GenotypeType)
#define GenotypeArrayDescriptor_CheckExact(op) ((op)->ob_type == &GenotypeArrayDescriptorType)

static int
genotype_init(GenotypeObject *self, PyObject *args, PyObject *kwds)
{
	UnphasedMarkerModelObject *model;
	Py_ssize_t index, allele1, allele2, allele_len, genotype_len;

	static char *kwlist[] = {"model", "allele1", "allele2", "index", NULL};

	self->allele1 = self->allele2 = self->index = 0;
	Py_CLEAR(self->model);

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "Onnn", kwlist,
		&model, &allele1, &allele2, &index))
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

	if(allele1 < 0 || allele1 >= allele_len ||
	   allele2 < 0 || allele2 >= allele_len)
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
	self->allele1 = allele1;
	self->allele2 = allele2;
	self->index   = index;

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
	self->ob_type->tp_free((PyObject*)self);
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
	PyObject *allele1;

	allele1 = PyList_GetItem(self->model->alleles, self->allele1);
	if(allele1)
		Py_INCREF(allele1);
	return allele1;

}

static PyObject *
genotype_allele2_get(GenotypeObject *self)
{
	PyObject *allele2;

	allele2 = PyList_GetItem(self->model->alleles, self->allele2);
	if(allele2)
		Py_INCREF(allele2);
	return allele2;
}

static PyMemberDef genotype_members[] = {
	{"model",   T_OBJECT_EX, offsetof(GenotypeObject, model),   RO, "UnphasedMarkerModel"},
	{"index",   T_ULONG,     offsetof(GenotypeObject, index),   RO, "index"},
	{NULL}  /* Sentinel */
};

static PyGetSetDef genotype_getsetlist[] = {
	{"allele1", (getter)genotype_allele1_get, NULL, NULL},
	{"allele2", (getter)genotype_allele2_get, NULL, NULL},
	{NULL, NULL, NULL, NULL},  /* Sentinel */
};

static PyObject *
genotype_alleles(GenotypeObject *self)
{
	PyObject *result = PyList_GetItem(self->model->genotuples, self->index);
	Py_INCREF(result);
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
	if( self->allele1 && self->allele2 && self->allele1 != self->allele2 )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_homozygote(GenotypeObject *self)
{
	if( self->allele1 && self->allele1 == self->allele2 )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_hemizygote(GenotypeObject *self)
{
	if( !!self->allele1 ^ !!self->allele2 )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static PyObject *
genotype_missing(GenotypeObject *self)
{
	if( !self->allele1 && !self->allele2 )
		Py_RETURN_TRUE;
	else
		Py_RETURN_FALSE;
}

static Py_ssize_t
genotype_length(GenotypeObject *self)
{
	Py_ssize_t n = 0;
	if(self->allele1) n++;
	if(self->allele2) n++;
	return n;
}

static PyObject *
genotype_item(GenotypeObject *self, Py_ssize_t item)
{
	PyObject *allele;
	if(item == 0)
		allele = PyList_GetItem(self->model->alleles, self->allele1);
	else if(item == 1)
		allele = PyList_GetItem(self->model->alleles, self->allele2);
	else
	{
		PyErr_SetString(PyExc_IndexError,"genotype index out of range");
		return NULL;
	}
	if(allele)
		Py_INCREF(allele);
	return allele;
}

static int
genotype_contains(GenotypeObject *self, PyObject *allele)
{
	int cmp;
	PyObject *allele1, *allele2;

	allele1 = PyList_GetItem(self->model->alleles, self->allele1);
	cmp = PyObject_RichCompareBool(allele, allele1, Py_EQ);
	if(cmp == 0)
	{
		allele2 = PyList_GetItem(self->model->alleles, self->allele2);
		cmp = PyObject_RichCompareBool(allele, allele2, Py_EQ);
	}
	return cmp;
}

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

static PyTypeObject GenotypeType = {
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
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,/* tp_flags          */
	"Genotype objects",			/* tp_doc            */
	(traverseproc)genotype_traverse,	/* tp_traverse       */
	(inquiry)genotype_clear,		/* tp_clear          */
	0,					/* tp_richcompare    */
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
	PyType_GenericAlloc,			/* tp_alloc          */
	PyType_GenericNew,			/* tp_new            */
	PyObject_GC_Del,			/* tp_free           */
};

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
	unsigned long initial_offset = 0;
	unsigned long *offset_data;
	int dims;
	static char *kwlist[] = {"models", "initial_offset", NULL};

	descr_clear(self);
	self->byte_size = self->bit_size = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|n:GenotypeArrayDescriptor", kwlist,
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
	offsets = (PyArrayObject*)PyArray_FromDims(1, &dims, PyArray_ULONG);
	if(!offsets) return -1;

	offset_data = (unsigned long *)PyArray_DATA(offsets);
	offset_data[0] = initial_offset;

	for(i=0; i<n; ++i)
	{
		model = (UnphasedMarkerModelObject *)PyList_GetItem(models, i);
		if(!model) goto error;
		if(!UnphasedMarkerModel_Check(model) || model->bit_width > 32)
		{
			PyErr_SetString(PyExc_TypeError,"invalid genotype model");
			goto error;
		}
		offset_data[i+1] = offset_data[i] + model->bit_width;
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
	self->ob_type->tp_free((PyObject*)self);
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

	if (!self->models)
	{
		PyErr_SetString(PyExc_ValueError,"invalid descriptor state");
		return -1;
	}

	n = PyList_Size(self->models);

	if (n == -1 && PyErr_Occurred())
		return -1;

	return n;
}

static PyMemberDef descr_members[] = {
	{"models",    T_OBJECT_EX, offsetof(GenotypeArrayDescriptorObject, models),    RO, "UnphasedMarkerModels"},
	{"offsets",   T_OBJECT_EX, offsetof(GenotypeArrayDescriptorObject, offsets),   RO, "offsets"},
	{"bit_size",  T_ULONG,     offsetof(GenotypeArrayDescriptorObject, bit_size),  RO, "bit_size"},
	{"byte_size", T_ULONG,     offsetof(GenotypeArrayDescriptorObject, byte_size), RO, "byte_size"},
	{NULL}  /* Sentinel */
};

static PyMethodDef descr_methods[] = {
	{NULL}  /* Sentinel */
};

static PySequenceMethods descr_as_sequence = {
	(lenfunc)descr_length,			/* sq_length */
	0,					/* sq_concat */
	0,					/* sq_repeat */
	0,					/* sq_item */
	0,					/* sq_slice */
	0,					/* sq_ass_item */
	0,					/* sq_ass_slice */
	0,					/* sq_contains */
	0,					/* sq_inplace_concat */
	0,					/* sq_inplace_repeat */
};

static PyTypeObject GenotypeArrayDescriptorType = {
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
	0,					/* tp_as_mapping     */
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	0,					/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,/* tp_flags          */
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
	PyType_GenericAlloc,			/* tp_alloc          */
	PyType_GenericNew,			/* tp_new            */
	PyObject_GC_Del,			/* tp_free           */
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
	PyObject *index = PyObject_CallMethod(self->alleles, "index", "(O)", allele);

	if (PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_ValueError))
	{
		PyErr_Clear();
		PyErr_SetObject(PyExc_KeyError, allele);
	}
	return index;
}

static Py_ssize_t
genomodel_add_allele_internal(UnphasedMarkerModelObject *self, PyObject *allele)
{
	Py_ssize_t result;
	PyObject *index = PyObject_CallMethod(self->alleles, "index", "(O)", allele);

	if (PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_ValueError))
	{
		Py_ssize_t n;

		PyErr_Clear();

		/* Note: alleles[0] is always the missing allele */
		n = PyList_Size(self->alleles);
		if(n==-1) return -1;

		if(n > self->max_alleles)
		{
			PyErr_SetString(PyExc_ValueError,"genotype model cannot accomodate additional alleles");
			return -1;
		}

		if(PyList_Append(self->alleles, allele) == -1)
			return -1;

		return n;
	}
	result = PyInt_AsSsize_t(index);
	Py_DECREF(index);
	return result;
}

static PyObject *
genomodel_add_allele(UnphasedMarkerModelObject *self, PyObject *allele)
{
	Py_ssize_t result = genomodel_add_allele_internal(self, allele);

	if(result == -1)
		return NULL;
	else
		return PyInt_FromSsize_t(result);
}

static PyObject *
genomodel_get_genotype(UnphasedMarkerModelObject *self,  PyObject *geno)
{
	PyObject *g;

	if(!PyTuple_CheckExact(geno) || PyTuple_GET_SIZE(geno) != 2)
	{
		PyErr_SetString(PyExc_ValueError,"genotype must be specified as a 2-tuple");
		return NULL;
	}

	g = PyDict_GetItem(self->genomap, geno);
	if(!g)
	{
		PyErr_SetObject(PyExc_KeyError, geno);
		return NULL;
	}

	Py_INCREF(g);
	return g;
}

static PyObject *
genomodel_add_genotype(UnphasedMarkerModelObject *self, PyObject *geno)
{
	PyObject *allele1, *allele2, *args, *g;
	Py_ssize_t index1, index2, n;
	int res;

	args = NULL;

	if(!PyTuple_CheckExact(geno) || PyTuple_GET_SIZE(geno) != 2)
	{
		PyErr_SetString(PyExc_ValueError,"genotype must be specified as a 2-tuple");
		return NULL;
	}

	/* Fast path: Genotype already exists */
	g = PyDict_GetItem(self->genomap, geno);
	if(g)
	{
		Py_INCREF(g);
		return g;
	}

	/* Slow path: Add alleles and index the genotype */
	allele1 = PyTuple_GET_ITEM(geno, 0);
	allele2 = PyTuple_GET_ITEM(geno, 1);

	if( !self->allow_hemizygote && ((allele1 == Py_None) ^ (allele2 == Py_None)) )
	{
		PyErr_SetString(PyExc_ValueError,"model does not allow hemizgote genotypes");
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

	/* Replace alleles with internal representations to avoid immortalizing them */
	allele1 = PyList_GET_ITEM(self->alleles, index1);
	if(!allele1) return NULL;
	allele2 = PyList_GET_ITEM(self->alleles, index2);
	if(!allele2) return NULL;

	/* Create new genotype */
	n = PyList_Size(self->genotypes);
	if(n==-1) return NULL;

	g = PyObject_CallFunction( (PyObject*)&GenotypeType, "(Onnn)", self, index1, index2, n);
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

static int
genomodel_init(UnphasedMarkerModelObject *self, PyObject *args, PyObject *kw)
{
	PyObject *missing, *geno;
	PyObject *allow_hemizygote = Py_False;
	unsigned short n = 2;
	static char *kwlist[] = {"allow_hemizygote","max_alleles", 0};

	genomodel_clear(self);

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|OH:UnphasedMarkerModel", kwlist,
		&allow_hemizygote, &n))
		return -1;

	if(n == 0)
	{
		PyErr_SetString(PyExc_ValueError,"marker model requires max_alleles greater than 0");
		return -1;
	}

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
		self->bit_width = ceil(log2((n+1)*(n+2)/2));
	else
		self->bit_width = ceil(log2(n*(n+1)/2 + 1));

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
	{"bit_width",        T_USHORT,    offsetof(UnphasedMarkerModelObject, bit_width),        RO, "bit width"},
	{NULL}  /* Sentinel */
};

PyDoc_STRVAR(genomodel_doc,
"UnphasedMarkerModelObject(allow_hemizygote=False, max_alleles=2)\n");

PyTypeObject UnphasedMarkerModelType = {
	PyObject_HEAD_INIT(&PyType_Type)
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
	0,					/* tp_as_mapping     */
	0,					/* tp_hash           */
	0,					/* tp_call           */
	0,					/* tp_str            */
	PyObject_GenericGetAttr,		/* tp_getattro       */
	0,					/* tp_setattro       */
	0,					/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC|
		Py_TPFLAGS_BASETYPE,		/* tp_flags          */
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
	PyType_GenericAlloc,			/* tp_alloc          */
	PyType_GenericNew,			/* tp_new            */
	PyObject_GC_Del,			/* tp_free           */
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

	if (!self->descriptor || !GenotypeArrayDescriptor_CheckExact(self->descriptor))
	{
		PyErr_SetString(PyExc_TypeError,"invalid descriptor");
		return -1;
	}

	models = self->descriptor->models;
	if(!models) return -1;

	n = PyList_Size(models);
	if (n == -1 && PyErr_Occurred())
		return -1;

	return n;
}

static PyObject *
genoarray_repr(GenotypeArrayObject *self)
{
	PyObject *junk, *result;

	junk = PySequence_GetSlice( (PyObject*)self, 0, genoarray_length(self));
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

	if (PyType_IS_GC(type))
		obj = _PyObject_GC_Malloc(size);
	else
		obj = (PyObject *)PyObject_MALLOC(size);

	if (obj == NULL)
		return PyErr_NoMemory();

	memset(obj, '\0', size);

	if (type->tp_flags & Py_TPFLAGS_HEAPTYPE)
		Py_INCREF(type);

	PyObject_INIT(obj, type);

	if (PyType_IS_GC(type))
		_PyObject_GC_TRACK(obj);
	return obj;
}

static PyObject *
genoarray_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
	unsigned long n;
	GenotypeArrayObject *self;
	PyObject *descr;
	GenotypeArrayDescriptorObject *descriptor = NULL;;
	PyObject *genos = NULL;
	static char *kwlist[] = {"descriptor","genos", 0};

	if (!PyArg_ParseTupleAndKeywords(args, kw, "O|O:GenotypeArray", kwlist, &descr, &genos))
		return NULL;

	if(GenotypeArrayDescriptor_CheckExact(descr))
		descriptor = (GenotypeArrayDescriptorObject*)descr;
	else if(GenotypeArray_Check(descr))
		descriptor = ((GenotypeArrayObject *)descr)->descriptor;
	else
	{
		PyErr_SetString(PyExc_TypeError,"invalid descriptor object specified");
		return NULL;
	}

	n = descriptor->byte_size;
	self = (GenotypeArrayObject *)type->tp_alloc(&GenotypeArrayType, n);
	if(!self) return NULL;

	self->descriptor = descriptor;
	Py_INCREF(self->descriptor);

	if( genos && genos!=Py_None )
	{
		if( PySequence_SetSlice( (PyObject*)self, 0, genoarray_length(self), genos) == -1)
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

PyDoc_STRVAR(getitem_doc,
"x.__getitem__(y) <==> x[y]");

static PyObject *genoarray_subscript(GenotypeArrayObject*, PyObject*);

static PyMethodDef genoarray_methods[] = {
	{"__getitem__", (PyCFunction)genoarray_subscript, METH_O|METH_COEXIST, getitem_doc},
	{NULL,		NULL}		/* sentinel */
};

PyDoc_STRVAR(genoarray_doc,
"genoarray(descriptor) -> new empty genotype array\n"
"genoarray(descriptor, genos) -> new genotype array initialized with genos");

static PyObject *
genoarray_inner_get(PyObject *models, const unsigned char *data, Py_ssize_t datasize,
                    unsigned long *offsets, Py_ssize_t i)
{
	Py_ssize_t k, offset1, offset2;
	PyObject *geno;
	UnphasedMarkerModelObject *model;
	char *status;

	model = (UnphasedMarkerModelObject*)PyList_GET_ITEM(models, i);

	if(!model || !UnphasedMarkerModel_Check(model) || !model->genotypes)
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

	geno = PyList_GET_ITEM(model->genotypes, k);
	Py_INCREF(geno);
	return geno;
}

static int
genoarray_checkstate(GenotypeArrayObject *self)
{
	if (!self->data || !self->descriptor
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
	PyObject *geno, *models;
	unsigned long *offsets;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	geno     = NULL;
	datasize = self->descriptor->byte_size;

	offsets  = (unsigned long *)PyArray_DATA(self->descriptor->offsets);
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

	if(item < 0 || item >= n)
	{
		PyErr_SetString(PyExc_IndexError,"genotype index out of range");
		goto done;
	}

	geno  = genoarray_inner_get(models, self->data, datasize, offsets, item);
	if(!geno) goto done;
done:
	return geno;
}

static PyObject *
genoarray_slice(GenotypeArrayObject *self, PySliceObject *slice)
{
	Py_ssize_t i, j, n, datasize;
	Py_ssize_t start, stop, step, slicelength;
	PyObject *result, *models, *geno;
	unsigned long *offsets;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	result   = NULL;
	datasize = self->descriptor->byte_size;

	offsets  = (unsigned long *)PyArray_DATA(self->descriptor->offsets);
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

	if (PySlice_GetIndicesEx(slice, n, &start, &stop, &step, &slicelength) < 0)
		goto done;

	if (slicelength <= 0)
	{
		result = PyList_New(0);
		goto done;
	}

	result = PyList_New(slicelength);
	if (!result) goto done;

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
genoarray_subscript(GenotypeArrayObject *self, PyObject *item)
{
	if (PyIndex_Check(item)) {
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if (i == -1 && PyErr_Occurred())
			return NULL;
		return genoarray_item(self, i);
	}
	else if (PySlice_Check(item))
		return genoarray_slice(self, (PySliceObject *)item);
	else {
		PyErr_SetString(PyExc_TypeError, "indices must be integers");
		return NULL;
	}
}

static int
genoarray_inner_set(PyObject *models, PyObject *geno, unsigned char *data, Py_ssize_t datasize, unsigned long *offsets, Py_ssize_t i)
{
	UnphasedMarkerModelObject *model;
	GenotypeObject *g;
	Py_ssize_t offset1, offset2;
	char *status;

	model = (UnphasedMarkerModelObject *)PyList_GET_ITEM(models, i);

	if(!model || !UnphasedMarkerModel_Check(model))
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype model");
		return -1;
	}

	if(Genotype_CheckExact(geno))
	{
		g = (GenotypeObject *)geno;
		if(g->model != model)
		{
			/* Handle foreign genotypes by looking them up by their alleles */
			PyObject *genotuple = genotype_alleles(g);
			if(!genotuple) return -1;
			g = (GenotypeObject *)genomodel_get_genotype( (UnphasedMarkerModelObject *)model, genotuple);
			Py_DECREF(genotuple);
			if(!g) return -1;
			Py_DECREF(g); /* safe to treat as a borrowed reference */
		}
	}
	else if(PyTuple_CheckExact(geno) && PyTuple_GET_SIZE(geno) == 2)
	{
		g = (GenotypeObject *)genomodel_get_genotype( (UnphasedMarkerModelObject *)model, geno);
		if(!g) return -1;
		Py_DECREF(g); /* safe to treat as a borrowed reference */
	}
	else
	{
		PyErr_SetString(PyExc_TypeError,"invalid genotype object");
		return -1;
	}

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
	unsigned long *offsets;

	if (!value)
	{
		/* delete item */
		PyErr_SetString(PyExc_TypeError, "genoarray objects do not support item deletion");
		return -1;
	}

	if( genoarray_checkstate(self) == -1 )
		return -1;

	datasize = self->descriptor->byte_size;

	offsets  = (unsigned long *)PyArray_DATA(self->descriptor->offsets);
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
	unsigned long *offsets;
	int ret = -1;

	if (!value) {
		/* delete slice */
		PyErr_SetString(PyExc_TypeError, "genoarray objects do not support item deletion");
		return -1;
	}

	if( genoarray_checkstate(self) == -1 )
		return -1;

	seq = NULL;
	datasize = self->descriptor->byte_size;

	offsets  = (unsigned long *)PyArray_DATA(self->descriptor->offsets);
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

	if (PySlice_GetIndicesEx(slice, n, &start, &stop, &step, &slicelength) < 0)
		goto error;

	if( PyList_CheckExact(value) || PyTuple_CheckExact(value) || (PyObject *)self == value )
	{
		seq = PySequence_Fast(value,"cannot convert slice into a sequence");
		if (!seq) goto error;

		m = PySequence_Fast_GET_SIZE(seq);

		if (m != slicelength) {
			PyErr_Format(PyExc_ValueError, "attempt to assign sequence of size %zd to slice of size %zd",
				     m, slicelength);
			goto error;
		}

		if (!slicelength) goto done;
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
		if (seq == NULL)
			goto error;
		iternext = *seq->ob_type->tp_iternext;

		for (i=start, j=0; j<slicelength; i+=step, j++)
		{
			PyObject *geno = iternext(seq);
			if(!geno)
			{
				if (PyErr_Occurred() && PyErr_ExceptionMatches(PyExc_StopIteration))
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
	if (PyIndex_Check(item)) {
		Py_ssize_t i = PyNumber_AsSsize_t(item, PyExc_IndexError);
		if (i == -1 && PyErr_Occurred())
			return -1;

		return genoarray_ass_item(self, i, value);
	}
	else if (PySlice_Check(item))
		return genoarray_ass_slice(self, (PySliceObject *)item, value);
	else {
		PyErr_SetString(PyExc_TypeError, "indices must be integers");
		return -1;
	}
}

static PyObject*
genoarray_descriptor_get(GenotypeArrayObject *self, void *closure)
{
	Py_INCREF(self->descriptor);
	return (PyObject*)self->descriptor;
}

static PyObject*
genoarray_data_get(GenotypeArrayObject *self, void *closure)
{
	PyArrayObject *data;
	int dims;

	if( genoarray_checkstate(self) == -1 )
		return NULL;

	dims = self->descriptor->byte_size;
	data = (PyArrayObject*)PyArray_FromDimsAndData(1,&dims,PyArray_UBYTE,self->data);
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

	if (!new_data)
	{
		PyErr_SetString(PyExc_ValueError,"cannot delete genoarray data");
		return -1;
	}

	if (!PyArray_CheckExact(new_data))
	{
		PyErr_SetString(PyExc_TypeError,"new data must be ndarray type");
		return -1;
	}

	if(PyArray_NBYTES(new_data) != self->descriptor->byte_size)
	{
		PyErr_SetString(PyExc_ValueError,"new data must be same size as current");
		return -1;
	}

	memcpy(self->data, PyArray_DATA(new_data), self->descriptor->byte_size);
	return 0;
}

static Py_ssize_t
genoarray_getreadbuf(GenotypeArrayObject *self, Py_ssize_t index, const void **ptr)
{
	if ( index != 0 ) {
		PyErr_SetString(PyExc_SystemError,
				"accessing non-existent genoarray segment");
		return -1;
	}
	*ptr = (void *)self->data;
	return genoarray_length(self);
}

static Py_ssize_t
genoarray_getwritebuf(GenotypeArrayObject *self, Py_ssize_t index, const void **ptr)
{
	PyErr_SetString(PyExc_TypeError,
			"Cannot use string as modifiable buffer");
	return -1;
}

static Py_ssize_t
genoarray_getsegcount(GenotypeArrayObject *self, Py_ssize_t *lenp)
{

	/* Return one segment and the length in bytes if lenp is not NULL */
	if(lenp)
	{
		Py_ssize_t len = genoarray_length(self);
		if(len == -1 && PyErr_Occurred())
		{
			PyErr_Clear();
			len = 0;
		}

		*lenp = len;
	}
	return 1;
}

static Py_ssize_t
genoarray_getcharbuf(GenotypeArrayObject *self, Py_ssize_t index, const char **ptr)
{
	PyErr_SetString(PyExc_TypeError,
			"Cannot use genoarray as a character buffer");
	return -1;
}

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

PyTypeObject GenotypeArrayType = {
	PyObject_HEAD_INIT(&PyType_Type)
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
	PyObject_GenericGetAttr,		/* tp_getattro       */
	0,					/* tp_setattro       */
	&genoarray_as_buffer,			/* tp_as_buffer      */
	Py_TPFLAGS_DEFAULT|Py_TPFLAGS_HAVE_GC,	/* tp_flags          */
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
	PyObject_GC_Del,			/* tp_free           */
};

static PyMethodDef genoarraymodule_methods[] = {
	{NULL}  /* Sentinel */
};

PyMODINIT_FUNC
init_genoarray(void)
{
	PyObject *m;

	import_array();

	if (PyType_Ready(&GenotypeType) < 0)
		return;

	if (PyType_Ready(&GenotypeArrayType) < 0)
		return;

	if (PyType_Ready(&UnphasedMarkerModelType) < 0)
		return;

	if (PyType_Ready(&GenotypeArrayDescriptorType) < 0)
		return;

	m = Py_InitModule3("_genoarray", genoarraymodule_methods,
		"Genotype array module");

	Py_INCREF(&GenotypeArrayType);
	PyModule_AddObject(m, "GenotypeArrayDescriptor",
	                       (PyObject *)&GenotypeArrayDescriptorType);
	PyModule_AddObject(m, "UnphasedMarkerModel",
	                       (PyObject *)&UnphasedMarkerModelType);
	PyModule_AddObject(m, "GenotypeArray", (PyObject *)&GenotypeArrayType);
	PyModule_AddObject(m, "Genotype",      (PyObject *)&GenotypeType);
}
