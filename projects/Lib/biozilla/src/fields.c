/*
File:          fields.c

Author:        Kevin Jacobs (jacobs@theopalgroup.com)

Created:       June 3, 2002

Purpose:       A faster and better database result object

Compatibility: Python 2.2+

Requires:

Revision:      $Id: fields.c 244 2006-06-22 21:26:39Z jacobske $

Copyright (c) 2001,2002 The OPAL Group.
Copyright (c) 2004 Kevin Jacobs.

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

#include "Python.h"
#include "abstract.h"
#include "stringobject.h"
#include "structmember.h"
#include "fields.h"
#include "row.h"

#define MEMOIZE_STRING_LOWER

#if PY_VERSION_HEX < 0x02030000
#  define FIXUP_PYTHON22_ORDERING
#endif

#if PY_VERSION_HEX < 0x02050000
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#define lenfunc              inquiry
#define ssizeargfunc         intargfunc
#define ssizessizeargfunc    intintargfunc
#define ssizessizeobjargproc intintobjargproc
#endif

/* (Illegally) Copied from typeobject.c */
typedef struct {
	PyTypeObject      type;
	PyNumberMethods   as_number;
#ifdef FIXUP_PYTHON22_ORDERING
	PySequenceMethods as_sequence;
	PyMappingMethods  as_mapping;
#else
	PyMappingMethods  as_mapping;
	PySequenceMethods as_sequence;
#endif
	PyBufferProcs     as_buffer;
	PyObject*         name;
	PyObject*         slots;
	PyMemberDef       members[1];
} etype;

typedef struct {
	PyObject_HEAD
	PyTypeObject* field_class;
	PyObject*     field_descr;
} fieldsobject;

#ifdef MEMOIZE_STRING_LOWER
static PyObject* (*lower)(PyObject*) = NULL;
#endif

static PyMemberDef fields_members[] = {
	{"__fieldclass__", T_OBJECT, offsetof(fieldsobject, field_class),        0, "fields class"},
	{"__fielddescr__", T_OBJECT, offsetof(fieldsobject, field_descr), READONLY, "fields descriptors"},
	{0}
};

/* Returns borrowed reference */
static inline PyTypeObject*
getfieldtype(fieldsobject* self)
{
	int i,n;
	PyObject* mro = self->ob_type->tp_mro;

	if(!mro || !PyTuple_Check(mro))
		return NULL;

	n = PyTuple_GET_SIZE(mro);

	for (i = 0; i < n; i++)
	{
		PyTypeObject* base = (PyTypeObject*)PyTuple_GET_ITEM(mro, i);
		if(base->tp_dict && PyDict_GetItemString(base->tp_dict, "__fieldnames__"))
			return base;
	}
	return NULL;
}

/* Returns borrowed reference */
static inline PyObject*
getslots(fieldsobject* self)
{
	etype* et;
	et = (etype*)self->field_class;
	if(!et || !et->slots || !PyTuple_Check(et->slots))
		return NULL;
	return et->slots;
}

/* Returns new reference */
static inline PyObject*
buildtuple(fieldsobject* self)
{
	Py_ssize_t n;
	Py_ssize_t i;
	PyObject* tuple;
	PyObject* slots;

	slots = getslots(self);

	n = slots? PyTuple_GET_SIZE(slots) : 0;

	tuple = PyTuple_New(n);

	if( !tuple )
		return NULL;

	for(i = 0; i < n; ++i)
	{
		PyObject* name  = PyTuple_GET_ITEM(slots, i);
		PyObject* value = PyObject_GetAttr( (PyObject*)self, name );
		if(!value)
		{
			value = Py_None;
			Py_INCREF(value);
		}
		PyTuple_SET_ITEM(tuple, i, value);
	}

	return tuple;
}

/* Returns new reference */
static inline PyObject*
builddescriptors(fieldsobject* self)
{
	Py_ssize_t n;
	Py_ssize_t i;
	PyObject* tuple = NULL;
	PyObject* slots;

	if(self->field_class->tp_dict)
		tuple = PyDict_GetItemString(self->field_class->tp_dict, "__fielddescr__");

	if(tuple)
	{
		Py_INCREF(tuple);
		return tuple;
	}

	slots = getslots(self);

	n = slots? PyTuple_GET_SIZE(slots) : 0;

	tuple = PyTuple_New(n);

	for(i = 0; i < n; ++i)
	{
		PyObject* name  = PyTuple_GET_ITEM(slots, i);
		PyObject* value = PyObject_GetAttr( (PyObject*)self->ob_type, name);
		if(!value)
		{
			Py_DECREF(tuple);
			return NULL;
		}
		PyTuple_SET_ITEM(tuple, i, value);
	}

	if(self->field_class->tp_dict)
		PyDict_SetItemString(self->field_class->tp_dict, "__fielddescr__", tuple);

	return tuple;
}

static int
fields_init(fieldsobject* self, PyObject* args, PyObject* kwds)
{
	Py_ssize_t i;
	Py_ssize_t nargs;
	Py_ssize_t nslots;
	PyObject* slots;
	PyObject* arg = NULL;

	if(!PyArg_ParseTuple(args, "|O", &arg))
		return -1;

	self->field_class = getfieldtype(self);

	if(!self->field_class)
		return -1;

	Py_INCREF(self->field_class);

	self->field_descr = builddescriptors(self);

	if(!self->field_descr)
		return -1;

	if(!arg || arg == Py_None)
		return 0;

	arg = PySequence_Tuple(arg);

	if(!arg)
		return -1;

	nargs = PyTuple_GET_SIZE(arg);

	slots = getslots(self);
	nslots = slots? PyTuple_GET_SIZE(slots) : 0;

	if(nargs > nslots)
	{
		PyErr_SetString(PyExc_ValueError, "incorrect number of row values");
		goto error;
	}

	for(i = 0; i < nargs; ++i)
	{
		PyObject* name  = PyTuple_GET_ITEM(slots, i);
		PyObject* value = PyTuple_GET_ITEM(arg, i);

		if( PyObject_SetAttr( (PyObject*)self, name, value) != 0 )
			goto error;
	}

#if SETMISSINGTONONE
	for(; i < nslots; ++i)
	{
		PyObject* slot  = PyTuple_GET_ITEM(slots, i);

		if( PyObject_SetAttr( (PyObject*)self, slot, Py_None) != 0 )
			goto error;
	}
#endif

	Py_DECREF(arg);
	return 0;

error:
	Py_DECREF(arg);
	return -1;

}

static void
fields_dealloc(fieldsobject* self)
{
	Py_XDECREF(self->field_class);
	self->field_class = NULL;
	Py_XDECREF(self->field_descr);
	self->field_descr = NULL;
	self->ob_type->tp_free( (PyObject*) self);
}

static int
fields_traverse(fieldsobject* self, visitproc visit, void* arg)
{
	int err;

#define VISIT(SLOT) \
	if (SLOT) { \
		err = visit((PyObject*)(SLOT), arg); \
		if (err) \
			return err; \
	}

	VISIT(self->field_class);
	VISIT(self->field_descr);

#undef VISIT

	return 0;
}

static PyObject*
fields_repr(fieldsobject* self)
{
	PyObject* tuple;
	PyObject* repr = NULL;

	tuple = buildtuple(self);

	if(tuple)
	{
		repr = PyObject_Repr(tuple);
		Py_DECREF(tuple);
	}
	return repr; /* new reference */
}

static PyObject*
fields_richcmp(PyObject* o1, PyObject* o2, int op)
{
	PyObject* result   = NULL;
	PyObject* o1_tuple = NULL;
	PyObject* o2_tuple = NULL;

	/* Assume that an error must already be set if either argument
	   if NULL or else this will default to an internal API error. */
	if(!o1 || !o2)
		goto error;

	if( PyFields_Check(o1) )
		o1_tuple = buildtuple( (fieldsobject*)o1 );
	else
		o1_tuple = PySequence_Tuple(o1);

	if( PyFields_Check(o2) )
		o2_tuple = buildtuple( (fieldsobject*)o2 );
	else
		o2_tuple = PySequence_Tuple(o2);

	if(o1_tuple && o2_tuple)
		result = PyObject_RichCompare(o1_tuple, o2_tuple, op);

	/* Clear errors that occured when converting to tuples,
	   but not from the tuple-richcomparison function */
	else if( PyErr_Occurred() )
		PyErr_Clear();

	if(!result)
		goto error;

	Py_XDECREF(o1_tuple);
	Py_XDECREF(o2_tuple);

	return result;

error:
	Py_XDECREF(o1_tuple);
	Py_XDECREF(o2_tuple);

	/* If an error occured other than converting from tuples,
	   then propogate it upward.  */
	if( PyErr_Occurred() )
		return NULL;

	/* Otherwise, the comparison is not implemented */
	Py_INCREF(Py_NotImplemented);
	return Py_NotImplemented;
}

/*************************************************************************/

static Py_ssize_t
fields_length(fieldsobject* self)
{
	PyObject* slots = getslots(self);
	return slots? PyTuple_GET_SIZE(slots) : 0;
}

staticforward PyObject* fields_item(fieldsobject* self, int i);

static PyObject*
fields_subscript(fieldsobject* self, PyObject* key)
{
	if(!PyInt_Check(key))
		return self->ob_type->tp_getattro( (PyObject*)self, key );
	return fields_item(self, PyInt_AsLong(key));
}

#ifdef FIXUP_PYTHON22_ORDERING
static PyObject*
fields_subscript_wrap(fieldsobject* self, PyObject* args)
{
	PyObject* key;
	if(!PyArg_ParseTuple(args, "O", &key))
		return NULL;
	return fields_subscript(self, key);
}
#endif

staticforward int fields_ass_item(fieldsobject* self, int i, PyObject* value);

static int
fields_ass_sub(fieldsobject* self, PyObject* key, PyObject* value)
{
	if(!value)
		value = Py_None;

	if(!PyInt_Check(key))
		return PyObject_SetAttr( (PyObject*)self, key, value);

	return fields_ass_item(self, PyInt_AsLong(key), value);
}

static PyMappingMethods fields_as_mapping = {
	(lenfunc)        fields_length,    /*mp_length       */
	(binaryfunc)     fields_subscript, /*mp_subscript    */
	(objobjargproc)  fields_ass_sub,   /*mp_ass_subscript*/
};

/*************************************************************************/

static PyObject*
fields_repeat(fieldsobject* self, Py_ssize_t count)
{
	PyObject* tuple;
	PyObject* result = NULL;

	tuple = buildtuple(self);

	if(tuple)
	{
		result = PySequence_Repeat(tuple, count);
		Py_DECREF(tuple);
	}

	return result; /* new reference */
}

static PyObject*
fields_item(fieldsobject* self, int i)
{
	int n;
	PyObject* slots;
	PyObject* descr;
	PyObject* result;
	descrgetfunc get_fun;

	slots = getslots(self);

	n = slots? PyTuple_GET_SIZE(slots) : 0;
	i = (i<0)? n+i : i;

	if( i < 0 || i >= n )
	{
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	descr   = PyTuple_GET_ITEM(self->field_descr, i);
	get_fun = descr->ob_type->tp_descr_get;
	result  = get_fun(descr, (PyObject*)self, (PyObject*)self->field_class);

#if 0
	if(!result)
	{
		PyErr_Clear();
		result = Py_None;
		Py_INCREF(result);
	}
#endif

	return result; /* already new reference */
}

static int
fields_ass_item(fieldsobject* self, int i, PyObject* value)
{
	int n;
	PyObject* slots;
	PyObject* descr;
	descrsetfunc set_fun;

	slots = getslots(self);

	n = slots? PyTuple_GET_SIZE(slots) : 0;
	i = (i<0)? n+i : i;

	if( i < 0 || i >= n )
	{
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return -1;
	}

	if(!value)
		value = Py_None;

	descr   = PyTuple_GET_ITEM(self->field_descr, i);
	set_fun = descr->ob_type->tp_descr_set;

	return set_fun(descr, (PyObject*)self, value);
}

static PyObject*
fields_slice(fieldsobject* self, int ilow, int ihigh)
{
	Py_ssize_t n;
	Py_ssize_t i;
	PyObject* slice;
	PyObject* slots;

	slots = getslots(self);

	n = slots? PyTuple_GET_SIZE(slots) : 0;

	ilow  = (ilow<0)?  0 : ilow;
	ihigh = (ihigh>n)? n : ihigh;

	slice = PyTuple_New(ihigh - ilow);

	for(i = ilow; i < ihigh; ++i)
	{
		PyObject* key   = PyTuple_GET_ITEM(slots, i);
		PyObject* value = PyObject_GetAttr( (PyObject*)self, key);

		if(!value)
		{
			value = Py_None;
			Py_INCREF(value);
		}
		PyTuple_SET_ITEM(slice, i - ilow, value);
	}

	return slice; /* new reference */
}

static int
fields_ass_slice(fieldsobject* self, int ilow, int ihigh, PyObject* slice)
{
	int i;
	int n;
	PyObject* slots;
	int nslice;

	slots = getslots(self);
	n = slots? PyTuple_GET_SIZE(slots) : 0;

	ilow  = (ilow<0)?  0 : ilow;
	ihigh = (ihigh>n)? n : ihigh;

	nslice = ihigh - ilow;

	if(slice)
	{
		slice = PySequence_Tuple(slice);
		nslice = slots? PyTuple_GET_SIZE(slice) : 0;
	}

	if( ihigh - ilow != nslice )
	{
		PyErr_SetString(PyExc_ValueError, "incorrect number of values");
		goto error;
	}

	for(i = 0; i < nslice; ++i)
	{
		PyObject* slot  = PyTuple_GET_ITEM(slots, ilow+i);
		PyObject* value = Py_None;

		if(slice)
			value = PyTuple_GET_ITEM(slice, i);

		if( PyObject_SetAttr( (PyObject*)self, slot, value) != 0 )
			goto error;
	}

	Py_XDECREF(slice);
	return 0;

error:
	Py_XDECREF(slice);
	return -1;
}

static int
fields_contains(fieldsobject* self, PyObject* key)
{
	int n;
	int i;
	int cmp;
	PyObject* slots;

	slots = getslots(self);
	n = slots? PyTuple_GET_SIZE(slots) : 0;

	for(i = 0; i < n; ++i)
	{
		PyObject* name  = PyTuple_GET_ITEM(slots, i);
		PyObject* value = PyObject_GetAttr( (PyObject*)self, name);
		if(value)
		{
			cmp = PyObject_RichCompareBool(key, value, Py_EQ);
			Py_DECREF(value);
			if(cmp > 0)
				return 1;
			else if(cmp < 0)
				return -1;
		}

	}
	return 0;
}

static PySequenceMethods fields_as_sequence = {
	0,                                      /* sq_length         */
	0,                                      /* sq_concat         */
	(ssizeargfunc)     fields_repeat,       /* sq_repeat         */
#ifdef FIXUP_PYTHON22_ORDERING
	0,                                      /* sq_item           */
#else
	(ssizeargfunc)     fields_item,         /* sq_item           */
#endif
	(ssizessizeargfunc)fields_slice,        /* sq_slice          */
	0,                                      /* sq_ass_item       */
	(ssizessizeobjargproc)fields_ass_slice, /* sq_ass_slice      */
	(objobjproc)       fields_contains,     /* sq_contains       */
	0,                                      /* sq_inplace_concat */
	0,                                      /* sq_inplace_repeat */
};

/*************************************************************************/

static PyObject*
fields_add(PyObject* o1, PyObject* o2)
{
	PyObject* result = NULL;
	PyObject* o1_tuple;
	PyObject* o2_tuple;

	if(!o1 || !o2)
		return result;

	if( PyFields_Check(o1) )
		o1_tuple = buildtuple( (fieldsobject*)o1 );
	else
		o1_tuple = PySequence_Tuple(o1);

	if( PyFields_Check(o2) )
		o2_tuple = buildtuple( (fieldsobject*)o2 );
	else
		o2_tuple = PySequence_Tuple(o2);

	if(o1_tuple && o2_tuple)
		result = PySequence_Concat(o1_tuple, o2_tuple);

	Py_XDECREF(o1_tuple);
	Py_XDECREF(o2_tuple);

	return result; /* new reference or NULL*/
}

PyNumberMethods fields_as_number = {
	(binaryfunc)fields_add,       /* nb_add */
};

/*************************************************************************/

static int
fields_setattro(fieldsobject* self, PyObject* key, PyObject* value)
{
	if(!value)
		value = Py_None;

	return PyObject_GenericSetAttr((PyObject*)self, key, value);
}


/***************************************************************************************/

static PyObject*
ifields_getattro(fieldsobject* self, PyObject* key)
{
	PyObject* result = NULL;

	if(!PyString_Check(key))
	{
		PyErr_SetString(PyExc_TypeError, "expected string argument");
		return NULL;
	}

#ifdef MEMOIZE_STRING_LOWER
	key = lower(key);
#else
	key = PyObject_CallMethod(key, "lower", NULL);
#endif

	if(key)
	{
		result = PyObject_GenericGetAttr((PyObject*)self, key);
		Py_DECREF(key);
	}

	return result; /* new reference or NULL */
}

static int
ifields_setattro(fieldsobject* self, PyObject* key, PyObject* value)
{
	int result = -1;

	if(!value)
		value = Py_None;

	if(!PyString_Check(key))
	{
		PyErr_SetString(PyExc_TypeError, "expected string argument");
		return result;
	}

#ifdef MEMOIZE_STRING_LOWER
	key = lower(key);
#else
	key = PyObject_CallMethod(key, "lower", NULL);
#endif

	if(key)
	{
		result = PyObject_GenericSetAttr((PyObject*)self, key, value);
		Py_DECREF(key);
	}

	return result; /* new reference or NULL */
}

/***************************************************************************************/


PyTypeObject PyFields_Type = {
	PyObject_HEAD_INIT(NULL)
	0,                              /* ob_size */
	"fields",                       /* tp_name */
	sizeof(fieldsobject),           /* tp_basicsize */
	0,                              /* tp_itemsize */
	/* methods */
	(destructor)fields_dealloc,     /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_compare */
	(reprfunc)fields_repr,          /* tp_repr */
	&fields_as_number,              /* tp_as_number */
	&fields_as_sequence,            /* tp_as_sequence */
	&fields_as_mapping,             /* tp_as_mapping */
	0,                              /* tp_hash */
	0,                              /* tp_call */
	0,                              /* tp_str */
	0,                              /* tp_getattro */
	(setattrofunc)fields_setattro,  /* tp_setattro */
	0,                              /* tp_as_buffer */
	/* tp_flags */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE
	                   | Py_TPFLAGS_CHECKTYPES,
	0,                              /* tp_doc */
	(traverseproc)fields_traverse,  /* tp_traverse */
	0,                              /* tp_clear */
	(richcmpfunc)fields_richcmp,    /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	0,                              /* tp_iter */
	0,                              /* tp_iternext */
	0,                              /* tp_methods */
	fields_members,                 /* tp_members */
	0,                              /* tp_getset */
	0,                              /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	(initproc)fields_init,          /* tp_init */
	0,                              /* tp_alloc */
	0,                              /* tp_new */
	0,                              /* tp_free */
};


/***************************************************************************************/

PyTypeObject PyIFields_Type = {
	PyObject_HEAD_INIT(NULL)
	0,                              /* ob_size */
	"ifields",                      /* tp_name */
	sizeof(fieldsobject),           /* tp_basicsize */
	0,                              /* tp_itemsize */
	/* methods */
	(destructor)fields_dealloc,     /* tp_dealloc */
	0,                              /* tp_print */
	0,                              /* tp_getattr */
	0,                              /* tp_setattr */
	0,                              /* tp_compare */
	(reprfunc)fields_repr,          /* tp_repr */
	&fields_as_number,              /* tp_as_number */
	&fields_as_sequence,            /* tp_as_sequence */
	&fields_as_mapping,             /* tp_as_mapping */
	0,                              /* tp_hash */
	0,                              /* tp_call */
	0,                              /* tp_str */
	(getattrofunc)ifields_getattro, /* tp_getattro */
	(setattrofunc)ifields_setattro, /* tp_setattro */
	0,                              /* tp_as_buffer */
	/* tp_flags */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE
	                   | Py_TPFLAGS_CHECKTYPES,
	0,                              /* tp_doc */
	(traverseproc)fields_traverse,  /* tp_traverse */
	0,                              /* tp_clear */
	(richcmpfunc)fields_richcmp,    /* tp_richcompare */
	0,                              /* tp_weaklistoffset */
	0,                              /* tp_iter */
	0,                              /* tp_iternext */
	0,                              /* tp_methods */
	fields_members,                 /* tp_members */
	0,                              /* tp_getset */
	&PyFields_Type,                 /* tp_base */
	0,                              /* tp_dict */
	0,                              /* tp_descr_get */
	0,                              /* tp_descr_set */
	0,                              /* tp_dictoffset */
	(initproc)fields_init,          /* tp_init */
	0,                              /* tp_alloc */
	0,                              /* tp_new */
	0,                              /* tp_free */
};

/***************************************************************************************/

#ifdef FIXUP_PYTHON22_ORDERING

static PyMethodDef fields_fixup[] = {
	{"__getitem__", (PyCFunction)fields_subscript_wrap, METH_VARARGS},
	{NULL,                  NULL}            /* Sentinel */
};

static PyMethodDef row_fixup[] = {
	{"__getitem__", (PyCFunction)db_row_subscript_wrap, METH_VARARGS},
	{NULL,                  NULL}            /* Sentinel */
};

int fixup_object(PyTypeObject* type, PyMethodDef* methods)
{
	PyObject* dict = type->tp_dict;

	for( ; methods && methods->ml_name; ++methods)
	{
		PyObject* method = PyDescr_NewMethod(type, methods);
		if(!method || PyDict_SetItemString(dict,methods->ml_name,method) < 0)
		{
			Py_XDECREF(method);
			return -1;
		}
		Py_DECREF(method);
	}
	return 0;
}

#endif

/***************************************************************************************/

static PyMethodDef PyFieldsModule_methods[] = {
	{NULL,                  NULL}            /* Sentinel */
};


DL_EXPORT(void)
initdb_rowc(void)
{
	PyObject* module;
	PyObject* module_dict;

	module = Py_InitModule("db_rowc", PyFieldsModule_methods);

	if(!module)
		return;

	module_dict = PyModule_GetDict(module);

	PyFields_Type.tp_new  = PyType_GenericNew;
	PyIFields_Type.tp_new = PyType_GenericNew;
	PyRow_Type.tp_new     = PyType_GenericNew;

	if(PyType_Ready(&PyFields_Type)  != 0)
		return;
	if(PyType_Ready(&PyIFields_Type) != 0)
		return;
	if(PyType_Ready(&PyRow_Type)     != 0)
		return;

#ifdef FIXUP_PYTHON22_ORDERING
	if( fixup_object(&PyFields_Type,  fields_fixup) != 0)
		return;
	if( fixup_object(&PyIFields_Type, fields_fixup) != 0)
		return;
	if( fixup_object(&PyRow_Type,     row_fixup)    != 0)
		return;
#endif

#ifdef MEMOIZE_STRING_LOWER
	{
		PyObject* method_obj = PyDict_GetItemString(PyString_Type.tp_dict, "lower");
		PyMethodDescrObject* lower_method = (PyMethodDescrObject*)method_obj;
		if(lower_method && lower_method->d_method->ml_meth)
			lower = (unaryfunc)lower_method->d_method->ml_meth;
	}
#endif

	PyModule_AddObject(module, "FieldsBase",  (PyObject*)&PyFields_Type);
	PyModule_AddObject(module, "IFieldsBase", (PyObject*)&PyIFields_Type);
	PyModule_AddObject(module, "RowBase",     (PyObject*)&PyRow_Type);
}
