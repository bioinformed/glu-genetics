/*
File:          row.c

Author:        Kevin Jacobs (jacobs@theopalgroup.com)

Created:       June 3, 2002

Purpose:       A faster and better database result object

Compatibility: Python 2.2+

Requires:

Revision:      $Id: row.c 244 2006-06-22 21:26:39Z jacobske $

Copyright (c) 2001,2002 The OPAL Group.

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
#include "structmember.h"
#include "fields.h"
#include "row.h"

#define DIRECT
#undef  DIRECT_CMP
#define DIRECT_CONCAT

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

typedef struct {
	PyObject_HEAD
	PyObject* fields;
} rowobject;

static PyMemberDef row_members[] = {
	{"fields", T_OBJECT, offsetof(rowobject, fields), 0, "the fields contained within this row"},
	{0}
};

static int
row_init(PyObject* self, PyObject* args, PyObject* kwds)
{
	rowobject* row    = (rowobject*)self;
	PyObject*  fields = Py_None;

	if (!PyArg_ParseTuple(args, "|O", &fields))
		return -1;

#if 0
	if (fields && fields != Py_None && !PyFields_Check(fields))
		return -1;
#endif

	if(!fields)
		fields = Py_None;

	Py_XINCREF(fields);
	row->fields = fields;
	return 0;
}

static void
row_dealloc(PyObject* self)
{
	rowobject* row = (rowobject*)self;

	Py_XDECREF(row->fields);
	row->fields = NULL;
	self->ob_type->tp_free(self);
}

static int
row_traverse(PyObject* self, visitproc visit, void* arg)
{
	rowobject* row = (rowobject*)self;
	int err;

#define VISIT(SLOT) \
	if (SLOT) { \
		err = visit((PyObject*)(SLOT), arg); \
		if (err) \
			return err; \
	}

	VISIT(row->fields);

#undef VISIT

	return 0;
}


static PyObject*
row_repr(PyObject* self)
{
	rowobject* row = (rowobject*)self;

	if (row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_repr(row->fields);
#else
		return PyObject_Repr(row->fields);
#endif
	else
		return PyString_FromString("<unbound row>");
}

static PyObject*
row_richcmp(rowobject* row, PyObject* other, int op)
{
	if(row->fields && PyFields_Check(row->fields))
	{
		if( PyRow_Check(other) )
		{
			rowobject* other_row = (rowobject*)other;
			if( other_row->fields && PyFields_Check(other_row->fields) )
				return PyObject_RichCompare(row->fields, other_row->fields, op);
		}
		else
			return PyObject_RichCompare(row->fields, other, op);
	}

	PyErr_SetString(PyExc_ValueError, "cannot compare unbound row object");
	return NULL;
}

/*************************************************************************/

static Py_ssize_t
row_length(rowobject* row)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_mapping->mp_length(row->fields);
#else
		return PyObject_Size(row->fields);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot determine length of unbound row object");
		return -1;
	}
}

static PyObject*
row_subscript(rowobject* row, PyObject* key)
{
	if(row->fields && PyFields_Check(row->fields))
	{
		PyObject* ret;
#ifdef DIRECT
		ret = row->fields->ob_type->tp_as_mapping->mp_subscript(row->fields, key);
#else
		ret = PyObject_GetItem(row->fields, key);
#endif

		if(ret == NULL && !PyInt_Check(key))
			PyErr_SetObject(PyExc_KeyError, key);

		return ret;
	}
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot lookup fields of unbound row object");
		return NULL;
	}
}

PyObject*
db_row_subscript_wrap(PyObject* row, PyObject* args)
{
	PyObject* key;
	if(!PyArg_ParseTuple(args, "O", &key))
		return NULL;
	return row_subscript( (rowobject*)row, key);
}

int
row_ass_sub(rowobject* row, PyObject* key, PyObject* value)
{
	if(row->fields && PyFields_Check(row->fields))
	{
		int ret;
#ifdef DIRECT
		ret = row->fields->ob_type->tp_as_mapping->mp_ass_subscript(row->fields, key, value);
#else
		if(value)
			ret = PyObject_SetItem(row->fields, key, value);
		else
			ret = PyObject_DelItem(row->fields, key);
#endif
		if(ret < 0 && !PyInt_Check(key))
			PyErr_SetObject(PyExc_KeyError, key);

		return ret;
	}
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot assign fields of unbound row object");
		return -1;
	}
}

static PyMappingMethods row_as_mapping = {
	(lenfunc)        row_length,    /*mp_length       */
	(binaryfunc)     row_subscript, /*mp_subscript    */
	(objobjargproc)  row_ass_sub,   /*mp_ass_subscript*/
};

/*************************************************************************/


static PyObject*
row_repeat(rowobject* row, Py_ssize_t count)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_sequence->sq_repeat(row->fields, count);
#else
		return PySequence_Repeat(row->fields, count);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot repeat unbound row object");
		return NULL;
	}
}

#ifndef FIXUP_PYTHON22_ORDERING
static PyObject*
row_item(rowobject* row, Py_ssize_t i)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_sequence->sq_item(row->fields, i);
#else
		return PySequence_GetItem(row->fields, i);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot lookup fields of unbound row object");
		return NULL;
	}
}
#endif

static PyObject*
row_slice(rowobject* row, Py_ssize_t ilow, Py_ssize_t ihigh)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_sequence->sq_slice(row->fields, ilow, ihigh);
#else
		return PySequence_GetSlice(row->fields, ilow, ihigh);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot slice unbound row object");
		return NULL;
	}
}

static int
row_ass_slice(rowobject* row, Py_ssize_t ilow, Py_ssize_t ihigh, PyObject* v)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_sequence->sq_ass_slice(row->fields, ilow, ihigh, v);
#else
		return PySequence_SetSlice(row->fields, ilow, ihigh, v);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot slice unbound row object");
		return -1;
	}
}

static int
row_contains(rowobject* row, PyObject* key)
{
	if(row->fields && PyFields_Check(row->fields))
#ifdef DIRECT
		return row->fields->ob_type->tp_as_sequence->sq_contains(row->fields, key);
#else
		return PySequence_Contains(row->fields, key);
#endif
	else
	{
		PyErr_SetString(PyExc_ValueError, "cannot search unbound row object");
		return -1;
	}
}

static PySequenceMethods row_as_sequence = {
	0,                                      /* sq_length         */
	0,                                      /* sq_concat         */
	(ssizeargfunc)     row_repeat,          /* sq_repeat         */
#ifdef FIXUP_PYTHON22_ORDERING
	0,                                      /* sq_item           */
#else
	(ssizeargfunc)     row_item,            /* sq_item           */
#endif
	(ssizessizeargfunc)row_slice,           /* sq_slice          */
	0,                                      /* sq_ass_item       */
	(ssizessizeobjargproc)row_ass_slice,    /* sq_ass_slice      */
	(objobjproc)       row_contains,        /* sq_contains       */
	0,                                      /* sq_inplace_concat */
	0,                                      /* sq_inplace_repeat */
};

/*************************************************************************/

static PyObject*
row_add(PyObject* o1, PyObject* o2)
{
	rowobject* row = NULL;

	if(!o1 || !o2)
		return NULL;

	if( PyRow_Check(o1) )
	{
		row = (rowobject*)o1;
		o1 = row->fields;
		if(!o1 || !PyFields_Check(o1))
			goto error;
	}

	if( PyRow_Check(o2) )
	{
		row = (rowobject*)o2;
		o2  = row->fields;
		if(!o2 || !PyFields_Check(o2))
			goto error;
	}

	return PyNumber_Add(o1, o2);

error:
	PyErr_SetString(PyExc_ValueError, "cannot concatinate unbound row object");
	return NULL;
}

PyNumberMethods row_as_number = {
	(binaryfunc)row_add,       /* nb_add */
};

/*************************************************************************/

PyTypeObject PyRow_Type = {
	PyObject_HEAD_INIT(NULL)
	0,                         /* ob_size */
	"row",                     /* tp_name */
	sizeof(rowobject),         /* tp_basicsize */
	0,                         /* tp_itemsize */
	/* methods */
	row_dealloc,               /* tp_dealloc */
	0,                         /* tp_print */
	0,                         /* tp_getattr */
	0,                         /* tp_setattr */
	0,                         /* tp_compare */
	row_repr,                  /* tp_repr */
	&row_as_number,            /* tp_as_number */
	&row_as_sequence,          /* tp_as_sequence */
	&row_as_mapping,           /* tp_as_mapping */
	0,                         /* tp_hash */
	0,                         /* tp_call */
	0,                         /* tp_str */
	0,                         /* tp_getattro */
	0,                         /* tp_setattro */
	0,                         /* tp_as_buffer */
	/* tp_flags */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE
	                   | Py_TPFLAGS_CHECKTYPES,
	0,                         /* tp_doc */
	row_traverse,              /* tp_traverse */
	0,                         /* tp_clear */
	(richcmpfunc)row_richcmp,  /* tp_richcompare */
	0,                         /* tp_weaklistoffset */
	0,                         /* tp_iter */
	0,                         /* tp_iternext */
	0,                         /* tp_methods */
	row_members,               /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	row_init,                  /* tp_init */
	0,                         /* tp_alloc */
	0,                         /* tp_new */
	0,                         /* tp_free */
};
