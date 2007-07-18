/*
 * pqueue - A priority-queue extension for Python using Fibonacci heaps.
 * Copyright (C) 1999 Andrew Snare <ajs@cs.monash.edu.au>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* Note to developers :-
 *   Fibonnaci heaps are nasty. Do _NOT_ attempt to debug this code
 *   unless you really really know what you're doing and have a copy
 *   of "Introduction to Algorithms" (Thomas Cormen et al) by your side.
 */

#include <Python.h>
#include <signal.h>
#include <assert.h>

#define MAXDEG	(64)			/* Maximum degree of any node */
                  /* NOTE: This allows us to handle ~ 2.3723E13 nodes */

#undef	DEBUG				/* Turns on debugging stuff */

#if PY_VERSION_HEX < 0x02050000
typedef int Py_ssize_t;
#define lenfunc inquiry
#endif

struct heapnodeRec {
	struct heapnodeRec *p;		/* Pointer to the parent */
	struct heapnodeRec *child;	/* Pointer to any of the children */
	struct heapnodeRec *left;	/* Pointer to the left sibling */
	struct heapnodeRec *right;	/* Pointer to the right sibling */
	int degree;			/* How many children we have */
	int mark;			/* Lost child since last made child */
	PyObject *priority;		/* Priority associated with node */
	PyObject *data;			/* Data associated with the node */
};
typedef struct heapnodeRec heapnode;

struct heapnodetrampRec {
	heapnode *ptr;			/* Trampoline pointer */
	Py_ssize_t refcount;		/* Reference counter */
};
typedef struct heapnodetrampRec heapnodetramp;

typedef struct {
	PyObject_HEAD
	heapnode *min;			/* Pointer to current minimum */
	Py_ssize_t n;			/* Number of nodes in the pq */
	PyObject *dict;			/* Dictionary of data */
} pqueueobject;

staticforward PyTypeObject PQueuetype;

/* PQueue -- some debugging routines */

#ifdef DEBUG

#define LONG(x)	(PyInt_AS_LONG((x)->priority))

static void
display_children(heapnode *child, int level, heapnode *parent)
{
	/* This is a debugging routine. It displays a child list recursively */
	if (child) {
		heapnode *w = child;
		do {
			PyObject* priority = PyObject_Repr(w->priority);
			PyObject* data = PyObject_Repr(w->data);

			printf("%*s%#x P(%#x) L(%#x) R(%#x) D(%d) M(%d) "
			       "(%s, %s) %#x",
			       level*4, "", w, w->p, w->left, w->right,
			       w->degree, w->mark,
			       PyString_AS_STRING(priority),
			       PyString_AS_STRING(data), w->child);
			Py_DECREF(priority);
			Py_DECREF(data);
			if (w->child ) {
				printf(" -> \n");
				display_children(w->child, level+1, w);
			} else
				printf("\n");
#ifdef DEBUG
			assert( w->p == parent );
#endif /* DEBUG */
			w = w->right;
		} while (child != w);
	}
}

static void
display_pqueue(pqueueobject *self)
{
	/* This is a debugging routine. Attempt to display the heap as best
	   as we can. */
	printf("PQueue at %#x with %d nodes.\n", self, self->n);
	if (self->min != NULL)
	{
		heapnode *x;
		printf("Min -> %#x\n", self->min);

		display_children(self->min, 1, NULL);
	}
}

static int
check_children(heapnode *child, int level)
{
	heapnode *w = child;
	heapnode *p = child->p;
	int count = 0;
	int degree = 0;
	do {
		int result;

		printf("%*sNow checking: %#x-[%ld]\n", level*4, "", w,LONG(w));

		/* Check the parent link is correct */
		if (p != w->p)
			printf("%*s%#x-[%ld]'s parent-link not intact -> %#x-[%ld]\n",
			       level*4, "", w, LONG(w), w->p, LONG(w->p));

		/* Check the left link */
		if (w->left->right != w)
			printf("%*s-L-> %#x-[%ld] -R-> %#x-[%ld]\n", level*4,"",
			       w->left, LONG(w->left),
			       w->left->right, LONG(w->left->right));

		/* Check the right link */
		if (w->right->left != w)
			printf("%*s-R-> %#x-[%ld] -L-> %#x-[%ld]\n", level*4,"",
			       w->right, LONG(w->right),
			       w->right->left, LONG(w->left->right));

		/* Check the heap condition */
		if (w->p != NULL) {
			PyObject_Cmp(w->priority, w->p->priority, &result);
			if (result < 0)
				printf("%*sHeap-condition violated: %#x-[%ld] > %#x-[%ld]\n",
				       level*4,"", w->p,LONG(w->p), w,LONG(w));
		}

		/* Recur to children if they're present */
		if (w->child) {
			printf("%*s(should have %d children -> %#x-[%ld])\n",
			       level*4, "", w->degree,w->child,LONG(w->child));
			if( w->child->p != w )
				printf("%*s(doesn't point back -> %#x-[%ld])\n",
				     level*4,"",w->child->p,LONG(w->child->p));

			count += check_children(w->child, level+1);
		}

		degree++;
		w = w->right;
	} while (w != child);

	/* Assert the degree was correct */
	if( w->p == NULL )
		printf("%*s(no parent, degree information unchecked)\n",
		       level*4, "");
	else if (degree != w->p->degree)
		printf("%*s(%d children, %d expected from parent)\n",
		       level*4, "", degree, w->p->degree);

	count += degree;

	return count;
}

static void
check_heap(pqueueobject *pqp)
{
	if (pqp->min != NULL) {
		int count = check_children(pqp->min,0);

		if (count == pqp->n)
			printf("Hmm... %d nodes expected and accounted for.\n",
			       pqp->n);
		else
			printf("Doh! %d nodes expected, only %d found.\n",
			       pqp->n, count);
	}
}

#endif /* DEBUG */

/* PQueue -- C API -- Constructor */

#define	is_pqueueobject(v)	((v)->ob_type == &PQueuetype)

static pqueueobject *
pqueue_new(void)
{
	pqueueobject *pqp;

	pqp = PyObject_NEW(pqueueobject, &PQueuetype);
	if (pqp == NULL)
		return NULL;

	/* Allocate the dictionary */
	pqp->dict = PyDict_New();
	if (pqp->dict == NULL)
		return NULL;

	/* No minimum to start off with (and also no nodes) */
	pqp->min = NULL;
	pqp->n = 0;

	return pqp;
}

/* PQueue methods */

static void
children_dealloc(heapnode *child)
{
	child->left->right = NULL;
	do {
		heapnode *x = child;
		if (x->child != NULL) {
			x->child->left->right = x->right;
			x->right = x->child;
		}
		Py_DECREF(x->priority);
		Py_DECREF(x->data);
		child = child->right;
		free(x);
	} while( child != NULL );
}

static void
pqueue_dealloc(pqueueobject *pqp)
{
	Py_DECREF(pqp->dict);
	if(pqp->min != NULL)
		children_dealloc(pqp->min);
	PyObject_Del(pqp);
}

static PyObject *
pqueue_insert(pqueueobject *self, PyObject *args)
{
	PyObject *priority, *data, *ptr;
	heapnode *x;
	heapnodetramp *tramp;
	int newmin;

	/* Check the parameters first */
	if (!PyArg_ParseTuple(args, "OO:insert", &priority, &data))
		return NULL;

	/* Retrieve the data if it already exists */
	ptr = PyDict_GetItem(self->dict, data);
	if ((ptr == NULL) && (PyErr_Occurred() != NULL))
		return NULL;

	/* Do the comparison early on to detect errors early */
	Py_INCREF(priority);
	Py_INCREF(data);
	if (self->min != NULL) {
		if(PyObject_Cmp(self->min->priority, priority, &newmin) == -1){
			PyErr_SetString(PyExc_ValueError,
					"unable to compare priority");
			Py_DECREF(priority);
			Py_DECREF(data);
			return NULL;
		}
	}

	/* Then try and allocate the node */
	if (NULL == (x = malloc(sizeof(heapnode)))) {
		PyErr_NoMemory();
		Py_DECREF(priority);
		Py_DECREF(data);
		return NULL;
	}

	/* Now make the CObject and put it in the dictionary (if need be) */
	if (ptr == NULL) {
		PyObject *cobject;
		tramp = malloc(sizeof(heapnodetramp));
		cobject = PyCObject_FromVoidPtr(tramp, free);
		if ((tramp == NULL) || (cobject == NULL) ||
		    (PyDict_SetItem(self->dict, data, cobject) == -1)) {
			/* If things failed, clean up and go home */
			if (tramp && !cobject)
				free(tramp);
			Py_XDECREF(cobject);
			Py_DECREF(priority);
			Py_DECREF(data);
			free(x);
			return NULL;
		}
		Py_DECREF(cobject);	/* Since PyDict_SetItem borrows */
		tramp->ptr = x;
		tramp->refcount = 1;
	} else {
		tramp = (heapnodetramp*)PyCObject_AsVoidPtr(ptr);
		/* CObject already exists, increment the counter */
		tramp->ptr = NULL;
		tramp->refcount++;
	}

	/* Initialize the node structure */
	x->degree = 0;
	x->p = NULL;
	x->child = NULL;
	x->mark = 0;
	x->priority = priority;
	x->data = data;

	/* Insert the node into the root list */
	if (self->min == NULL) {
		self->min = x->left = x->right = x;
	} else {
		x->right = self->min;
		x->left = self->min->left;
		self->min->left->right = x;
		self->min->left = x;
		if (newmin > 0) {
			self->min = x;
		}
	}
	self->n++;

	/* We return None to indicate success */
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject *
pqueue_peek(pqueueobject *self)
{
	/* If empty, return an error */
	if (self->min == NULL) {
		PyErr_SetString(PyExc_IndexError, "nothing in the queue");
		return NULL;
	}

	/* Return a tuple of the current min node */
	return Py_BuildValue("OO", self->min->priority, self->min->data);
}

static inline void
consolidate(pqueueobject *self)
{
#ifdef DEBUG
	printf("Starting consolidate()\n");
#endif /* DEBUG */
	if (self->min != NULL)
	{
		int i;
		heapnode *A[MAXDEG];
		memset(A, 0, sizeof(heapnode*)*MAXDEG);

		/* We break the link to detect when we've gone through all the
		   nodes. This can be done since we only remove nodes and
		   never insert them. */
		self->min->left->right = NULL;
		do {
			heapnode *x = self->min;
			int d = x->degree;

			/* Advance to the next one while we can */
			self->min = self->min->right;
#ifdef DEBUG
			printf("Looking at %#x-[%ld] (d=%d).\n", x, LONG(x), d);
#endif /* DEBUG */
			while( A[d] != NULL ) {
				int cmpresult;
				heapnode *y = A[d];
#ifdef DEBUG
				printf("Doh! %#x-[%ld] already has d=%d.\n",
				       y, LONG(y), d);
#endif /* DEBUG */
				/* This _should_ never fail? */
				PyObject_Cmp(x->priority,
					     y->priority, &cmpresult);
				if (cmpresult > 0) {
					heapnode *t = x;
					x = y;
					y = t;
#ifdef DEBUG
					printf("(and we're bigger)\n");
#endif /* DEBUG */
				}

				/* Make y a child of x */
#ifdef DEBUG
				printf("Making %#x-[%ld] a child of %#x-[%ld].\n",
				       y, LONG(y), x, LONG(x));
#endif /* DEBUG */
				y->p = x;
				if (x->child == NULL)
					x->child = y->left = y->right = y;
				else {
					y->right = x->child;
					y->left = x->child->left;
					x->child->left->right = y;
					x->child->left = y;
				}
				x->degree++;
				/* Mark y as having been made a child */
				y->mark = 0;

				A[d++] = NULL;
			}
			A[d] = x;
#ifdef DEBUG
			printf("Storing %#x-[%ld] at d=%d.\n", x, LONG(x), d);
#endif /* DEBUG */
		} while (self->min != NULL);

#ifdef DEBUG
		printf("Time to build the root list....\n");
#endif /* DEBUG */
		for(i=0; i<MAXDEG; i++)
			if (A[i] != NULL) {
#ifdef DEBUG
				printf("Oooh.. found %#x-[%ld] with d=%d.\n",
				       A[i], LONG(A[i]), i);
#endif /* DEBUG */
				/* Insert the node into the root list */
				if (self->min == NULL) {
					self->min = A[i]->left =
						A[i]->right = A[i];
				} else {
					int newmin;
					A[i]->right = self->min;
					A[i]->left = self->min->left;
					self->min->left->right = A[i];
					self->min->left = A[i];

					/* Check to see if we have a new min */
					PyObject_Cmp(self->min->priority,
						     A[i]->priority, &newmin);
					if (newmin > 0) {
						self->min = A[i];
					}
				}
			}
	}
}

static PyObject *
pqueue_pop(pqueueobject *self, PyObject *args)
{
	PyObject *ret;
	heapnode *min;
	heapnode *child;
	heapnodetramp *tramp;

	/* Check the parameters first */
	if (!PyArg_ParseTuple(args, ":pop"))
		return NULL;

	/* Enter the debugger... */
	/* raise(SIGTRAP); */

#ifdef DEBUG
	check_heap(self);
#endif /* DEBUG */

	min = self->min;

	/* If empty, return an error */
	if (min == NULL) {
		PyErr_SetString(PyExc_IndexError, "nothing in the queue");
		return NULL;
	}

#ifdef DEBUG
	printf("Removing min from the root list...\n");
#endif /* DEBUG */

	/* Remove the children of the min-node and put them on root list */
	child = min->child;
	if (child != NULL) {
		heapnode *c = child;
		/* Set each child to have no parent */
		do {
			c->p = NULL;
			c = c->right;
		} while(c != child);

		/* Now put the child-list on the root list */
		min->left->right = child;
		child->left->right = min;
		c = child->left;
		child->left = min->left;
		min->left = c;
	}

	/* Now pull the min-node out of the root list */
	min->left->right = min->right;
	min->right->left = min->left;

	/* Check if we're the only node in the root list */
	if (min == min->right)
		self->min = NULL;
	else {
		self->min = min->right;
#ifdef DEBUG
		check_heap(self);
#endif /* DEBUG */
		consolidate(self);
	}
	self->n--;

	/* Reduce the refcount for the data */
	tramp = PyCObject_AsVoidPtr(PyDict_GetItem(self->dict, min->data));
	if (--(tramp->refcount) == 0) {
		PyDict_DelItem(self->dict, min->data);
	}

	/* Build the return tuple, and de-allocate the node */
	ret = Py_BuildValue("OO", min->priority, min->data);
	Py_DECREF(min->priority);
	Py_DECREF(min->data);
	free(min);
	return ret;
}

#ifdef DEBUG
static PyObject *
pqueue_display(pqueueobject *self, PyObject *args)
{
	/* Check the parameters first */
	if (!PyArg_ParseTuple(args, ":display"))
		return NULL;

	display_pqueue(self);

	Py_INCREF(Py_None);
	return Py_None;
}
#endif /* DEBUG */

static Py_ssize_t
pqueue_length(pqueueobject *pqp)
{
	return pqp->n;
}

static PyObject *
pqueue_subscript(pqueueobject *pqp, PyObject *key)
{
	heapnode *hp;
	heapnodetramp *tramp;
	PyObject *cobject = PyDict_GetItem(pqp->dict, key);

	if ((cobject == NULL) ||
	    ((tramp = PyCObject_AsVoidPtr(cobject))->ptr == NULL)) {
		PyErr_SetObject(PyExc_KeyError, key);
		return NULL;
	}

	hp = tramp->ptr;

	Py_INCREF(hp->priority);
	return hp->priority;
}

static int
pqueue_contains(pqueueobject *pqp, PyObject *key)
{
	heapnodetramp *tramp;
	int ok;
	PyObject *cobject = PyDict_GetItem(pqp->dict, key);
	ok = cobject != NULL && (tramp = PyCObject_AsVoidPtr(cobject))->ptr != NULL;
	return ok;
}

static PyObject *
pqueue_has_key(pqueueobject *pqp, PyObject *key)
{
	int ok = pqueue_contains(pqp, key);
	if(ok == -1)
		return NULL;
	return PyBool_FromLong(ok);
}

static inline void
cut(pqueueobject *pqp, heapnode *x, heapnode *y)
{
#ifdef DEBUG
	printf("Starting cut()\n");
#endif /* DEBUG */
	/* Remove x from the child list of y */
	if (x->right == x)	/* Only child */
		y->child = NULL;
	else {
		if (y->child == x)	    /* Is it worth doing this test? */
			y->child = x->right;
		x->right->left = x->left;
		x->left->right = x->right;
	}
	y->degree--;
	/* Put x on the root list */
	x->left = pqp->min->left;
	x->right = pqp->min;
	pqp->min->left->right = x;
	pqp->min->left = x;
	x->p = NULL;
	x->mark = 0;
}

static void
cascading_cut(pqueueobject *pqp, heapnode *y)
{
	heapnode *z = y->p;
#ifdef DEBUG
	printf("Starting cascading_cut()\n");
#endif /* DEBUG */
	if (z != NULL) {
		if (y->mark == 0)
			y->mark = 1;
		else {
			cut(pqp, y, z);
			cascading_cut(pqp, z);
		}
	}
}

static int
decrease_key(pqueueobject *pqp, heapnode *x, PyObject *priority)
{
	/* Assume we've already checked that x->priority <= priority */
	int result = -1;
	heapnode *y = x->p;
#ifdef DEBUG
	printf("Starting decrease_key()\n");
#endif /* DEBUG */
	if (y != NULL) {
		if ((priority != NULL) &&
		    (PyObject_Cmp(priority, y->priority, &result) == -1)) {
			Py_DECREF(priority);
			PyErr_SetString(PyExc_ValueError,
					"unable to compare value");
			return -1;
		}
	}

	/* Throw away the old priority now */
	Py_DECREF(x->priority);
	x->priority = priority;

	/* If we need to move the node up the tree, do it */
	if ((y != NULL) && (result < 0)) {
		cut(pqp,x,y);
		cascading_cut(pqp,y);
	}

	/* If we have a new minimum, make note of it */
	if (priority != NULL)
		PyObject_Cmp(x->priority, pqp->min->priority, &result);
	if (result < 0)
		pqp->min = x;

	return 0;
}

static inline int
delete_key(pqueueobject *pqp, heapnode *x)
{
	PyObject *min;
	/* Now do a reduce on it */
	decrease_key(pqp, x, NULL);		/* NULL == -infinity */
	/* Now put in the Py_None for the priority */
	Py_INCREF(Py_None);
	x->priority = Py_None;
	/* Now do the extract-min to remove the same element */
	min = pqueue_pop(pqp, PyTuple_New(0));
	if (min == NULL)
		return -1;
	/* Throw away the result */
	Py_DECREF(min);
	return 0;
}

static int
pqueue_ass_sub(pqueueobject *pqp, PyObject *data, PyObject *priority)
{
	int result;
	heapnode *hp;
	heapnodetramp *tramp;
	PyObject *cobject = PyDict_GetItem(pqp->dict, data);

	/* Check we could find the node they're referring to */
	if ((cobject == NULL) ||
	    ((tramp = PyCObject_AsVoidPtr(cobject))->ptr == NULL)) {
		if (priority == NULL) {	     /* Deleting non-existent node */
			PyErr_SetObject(PyExc_KeyError, data);
			return -1;
		} else {		     /* Setting non-existent node */
			/* Turn the set into an insert() */
			PyObject *ret =
				pqueue_insert(pqp,
					      Py_BuildValue("OO", priority,
							    data));
			if (ret == NULL)
				return -1;
			Py_DECREF(ret);
			return 0;
		}
	}

	/* Find the node they're talking about */
	hp = tramp->ptr;

	/* Check if we're doing a deletion */
	if (priority == NULL) {
		return delete_key(pqp, hp);
	}

	/* Next check they're reducing the key */
	if (PyObject_Cmp(priority, hp->priority, &result) == -1) {
		PyErr_SetString(PyExc_ValueError, "unable to compare value");
		return -1;
	} else if (result > 0) {
		/* New key is greater than old. Do a delete/insert. */
		int ret = delete_key(pqp, hp);
		if (ret != 0)
			return ret;
		else {
			PyObject *ret =
				pqueue_insert(pqp,
					      Py_BuildValue("OO", priority,
							    data));
			if (ret == NULL)
				return -1;
			Py_DECREF(ret);
			return 0;
		}
	}
#ifdef DEBUG
	check_heap(pqp);
#endif /* DEBUG */

	/* Take ownership of the new priority, and assign it to the node */
	Py_INCREF(priority);
	return decrease_key(pqp, hp, priority);
}

static long
pqueue_nohash(PyObject *self)
{
	PyErr_SetString(PyExc_TypeError, "pqueue objects are unhashable");
	return -1;
}

static PyObject *
pqueue_iter(PyObject *pqp)
{
	return PyObject_GetIter(((pqueueobject*)pqp)->dict);
}

static PyMappingMethods pqueue_as_mapping = {
	(lenfunc)pqueue_length,		/* mp_length */
	(binaryfunc)pqueue_subscript, 	/* mp_subscript */
	(objobjargproc)pqueue_ass_sub,	/* mp_ass_subscript */
};

static PySequenceMethods pqueue_as_sequence = {
	0,					/* sq_length */
	0,					/* sq_concat */
	0,					/* sq_repeat */
	0,					/* sq_item */
	0,					/* sq_slice */
	0,					/* sq_ass_item */
	0,					/* sq_ass_slice */
	(objobjproc)pqueue_contains,		/* sq_contains */
	0,					/* sq_inplace_concat */
	0,					/* sq_inplace_repeat */
};

static PyMethodDef pqueue_methods[] = {
	{"__contains__",(PyCFunction)pqueue_has_key,	METH_O | METH_COEXIST},
	{"insert",	(PyCFunction)pqueue_insert, 	METH_VARARGS},
	{"peek",	(PyCFunction)pqueue_peek,	METH_NOARGS},
	{"pop",		(PyCFunction)pqueue_pop,	METH_VARARGS},
#ifdef DEBUG
	{"display",	(PyCFunction)pqueue_display,	METH_VARARGS},
#endif /* DEBUG */
	{NULL, NULL}			/* Sentinel */
};

static char PQueuetype_Type__doc__[] =
"Priority queues are used as a FIFO buffer, but with the difference that\n\
items in the queue have a priority associated with them. The item with the\n\
lowest priority is always extracted from the list first.\n";

static PyTypeObject PQueuetype = {
	PyObject_HEAD_INIT(NULL)
	0,			/* ob_size */
	"pqueue",		/* tp_name */
	sizeof(pqueueobject),	/* tp_basicsize */
	0,			/* tp_itemsize */
	/* methods */
	(destructor)pqueue_dealloc,	/* tp_dealloc */
	0,			/* tp_print */
	0,			/* tp_getattr */
	0,			/* tp_setattr */
	0,			/* tp_compare */
	0,			/* tp_repr */
	0,			/* tp_as_number */
	&pqueue_as_sequence,	/* tp_as_sequence */
	&pqueue_as_mapping,	/* tp_as_mapping */
	pqueue_nohash,		/* tp_hash */
	0,			/* tp_call */
	0,			/* tp_str */
	0,			/* tp_getattro */
	0,			/* tp_setattro */
	0,			/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE,	/* tp_flags */
	PQueuetype_Type__doc__,	/* tp_doc */
	0,			/* tp_traverse */
	0,			/* tp_clear */
	0,			/* tp_richcompare */
	0,			/* tp_weaklistoffset */
	pqueue_iter,		/* tp_iter */
	0,			/* tp_iternext */
	pqueue_methods,		/* tp_methods */
	0,			/* tp_members */
	0,			/* tp_getset */
	0,			/* tp_base */
	0,			/* tp_dict */
	0,			/* tp_descr_get */
	0,			/* tp_descr_set */
	0,			/* tp_dictoffset */
	0,			/* tp_init */
	0,			/* tp_alloc */
	0,			/* tp_new */
	0,		        /* tp_free */
};

/* PQueue -- Python API -- Constructor */

static PyObject *
pqueue_PQueue(PyObject *self, PyObject *args)
{
	if (!PyArg_ParseTuple(args, ""))
		return NULL;

	return (PyObject *)pqueue_new();
}

static PyMethodDef PQueueMethods[] = {
	{"PQueue",	(PyCFunction)pqueue_PQueue, METH_VARARGS},
	{NULL, NULL}			/* Sentinel */
};

void
initpqueue(void)
{
	PQueuetype.ob_type = &PyType_Type;
	PQueuetype.tp_getattro = PyObject_GenericGetAttr;
	Py_InitModule("pqueue", PQueueMethods);
}
