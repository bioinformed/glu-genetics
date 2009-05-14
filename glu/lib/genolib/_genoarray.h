/*
File:          _genoarray.h

Authors:       Kevin Jacobs (jacobske@mail.nih.gov)

Created:

Abstract:      fast implementation of a bit-packed genotype array

Requires:      Python 2.5, glu

Revision:      $Id$

Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.
See GLU license for terms by running: glu license
*/

#ifndef GLU_GENOLIB_GENOARRAY
#define GLU_GENOLIB_GENOARRAY

#include "Python.h"
#include "numpy/arrayobject.h"

/* Genotype Array objects declarations */

typedef struct {
	PyObject_HEAD
	PyObject      *models;
	PyArrayObject *offsets;
	unsigned int  bit_size;
	unsigned int  max_bit_size;
	unsigned int  byte_size;
	unsigned int  homogeneous;
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
	unsigned short bit_size;
} UnphasedMarkerModelObject;

typedef struct {
	PyObject_HEAD
	UnphasedMarkerModelObject *model;
	unsigned int   index;
	unsigned int   allele1_index;
	unsigned int   allele2_index;
} GenotypeObject;

typedef GenotypeObject * (*cmergefunc_t)(UnphasedMarkerModelObject *model, PyObject *genos, Py_ssize_t *status);

typedef struct {
	PyObject_HEAD
	cmergefunc_t	cmergefunc;
	PyObject	*pymergefunc;
	PyObject	*samplestats;
	PyObject	*locusstats;
	int		trackstats;
} GenotypeMergerObject;

/* Forward declaration */
PyTypeObject GenotypeArrayDescriptorType;
PyTypeObject UnphasedMarkerModelType;
PyTypeObject GenotypeArrayType;
PyTypeObject GenotypeType;
PyTypeObject GenotypeMergerType;

#define GenotypeArray_Check(op)                PyObject_TypeCheck(op, &GenotypeArrayType)
#define GenotypeArray_CheckExact(op)           ((op)->ob_type == &GenotypeArrayType)
#define UnphasedMarkerModel_Check(op)         (((op)->ob_type == &UnphasedMarkerModelType) || PyObject_TypeCheck(op, &UnphasedMarkerModelType))
#define UnphasedMarkerModel_CheckExact(op)     ((op)->ob_type == &UnphasedMarkerModelType)
#define Genotype_CheckExact(op)                ((op)->ob_type == &GenotypeType)
#define GenotypeArrayDescriptor_CheckExact(op) ((op)->ob_type == &GenotypeArrayDescriptorType)
#define GenotypeMerger_CheckExact(op)          ((op)->ob_type == &GenotypeMergerType)

typedef int (*geno_foreach)(Py_ssize_t i, GenotypeObject *geno, void *state);

int for_each_genotype_genoarray(GenotypeArrayObject *genos, geno_foreach func, void *state);
int for_each_genotype(PyObject *genos, geno_foreach func, void *state);

PyObject *count_haplotypes(PyObject *self, PyObject *args);
PyObject *estimate_ld(PyObject *self, PyObject *args);

/* Exceptions */
PyObject *GenotypeLookupError;
PyObject *GenotypeRepresentationError;

#endif
