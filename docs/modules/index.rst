.. _modules-index:

#############################
 GLU Module Reference Manual
#############################

:Release: |version|
:Date: |today|

Genotype Library and Utilities (GLU) is a suite of tools for statistical
geneticists, epidemiologists, statisticicians and analysis working with
large sets of Single Nucleotide Polymorphism (SNP) in order to find
associations between specific variants and specific phenotypes.

GLU provides tools to:

* manage of SNP genotype data sets (designed to scale to 10 billion
  genotypes and beyond...)

* many methods to verify genotype data quality

* tests for association between SNP markers and continuous or discrete trait
  phenotypes.

While the :ref:`user_manual-index` describes much of the high-level
functionality of GLU, this reference guide describes the standard modules
that are distributed with GLU.  These modules are programs that perform
specific functions like data import, transformation, data quality checking,
and analysis.

.. toctree::
  :maxdepth: 2

  intro
  datamanagement
  convert/index
  genedb/index
  qc/index
  assoc/index
  struct/index
  tagzilla/index
  util/index
