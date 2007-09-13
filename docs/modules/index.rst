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

  datamanagement.rst
  qc/index.rst
  assoc/index.rst

General modules:

 * intro : Basic introductory information
 * intro.overview : More detailed overview of GLU
 * intro.formats : Basic information on GLU file formats
 * intro.quickstart : Just the basics for how to get started
 * list : List major modules

Data management:

 * transform : Genotype data file manipulation
 * split : Split a genotype file into subsets

Quality control:

 * qc.completion : Assay and sample completion
 * qc.dupcheck : Duplicate sample analysis
 * qc.concordance : Genotype concordance between data sets
 * qc.hwp : Test for deviations from Hardy-Weinberg proportions

Genotype-phenotype association testing:

 * assoc.logit1 : Single-locus association tests of dichotomous and unordered polytomous outcomes with optional covariates
 * assoc.linear1 : Single-locus association tests of continuous outcomes with optional covariates
