.. _user-manual-intro:

+++++++++++++++++
User Introduction
+++++++++++++++++

Introduction
============

Whole-genome association studies are generating unprecedented amounts of
genotype data, frequently billions of genotypes per study, and require new
and scalable computational approaches to address the challenges of storage,
management, quality control, and genetic analysis. GLU, a framework and a
software package, was designed around a set of novel conceptual approaches
to address the need for general and powerful tools that can scale to
effectively handle trillions of genotypes.

Key innovations include:

* compressed binary genotype storage

* use of streaming and stackable data transformations that avoid fully
  materializing data sets in main memory

* integration with a high-level scripting language for easy customization and extension

* *support for parallel processing and distributed computing (not yet released)*

GLU's data management features include:

* import, export, merge, and split genotype data among several common
  formats and standards

* filter based on powerful criteria for inclusion and exclusion

* rename and adjust sample and locus metadata

GLU includes descriptive tools for genotype quality assurance, including:

* estimation of assay completion

* calculation of reproducibility and concordance

* verification of known duplicate samples, and detection of unknown
  duplicate samples

* empirical sex determination

* testing for deviations from Hardy-Weinberg proportions, Mendelian
  inheritance patterns, and non-random patterns of missing data.

Analytic tools are provided, including:

* linkage disequilibrium (LD) estimation

* fitting generalized linear models to test for phenotype/genotype association

* a high-performance tagger-like application with the ability to augment
  SNPs from set panels with optimal tag SNPs using flexible criteria
  including design scores from major genotyping vendors.

Viewed as a library or framework, GLU is designed to be highly extensible,
so that it can be easily augmented and customized, and can serve as a
foundation for the rapid development of new applications.
