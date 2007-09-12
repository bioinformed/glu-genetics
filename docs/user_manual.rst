+++++++++++
User Manual
+++++++++++

GLU is a suite of tools for statistical geneticists, epidemiologists,
statisticicians and analysis working with large sets of Single Nucleotide
Polymorphism (SNP) in order to find associations between specific variants
and specific phenotypes.

GLU provides tools to:

* manage of SNP genotype data sets (designed to scale to 10 billion
  genotypes and beyond...)

* many methods to verify genotype data quality

* tests for association between SNP markers and continuous or discrete trait
  phenotypes.

GLU is free and open-source, and is constantly being enhanced.

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

File formats
============

GLU recognizes several tab-delimited ASCII file formats for specifying
genotypes and other data used in GLU operations. The *transform* module
can convert among the multiple formats described below.

Summary
-------

Genotype data:

  triple
    file containing a single sample, locus, and genotype per line

  ldat
    file containing a matrix of genotypes with a locus per row and a sample per column

  sdat
    file containing a matrix of genotypes with a sample per row and a locus per column

Other data:

  list file
    file containing a list of values in a single column

  map file
    file containing two columns of values, used to indicate a mapping from
    the value in one column to the value in the other column

  counts file
    file containing genotype counts (homozygous major allele, heterozygous,
    homozygous minor allele), by locus or by sample - used in the hwp, maf,
    and hets modules.

Triple format
-------------

  The GLU triple format is a tab-delimited ASCII format with exactly three
  columns of data and one row (line) of data per genotype to be stored.
  Rows may appear in any order and no header line is allowed.

  Columns:

  1. Sample name: User specified sample name with no length limit or format restriction.
  2. Locus name : User specified locus name with no length limit or format restriction.
  3. Genotype   : A two character genotype with alleles ' ACGT'. Genotypes
                  containing one blank are treated as hemizygous, while
                  genotpes composed of an empty string or two blanks are
                  treated as missing.

  Example::

    S1  rs12345  AA
    S2  rs54321  GG
    S3  rs12345  AT
    S1  rs54321  AG

LDAT format
-----------

  The GLU LDAT format is a tab-delimited ASCII format that contains a
  matrix of genotypes with a locus per row and a sample per column:

  Row 1: 'ldat' in the first column, each sample name in subsequent columns

  Row 2: locus name in the first column, genotypes for that locus for each sample in the subsequent columns

  Column 1: 'ldat' in the first row, each locus name in subsequent rows

  Column 2: sample name in the first row, genotypes for that sample for each locus in subsequent rows

  Genotypes are coded as a two character string with alleles ' ACGT'.
  Genotypes containing one blank are treated as hemizygous, while
  genotypes composed of an empty string or two blanks are treated as
  missing.

  Example::

    ldat   S1  S2  S3
    rs123  AA  AT
    rs321  AG  GG  AA
    rs555  CC  CC  CC


SDAT format
-----------

  The GLU SDAT format is a tab-delimited ASCII format that contains a matrix of genotypes with a sample per row and a locus per column.

  Row 1: 'sdat' in the first column, each locus name in subsequent columns

  Row 2: sample name in the first column, genotypes for that sanple for each locus in the subsequent columns

  Column 1: 'sdat' in the first row, each sample name in subsequent rows

  Column 2: locus name in the first row, genotypes for that locus for each sample in subsequent rows

  Genotypes are coded as a two character string with alleles ' ACGT'. Genotypes containing one blank are treated as hemizygous, while genotpes composed of an empty string or two blanks are treated as missing.

  Example::

    sdat  rs123  rs321  rs555
    S1  AA  AG  CC
    S2  AT  GG  CC
    S3    AA  CC

GLU Modules
===========

Most GLU functionality is made available via a module, which can be thought
of as a miniture program.  Each modules takes command line parameters and
options that are specific to the intended purpose, although with consistant
naming and semantics between modules.  Thus it is typically straightforward
to switch among similar modules, as most of the parameters will be
identical.

Data management
---------------

* **transform:** genotype data file manipulation

  The **transform** module is a utility for converting genotype data files
  from one format to another. Supported input file formats are hapmap,
  ldat, sdat, triple, or genotriple. For output, all formats except hapmap
  are supported. Genotypes can be merged, when genotypes are discordant,
  consensus genotypes can be determined using several algorithms.
  Concordance statistics can be generated by locus or by sample.

* **split:** split a genotype file into subsets

  The **split** module splits an input matrix into multiple output files
  by row and/or column groupings, based on mappings provided as input
  parameters. The genotype data input format can be hapmap, ldat, sdat,
  triple, or genotriple. Genotype representation is indicated by an input
  parameter.

Genotype quality control
------------------------

* **completion:** assay and sample completion

  The **completion** module looks at genotype data by locus and by sample,
  and determines the percentage of missing values. The user can indicate
  the number of rows and columns previously dropped from the dataset, so
  that overall completion can be computed accurately.

* **dupcheck:** duplicate sample analysis
  The **dupcheck** module looks at genotype data and identifies three categories of duplicates:

  - expected duplicates

  - unexpected duplicates

  - expected duplicates not found.

  The input can be in hapmap, ldat, sdat, triple, or genotriple format.
  The user can indicate a threshold value (percent) for determining
  identity between two individuals, and the minimum number of genotypes to
  be considered informative.

* **concordance:** genotype concordance between data sets

  The **concordance** module can calculate concordance by sample or by
  locus, and between research data sets or between one research data set
  and a reference data set. Where allele mapping is unknown (i.e. positive
  or negative strand), the allele mapping can be identified by finding the
  greatest concordance resulting from different mappings.

  Many input formats are supported. Default format for reference data is
  ldat, and for comparison data is hapmap.

* **hwp:** test for deviations from Hardy-Weinberg proportions

  The **hwp** module measures deviation from Hardy-Weinberg proportion
  based on data from the following input formats: ldat, sdat, triple, or
  counts. Parameters allow the user to select subsets of samples or loci
  for analysis.

Genotype-phenotype association testing
--------------------------------------

* **assoc.logit1:** Single-locus association tests of dichotomous and
  unordered polytomous outcomes with optional covariates

* **assoc.linear1:** Single-locus association tests of continuous outcomes
  with optional covariates
