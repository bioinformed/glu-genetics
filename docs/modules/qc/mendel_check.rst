===========================================================================
:mod:`qc.mendel_check` --- Check family data for non-Mendelian transmission
===========================================================================

.. module:: qc.mendel_check
   :synopsis: Check family data for non-Mendelian transmission

.. module:: mendel_check
   :synopsis: Check family data for non-Mendelian transmission

The :mod:`qc.mendel_check` module checks genotype data from related
individuals and detects patterns of non-Mendelian transmission.

Usage::

  glu qc.mendel_check [options] genotypes

Options:

  -h, --help            show this help message and exit
  -f NAME, --informat=NAME
                        Input genotype format
  -g REP, --ingenorepr=REP
                        Input genotype representation
  -l FILE, --loci=FILE  Locus description file and options
  -p FILE, --pedigree=FILE
                        Pedigree description file and options
  --filtermissing       Filters out the samples or loci with missing genotypes
  --includesamples=FILE
                        List of samples to include
  --includeloci=FILE    List of loci to include
  --excludesamples=FILE
                        List of samples to exclude
  --excludeloci=FILE    List of loci to exclude
  -o FILE, --output=FILE
                        Output Mendelian concordance by sample
  -O FILE, --locout=FILE
                        Output Mendelian concordance by locus

Methods
=======

Only parent-offspring groups are evaluated in this version.  These groups
are composed of a single child with either one or two parents, all with
non-missing genotypes. More comprehensive checks are forthcoming in a future
release that will utilize full pedigree data to detect inconsistent
genotypes among arbitrary sets of relatives.

Output
======

Locus output
------------

The optional locus output file, obtained from the -O option, contains
marginal counts and statistics for each locus across all informative
parent-offspring pairs.

======================= ===================================================================
Column                  Description
======================= ===================================================================
LOCUS                   locus name
CONCORDANT              count of consistent genotypes among parent/offspring groups
TOTAL                   total informative parent/offspring groups
RATE                    genotype concordance rate
======================= ===================================================================

Sample output
-------------

The optional sample output file, obtained from the -o option, contains
marginal counts and statistics for each family across all loci.

======================= ===================================================================
Column                  Description
======================= ===================================================================
CHILD                   child name
PARENT1                 parent 1 name
PARENT2                 parent 2 name
CONCORDANT              count of consistent genotypes
TOTAL                   total informative loci for the given group
RATE                    genotype concordance rate
======================= ===================================================================

Example
-------

Run::

    qc.mendel_check mydat.ldat -p pedigree.tsv -o mend_sample.out -O mend_locus.out

.. seealso::

  :mod:`qc.summary`
    Genotype summary statistics

  :mod:`qc.dupcheck`
    Find expected and unexpected duplicates samples

