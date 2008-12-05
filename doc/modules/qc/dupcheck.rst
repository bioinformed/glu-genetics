========================================================================
:mod:`qc.dupcheck` --- Find duplicate samples within a genotype data set
========================================================================

.. module:: qc.dupcheck
   :synopsis: Find duplicate samples within a genotype data set

.. module:: dupcheck
   :synopsis: Find duplicate samples within a genotype data set

The :mod:`qc.dupcheck` module compares all possible pairs of samples and
classifies them into categories based on the pairwise genotype concordance
rate:

  * expected duplicates

  * unexpected duplicates

  * expected duplicates not found

  * uninformative

Only comparisons between samples with greater than a minimum number of
concordant genotypes are considered to be informative.  Pairs of samples
that exceed a given concordance threshold are considered duplicates.  Pairs
of samples are expected or unexpected duplicates if the user specifies a
mapping between sample identifiers and subject identifiers.  Expected
duplicates map to the same subject identifier, while unexpected duplicates
map to different subjects or do not appear in the mapping.

The ``--checkexp`` option is an optimization that only checks expected
duplicates indicated by the ``-e``/``--duplicates`` option.  It should only
be used in specialized situations, since unexpected duplicates will not be
detected and these can be a sign of serious data quality problem.

The input can be any GLU input format, though sdat/sbat/PLINK ped formats
are optimal.

Usage::

  glu qc.dupcheck [options] genotypes

Options:

  -h, --help            show this help message and exit
  -f NAME, --format=NAME
                        Input genotype format
  -g REP, --genorepr=REP
                        Input genotype representation
  -l FILE, --loci=FILE  Locus description file and options
  -p FILE, --pedigree=FILE
                        Pedigree description file and options
  -e FILE, --duplicates=FILE
                        Mapping from sample identifier to subject identifier
  --checkexp            Check only expected duplicate pairs
  -o FILE, --output=FILE
                        Output of duplicate check report
  -t N, --threshold=N   Minimum proportion genotype concordance threshold for
                        duplicates (default=0.8)
  -m N, --mingenos=N    Minimum number of concordant genotypes to be
                        considered informative (default=20)

Methods
=======

This algorithm is relatively expensive in terms of both computer memory and
computing time.  Unlike many other modules, all genotypes are materialized
and must fit in the available RAM or virtual memory space.  GLU uses
efficient representations, using as little as two bits per biallelic
autosomal SNP, but for large data sets the amount of memory required can
still be quite large.

This algorithm compares all possible pairs of samples, which results in
quadratic growth in processing time in the number of samples.  Increasing
the number of loci increases the amount of memory requires and processing
time linearly in the number of loci.

Genotype comparisons are only considered between pairs of non-missing
genotypes and concordance requires that both alleles match.  Thus genotype
concordance rates are estimates of the probability of the pair sharing two
alleles identical by state (IBS).

Output
======

======================= ===================================================================
Column                  Description
======================= ===================================================================
SAMPLE1                 name of the first sample
SAMPLE2                 name of the second sample
CONCORDANT_GENOTYPES    count of concordant genotypes
COMPARISONS             count of informative genotype comparisons
CONCORDANCE_RATE        genotype concordance rate
EXPECTED_DUPLICATE      indicator if the pair is expected to be a duplicate
OBSERVED_DUPLICATE      indicator if the pair is found to be a duplicate (based on the
                        specified threshold)
======================= ===================================================================


Example
=======

Run::

    glu qc.dupcheck mydat.sbat -o dupcheck.out

.. seealso::

  :mod:`qc.summary`
    Genotype summary statistics

  :mod:`qc.concordance`
    Compute concordance between two sets of genotypes
