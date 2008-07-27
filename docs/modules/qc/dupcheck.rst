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

  * expected duplicates not found.

  * uninformative

Only comparisons between samples with greater than a minimum number of
concordant genotypes are considered to be informative.  Pairs of samples
that exceed a given concordance threshold are considered duplicates.  Pairs
of samples are expected or unexpected duplicates if the user specifies a
mapping between sample identifiers and subject identifiers.  Expected
duplicates map to the same subject identifier, while unexpected duplicates
map to different subjects or do not appear in the mapping.

The input can be any GLU input format, though sdat/sbat/PLINK ped formats
are optimal.

Usage::

  glu qc.dupcheck [options] file

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
  -o FILE, --output=FILE
                        Output of duplicate check report
  -t N, --threshold=N   Minimum proportion genotype concordance threshold for
                        duplicates (default=0.8)
  -m N, --mingenos=N    Minimum number of concordant genotypes to be
                        considered informative (default=20)
