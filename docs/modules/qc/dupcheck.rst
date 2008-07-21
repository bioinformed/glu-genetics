========================================================================
:mod:`qc.dupcheck` --- Find duplicate samples within a genotype data set
========================================================================

.. module:: qc.dupcheck
   :synopsis: Find duplicate samples within a genotype data set

.. module:: dupcheck
   :synopsis: Find duplicate samples within a genotype data set

The :mod:`qc.dupcheck` module looks at genotype data and identifies three
categories of duplicates:

  * expected duplicates

  * unexpected duplicates

  * expected duplicates not found.

The input can be in hapmap, ldat, sdat, triple, or genotriple format. The
user can indicate a threshold value (percent) for determining identity
between two individuals, and the minimum number of genotypes to be
considered informative.

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
  -T N, --threshold=N   Threshold for the percentage of identity shared
                        between two individuals (default=85)
  -m N, --mingenos=N, --mincount=N
                        Minimum concordant genotypes to be considered
                        informative for duplicate checking