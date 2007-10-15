========================================================================
:mod:`qc.dupcheck` --- Find duplicate samples within a genotype data set
========================================================================

.. module:: qc.dupcheck
   :synopsis: Find duplicate samples within a genotype data set

.. module:: dupcheck
   :synopsis: Find duplicate samples within a genotype data set

The dupcheck module looks at genotype data and identifies three categories
of duplicates:

- expected duplicates

- unexpected duplicates

- expected duplicates not found.

The input can be in hapmap, ldat, sdat, triple, or genotriple format. The
user can indicate a threshold value (percent) for determining identity
between two individuals, and the minimum number of genotypes to be
considered informative.
