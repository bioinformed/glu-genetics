*********************
What's New in GLU 1.0
*********************

Version 1.0a2 (2008-07-27)
==========================

* Rate parateters for :mod:`tagzilla` and :mod:`qc.dupcheck` now take
  decimal rates and not integer percentages.

* Fixed a missing import that prevented :mod:`qc.dupcheck` from running.

* Corrected a metadata sequencing bug when recoding or merging genotriple
  streams.

* Corrected a bug in code that selects the optimal genotype merge
  algorithm that affected merging genotriple files (tdat/tbat/PrettyBase).

Version 1.0a1 (2008-07-23)
==========================

* Update version numbers and tag release

* :mod:`assoc.logit1` and :mod:`assoc.linear1` are smarter about dropping
  records with missing data.  Only columns used in the model are checked for
  missing values, which allows use of phenotype files with many more
  variables than will be used in a given analysis.  In addition, the subject
  ID and phenotype columns are now configurable.

* Refactored genedb and related code to search for database files based on
  an optional database name and search path. If not specified, a series of
  standard database names and paths will be explored.

  The following modules no longer take the database name as the first argument:

    * :mod:`genedb.find_snps`
    * :mod:`genedb.find_regions`
    * :mod:`genedb.annotate`

  Instead, a '-g/--genedb' option is provided.  E.g.::

    > glu genedb.annotate -g genome36.3 assoc.txt -o assoc_annotated.txt

  This will look for the genome36.3.db file in the standard GLU genedb paths
  (places like /usr/local/share/genedb/).  Absolute paths are also allowed::

    > glu genedb.annotate -g /path/to/genome36.3.db assoc.txt -o assoc_annotated.txt

* Many documentation improvements

* Minor bug fixes, including an internal issue with the genotype counts in
  :mod:`qc.summary` (r725) and to the PLINK genotype writers (r724,r741).
