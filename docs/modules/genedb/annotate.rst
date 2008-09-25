==========================================================
:mod:`genedb.annotate` --- Add genomic annotation to SNPs
==========================================================

.. module:: genedb.annotate
   :synopsis: Add genomic annotation to SNPs

.. module:: annotate
   :synopsis: Add genomic annotation to SNPs

Add columns of genomic annotation to a file containing a list of SNP names.

Usage::

  glu genedb.annotate [options] file

Options:
  -h, --help            show this help message and exit
  -g NAME, --genedb=NAME
                        Genedb genome annotation database name or file
  -c COLUMN, --column=COLUMN
                        Column name or number in which to find SNPs
  -u N, --upstream=N    upstream margin in bases (default=20000)
  -d N, --downstream=N  the downstream margin in bases (default=10000)
  -o FILE, --output=FILE
                        name of the output file, '-' for standard out
