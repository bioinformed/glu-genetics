==========================================================
:mod:`genedb.annotate` --- Add genomic annotation to SNPs
==========================================================

.. module:: genedb.annotate
   :synopsis: Add genomic annotation to SNPs

.. module:: annotate
   :synopsis: Add genomic annotation to SNPs

Add columns of genomic annotation to a file containing a list of SNP names.

Usage::

  glu genedb.annotate [options] genome_database file

Options:

  -h, --help            show this help message and exit
  -c COLUMN, --column=COLUMN
                        Column name or number in which to find SNPs
  -u N, --upstream=N    upstream margin in bases
  -d N, --downstream=N  the downstream margin in bases
  -o FILE, --outfile=FILE
                        name of the output file, '-' for standard out
