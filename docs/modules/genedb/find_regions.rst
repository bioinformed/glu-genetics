========================================================================
:mod:`genedb.find_regions` --- Find metadata on SNPs, genes, and regions
========================================================================

.. module:: genedb.find_regions
   :synopsis: Find metadata on SNPs, genes, and regions


.. module:: find_regions
   :synopsis: Find metadata on SNPs, genes, and regions

Resolve genomic metadata given feature names for SNPs, genes, and bounded
regions

Usage::

  glu genedb.find_regions [options] genome_database file

Options:

  -h, --help            show this help message and exit
  -u N, --upbases=N     upstream margin in bases
  -d N, --downbases=N   downstream margin in bases
  -U N, --upsnps=N      maximum number of upstream SNPs
  -D N, --downsnps=N    maximum number of downstream SNPs
  -o FILE, --outfile=FILE
                        output file name, '-' for standard out
