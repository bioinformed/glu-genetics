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

  glu genedb.find_regions [options] file

Options:
  -h, --help            show this help message and exit
  -g NAME, --genedb=NAME
                        Genedb genome annotation database name or file
  -u N, --upbases=N     upstream margin in bases (default=20000)
  -d N, --downbases=N   downstream margin in bases (default=10000)
  -U N, --upsnps=N      maximum number of upstream SNPs (default=0 for no
                        limit)
  -D N, --downsnps=N    maximum number of downstream SNPs (default=0 for no
                        limit)
  -o FILE, --outfile=FILE
                        output file name, '-' for standard out
