=================================================================
:mod:`genedb.find_snps` --- Find SNPs near SNPs, genes or regions
=================================================================

.. module:: genedb.find_snps
   :synopsis: Find SNPs near SNPs, genes or regions


.. module:: find_snps
   :synopsis: Find SNPs near SNPs, genes or regions

Find SNPs near SNPs, genes or regions

Usage::

  glu genedb.find_snps [options] file

Options:
  -h, --help            show this help message and exit
  -g NAME, --genedb=NAME
                        Genedb genome annotation database name or file
  --includeloci=FILE    List of loci to include
  --excludeloci=FILE    List of loci to exclude
  -u N, --upbases=N     upstream margin in bases (default=20000)
  -d N, --downbases=N   downstream margin in bases (default=10000)
  -U N, --upsnps=N      maximum number of upstream SNPs (default=0 for no
                        limit)
  -D N, --downsnps=N    maximum number of downstream SNPs (default=0 for no
                        limit)
  -o FILE, --output=FILE
                        the name of the output file, '-' for standard out
