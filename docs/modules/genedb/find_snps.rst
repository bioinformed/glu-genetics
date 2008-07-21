=================================================================
:mod:`genedb.find_snps` --- Find SNPs near SNPs, genes or regions
=================================================================

.. module:: genedb.find_snps
   :synopsis: Find SNPs near SNPs, genes or regions


.. module:: find_snps
   :synopsis: Find SNPs near SNPs, genes or regions

Find SNPs near SNPs, genes or regions

Usage::

  glu genedb.find_snps [options] genome_database file

Options:

  -h, --help            show this help message and exit
  -u N, --upbases=N     upstream margin in bases
  -d N, --downbases=N   downstream margin in bases
  -U N, --upsnps=N      maximum number of upstream SNPs
  -D N, --downsnps=N    maximum number of downstream SNPs
  -o FILE, --outfile=FILE
                        the name of the output file, '-' for standard out
  -i FILE, --infile=FILE
                        the name of the feature file containing list of
                        features, '-' for standard in
