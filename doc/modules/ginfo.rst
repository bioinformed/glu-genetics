:mod:`ginfo` --- Get information about genotype files
=====================================================

.. module:: ginfo
   :synopsis: Get information about genotype files

:mod:`ginfo` can query information from genotypes files readable by GLU.
Output can include number of loci, number of samples, and detailed metadata
about loci and samples.

Usage::

  glu ginfo [options] genotypes

Options:

  -h, --help            show this help message and exit
  -f NAME, --format=NAME
                        Input genotype format
  -g REP, --genorepr=REP
                        Input genotype representation
  -l FILE, --loci=FILE  Locus description file and options
  -p FILE, --pedigree=FILE
                        Pedigree description file and options
  --filtermissing       Filters out the samples or loci with missing genotypes
  --includesamples=FILE
                        List of samples to include
  --includeloci=FILE    List of loci to include
  --excludesamples=FILE
                        List of samples to exclude
  --excludeloci=FILE    List of loci to exclude
  -z, --lazy            Be lazy and never materialize the genotypes.  Some
                        results may come back unknown
  -o FILE, --output=FILE
                        Output results (default is "-" for standard out)
  --outputloci=FILE     Output the list of loci to FILE
  --outputsamples=FILE  Output the list of samples to FILE
