==========================================================
:mod:`qc.summary` --- Genotype summary statistics
==========================================================

.. module:: qc.summary
   :synopsis: Genotype summary statistics

.. module:: summary
   :synopsis: Genotype summary statistics

The :mod:`qc.summary` module provides many simple count-based summary
statistics for genotype data.  These counts of missing and non-missing
genotypes per sample and per locus, missing data rates, heterozygosity, and
tests for deviations from Hardy-Weinberg proportions.

Usage::

  glu qc.summary [options] file

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
  --renamesamples=FILE  Rename samples from a file containing rows of original
                        name, tab, new name
  --renameloci=FILE     Rename loci from a file containing rows of original
                        name, tab, new name
  --hwp                 Test for deviation from Hardy-Weinberg proportions
  -s FILE, --summaryout=FILE
                        Summary output file name
  -o FILE, --locusout=FILE
                        Locus output table file name
  -O FILE, --sampleout=FILE
                        Sample output table file name
