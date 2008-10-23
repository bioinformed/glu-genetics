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

  glu qc.summary [options] genotypes

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

Methods
=======

The overall genotyping missing rates, which are computed under three
categories:

* Attempted: the number of samples or loci in the genotype file, plus (or
             minus) the number of samples/loci in the include (or exclude)
             list

* Observed:  the actual number of samples or loci present in the genotype
             file

* Non-empty: the number of samples or loci remaining after excluding those
             with only missing data.  These typically represent lack of
             attempted genotypes or complete genotyping failures for given
             samples or loci

Output
======

Summary output
--------------

The optional summary output file, obtained from the from option
-s/--summaryout option, contains counts and statistics across all samples
and loci.  These include counts of the attempted, observed and non-empty
loci and samples, plus completion rates all attempted, observed, and
non-empty samples and loci.

Locus output
------------

The optional locus output file, obtained from the from option -o/--locusout
option, contains marginal counts and statistics for each locus across all
samples.

======================= ===================================================================
Column                  Description
======================= ===================================================================
LOCUS                   locus name
CHROMOSOME              chromosome
LOCATION                chromosome location
STRAND                  allelic strand (+ = forward, - = reverse)
NUM_ALLELES             total number of alleles observed
ALLELES                 alleles observed
ALLELE_COUNTS           count of each allele
MAF                     minor allele frequency for biallelic loci
NUM_GENOTYPES           number of genotypes observed
GENOTYPES               genotypes observed
GENOTYPE_COUNTS         count of each genotype
MISSING_COUNT           count of missing genotypes for all attempted samples
INFORMATIVE_COUNT       counts of observed non-missing genotypes
ATTEMPTED_MISSING_RATE  missing rates calculated for attempted samples
OBSERVED_MISSING_RATE   missing rates calculated for observed samples
NONEMPTY_MISSING_RATE   missing rates calculated for non-empty samples
HW_PVALUE               p-value of exact test for deviation from Hardy-Weinberg proportions
======================= ===================================================================

Sample output
-------------

The optional sample output file, obtained from the from option -O/--sampleout
option, contains marginal counts and statistics for each sample across all
loci.

======================= ===================================================================
Column                  Description
======================= ===================================================================
SAMPLE                  sample name
MISSING_COUNT           count of missing genotypes for all attempted loci
HEMIZYGOTE_COUNT        count of hemizygote genotypes
HOMOZYGOTE_COUNT        count of homozygote genotypes
HETEROZYGOTE_COUNT      count of heterozygote genotypes
INFORMATIVE_COUNT       counts of obeserved genotypes
ATTEMPTED_MISSING_RATE  missing rates calculated for attempted loci
OBSERVED_MISSING_RATE   missing rates calculated for observed loci
NONEMPTY_MISSING_RATE   missing rates calculated for non-empty loci
HETEROZYGOSITY          ratio of heterozygote genotypes to total obeserved genotypes
======================= ===================================================================

Example
-------

Run::

    glu qc.summary mydata.lbat -o locus_summary.out -O sample_summary.out -s summary.out

.. seealso::

  :mod:`qc.dupcheck`
    Find expected and unexpected duplicates samples

  :mod:`qc.concordance`
    Verify genotype concordance with another reference set of genotypes
