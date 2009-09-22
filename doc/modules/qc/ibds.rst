========================================================================
:mod:`qc.ibds` --- Compute IBS and IBD sharing for pairs of samples
========================================================================

.. module:: qc.ibds
   :synopsis: Compute IBS and IBD sharing for pairs of samples

.. module:: ibds
   :synopsis: Compute IBS and IBD sharing for pairs of samples

The :mod:`qc.ibds` module compares all possible pairs of samples and
computes allele sharing identical by state (IBS) and identical by descent
(IBD) assuming a homogeneous population using a method of moments
approximation.  These statistics are useful for testing for duplicates and
close relatives.

The input can be any GLU input format, though sdat/sbat/PLINK ped formats
are optimal.

Usage::

  glu qc.ibds [options] genotypes

Options:

  -h, --help            show this help message and exit
  -f NAME, --informat=NAME
                        Input genotype format
  -g REP, --ingenorepr=REP
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
  --filterfounders      Excludes founders
  --filternonfounders   Excludes non-founders
  --frequencies=FILE    Optional genotype file to estimate allele frequencies
  --includetest=FILE    List of samples to test
  --excludetest=FILE    List of samples not to test
  --testpairs=FILE      File containing a list of pairs to test
  -t N, --threshold=N   Output only pairs with estimated IBD0 sharing less
                        than N (default=0.90)
  -o FILE, --output=FILE
                        output table file name
  -P, --progress        Show analysis progress bar, if possible

Methods
=======

IBS estimation counts each pair of non-missing genotypes from each sample at
each locus.  Pairs of genotypes that share no alleles by state (by allele
name) are counted as IBS0, those that share 1 allele by state are counted as
IBD1, and identical genotypes are counted as IBS2.  The resulting IBS counts
are divided by the total number of non-missing comparisons made and reflect
the proportion of informative genotypes comparisons that share n alleles IBS.

IBD computation assumes knowledge of genotypes frequencies (see the
--frequencies option) or population homogeneity so that genotype
frequencies may be estimated empirically from the genotype data provided.
The resulting IBD probabilties reflect the proportion of informative
genotypes that share n alleles identical by descent.

This algorithm is relatively expensive in terms of both computer memory and
computing time.  Unlike many other modules, all genotypes are materialized
and must fit in the available RAM or virtual memory space.  GLU uses
efficient representations, using as little as two bits per biallelic
autosomal SNP, but for large data sets the amount of memory required can
still be quite large.

This algorithm compares all possible pairs of samples, which results in
quadratic growth in processing time in the number of samples.  Increasing
the number of loci increases the amount of memory requires and processing
time linearly in the number of loci.

Output
======

======================= ===================================================================
Column                  Description
======================= ===================================================================
SAMPLE1                 name of the first sample
SAMPLE2                 name of the second sample
COMPARISONS             count of informative genotype comparisons
IBS0                    probability of sharing 0 alleles identical by state
IBS1                    probability of sharing 1 alleles identical by state
IBS2                    probability of sharing 2 alleles identical by state
IBD0                    probability of sharing 0 alleles identical by descent
IBD1                    probability of sharing 1 alleles identical by descent
IBD2                    probability of sharing 2 alleles identical by descent
PIHAT                   average proportion of alleles shared identical by
                        descent (IBD1/2+IBD2)
======================= ===================================================================


Example
=======

Run::

    glu qc.ibds mydat.sbat -o ibds.out

.. seealso::

  :mod:`qc.summary`
    Genotype summary statistics

  :mod:`qc.dupcheck`
    Find duplicate samples within a genotype data set

  :mod:`qc.concordance`
    Compute concordance between two sets of genotypes
