============================================================
:mod:`ld.matrix` --- Generate a matrix of pairwise LD values
============================================================

.. module:: ld.surrogates
   :synopsis: Generate a matrix of pairwise LD values

.. module:: surrogates
   :synopsis: Generate a matrix of pairwise LD values

Introduction
============

This user's guide explains the use of the :mod:`ld.matrix`.  The options for
this module are very similar to that of TagZilla (:mod:`ld.tagzilla`),
except that LD values are printed, no binning is performed, and obligate
includes are not accepted.

Usage::

  glu ld.matrix [options] genotypes...

Input options:

    -f NAME, --informat=NAME
                        Input genotype format
    -g REP, --ingenorepr=REP
                        Input genotype representation
    -l FILE, --loci=FILE
                        Locus description file and options
    -p FILE, --pedigree=FILE
                        Pedigree description file and options
    --filtermissing     Filters out the samples or loci with missing genotypes
    --includesamples=FILE
                        List of samples to include
    --includeloci=FILE  List of loci to include
    --excludesamples=FILE
                        List of samples to exclude
    --excludeloci=FILE  List of loci to exclude
    --filterfounders    Excludes founders
    --filternonfounders
                        Excludes non-founders
    -S FILE, --ldsubset=FILE
                        File containing loci within the region these loci LD
                        will be analyzed (see -d/--maxdist)
    -R RANGE, --range=RANGE
                        Ranges of genomic locations to analyze, specified as a
                        comma seperated list of start and end coordinates
                        "S-E".  If either S or E is not specified, then the
                        ranges are assumed to be open.  The end coordinate is
                        exclusive and not included in the range.

Output options:

    -o FILE, --output=FILE
                        Output file for formatted data
    -M MEASURE, --measure=MEASURE
                        Measure of LD: r2 (default) or D'

Genotype and LD estimation options:

    -a FREQ, --minmaf=FREQ
                        Minimum minor allele frequency (MAF) (default=0.05)
    -c N, --mincompletion=N
                        Drop loci with less than N valid genotypes. Default=0
    --mincompletionrate=N
                        Drop loci with completion rate less than N (0-1).
                        Default=0
    -m D, --maxdist=D   Maximum inter-marker distance in kb for LD comparison
                        (default=200)
    -P p, --hwp=p       Filter out loci that fail to meet a minimum
                        signficance level (pvalue) for a test Hardy-Weinberg
                        proportion (no default)

LD threshold options:

    -d DPRIME, --dthreshold=DPRIME
                        Minimum d-prime threshold to output (default=0)
    -r N, --rthreshold=N
                        Minimum r-squared threshold to output (default=0)
