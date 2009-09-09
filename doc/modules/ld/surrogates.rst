==========================================================
:mod:`ld.surrogates` --- Find LD surrogates for SNPs
==========================================================

.. module:: ld.surrogates
   :synopsis: Find LD surrogates for SNPs

.. module:: surrogates
   :synopsis: Find LD surrogates for SNPs

Introduction
============

This user's guide explains the use of the :mod:`ld.surrogates`.  The options
for this module are very similar to that of TagZilla (:mod:`ld.tagzilla`),
except that the best LD surrogate is detected for all excludes (-e) and
design entries (-D) that do not meet the required minimum threshold.  As no
binning is performed, obligate includes are not accepted.

Usage::

  glu ld.surrogates [options] genotypes...

Options:

  -h, --help            show this help message and exit

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
    -e FILE, --excludetag=FILE
                        File containing loci that are excluded from being a
                        tag
    -s FILE, --subset=FILE
                        File containing loci that define the subset to be
                        analyzed of the loci that are read
    -R RANGE, --range=RANGE
                        Ranges of genomic locations to analyze, specified as a
                        comma seperated list of start and end coordinates
                        "S-E".  If either S or E is not specified, then the
                        ranges are assumed to be open.  The end coordinate is
                        exclusive and not included in the range.
    -D FILE, --designscores=FILE
                        Read in design scores or other weights to use as
                        criteria to choose the optimal tag for each bin
    --designdefault=N   Default design score for any locus not found in a
                        design file
    -L N, --limit=N     Limit the number of loci considered to N for testing
                        purposes (default=0 for unlimited)

Output options:

    -o FILE, --output=FILE
                        Output tabular LD information for bins to FILE ('-'
                        for standard out)

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
                        Minimum r-squared threshold to output (default=0.95)
