================================================================================
:mod:`ld.filter` --- Sequentially filter a list of SNPs based on an LD threshold
================================================================================

.. module:: ld.filter
   :synopsis: Sequentially filter a list of SNPs based on an LD threshold

.. module:: filter
   :synopsis: Sequentially filter a list of SNPs based on an LD threshold

Sequentially filter a list of SNPs based on an LD threshold.

Usage::

  glu ld.filter [options] genotypes...

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
  -r N, --r2threshold=N
                        Minimum r-squared threshold (default=0.80)
  -m BASES              Maximum distance in bases between loci to apply LD
                        check.  default=200000
  --lheader=LHEADER     Locus header column name or number (default=Locus)
  -L N, --limit=N       Filter the top N loci (default=0 for unlimited)
  -o FILE, --output=FILE
                        Output LD filter results to FILE
