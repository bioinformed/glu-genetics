====================================================================
:mod:`struct.admix` --- Estimate genetic admixture proportions
====================================================================

.. module:: struct.admix
   :synopsis: Estimate genetic admixture proportions

.. module:: admix
   :synopsis: Estimate genetic admixture proportions

This module estimates admixture proportions from a series of presumed
ancestral populations by maximizing a simplified admixture likelihood that
assumes independent loci and that each genotype appears once at minimum (to
avoid trivially zero likelihoods).

This method gives very similar results as STRUCTURE (Pritchard, Stephens &
Donnelly, 2000) for this specific problem, but requires only a small
fraction of the computational time.

Usage::

  glu struct.admix [options] test_genotypes pop1_genotypes pop2_genotypes [pop3_genotypes...]

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
  --labels=LABELS       Population labels (specify one per population
                        separated with commas)
  --model=MODEL         Model for genotype frequencies.  HWP to assume Hardy-
                        Weinberg proportions, otherwise GENO to fit genotypes
                        based on frequency.  (Default=HWP)
  -t N, --threshold=N   Imputed ancestry threshold (default=0.80)
  -o FILE, --output=FILE
                        output table file name
  -P, --progress        Show analysis progress bar, if possible

Methods
=======

Admixture proportions are determined by maximizing the likelihood of
observing each sample in the test set assuming all loci are independent and
that each genotype appears at least once in each assumed ancestral
population.  Genotype frequencies are used directly in the likelihood,
without assuming Hardy-Weinberg proportions.  Thus the penalty for having a
genotype that is not observed in an ancestral population is proportional to
the total number of individuals observed in that population.  In small
samples, this also tends to downwardly bias the estimated genotype
frequencies.

An individual is considered of a given ancestry based on the supplied labels
and estimated admixture coefficients if their coefficient is greater than a
given threshold.

Otherwise, an individual who has no single estimated admixture coefficient
that meets the specified threshold then one of two behaviors result.  If
only one population group exceeds 1-threshold then the ancestry is deemed
'ADMIXED' for that population.  Otherwise, a list of populations with
estimated admixture above 1-threshold is returned.

Output
======

======================= ===================================================================
Column                  Description
======================= ===================================================================
SAMPLE                  sample name
*label 1*               Admixture proportion for the first population (heading
                        name specified by --label)
*label 2*               Admixture proportion for the second population (heading
                        name specified by --label)
...
*label n*               Admixture proportion for the n'th population (heading
                        name specified by --label)
IMPUTED_ANCESTRY        Classification of imputed ancestry based on admixture
                        proportions and minimum threshold (specified by --threshold)
======================= ===================================================================


Example
=======

Run::

    glu struct.admix -P mysamples.sbat CEU.lbat YRI.lbat ASA.lbat --labels=CEU,YRI,ASA -o admix.out

.. seealso::

  :mod:`struct.pca`
    Principle Components Analysis on genotype data
