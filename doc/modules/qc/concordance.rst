===========================================================================
:mod:`qc.concordance` --- Compute concordance between two sets of genotypes
===========================================================================

.. module:: qc.concordance
   :synopsis: Compute concordance between two sets of genotypes

.. module:: concordance
   :synopsis: Compute concordance between two sets of genotypes

The :mod:`qc.concordance` module can calculate concordance by sample or by
locus, and between research data sets or between one research data set and a
reference data set. Where allele mapping is unknown (i.e. positive or
negative strand), the allele mapping can be identified by finding the
greatest concordance resulting from different mappings.

Usage::

  glu qc.concordance [options] reference comparison...

Options:

  -h, --help            show this help message and exit
  -f FILE, --refformat=FILE
                        The file format for reference genotype data
  -F FILE, --compformat=FILE
                        The file format for other(comparison) genotype data
  -r FILE, --remap=FILE
                        Determine and output the optimal allele mapping based
                        on greatest concordance
  -a FILE, --allelemap=FILE
                        A list of loci to remap the comparison data alleles to
                        the reference data alleles
  -o FILE               Output the concordance statistics by sample to FILE
  -O FILE               Output the concordance statistics by locus to FILE
  --samplemap=FILE      Map the sample ids for the comparison data to the set
                        of ids in the sample equivalence map
  --locusmap=FILE       Map the locus ids for the comparison data to the set
                        of ids in the locus equivalence map
  --sampleeq=FILE       Equivalence mapping between the sample ids from the
                        comparison data and the reference data
  --locuseq=FILE        Equivalence mapping between the locus ids from the
                        comparison data and the reference data

Methods
=======

Compute the overall genotype discordance rates between a set of reference
data and comparison set of data.

Output
======

Locus output
------------

The optional locus output file, obtained from the -O option, contains
marginal counts and statistics for each locus across all samples that appear
in both genotype sets.

======================= ===================================================================
Column                  Description
======================= ===================================================================
REFKEY		        locus name of reference genotype
COMPKEY                 locus name of comparison genotype
CONCORD                 count of concordant genotypes
DISCORD_HET_HET         count of discordant genotypes where both reference and comparison are heterozygous
DISCORD_HET_HOM         count of discordant genotypes where the reference is heterozygous and comparison is homozygous
DISCORD_HOM_HET         count of discordant genotypes where the reference is homozygous and comparison is heterozygous
DISCORD_HOM_HOM         count of discordant genotypes where both reference and comparison are homozygous
CONCORDANCE_RATE        genotype concordance rate between reference and comparison for samples
REF_HWP                 p-value of exact test for deviation from Hardy-Weinberg proportions of reference genotypes
COMP_HWP                p-value of exact test for deviation from Hardy-Weinberg proportions of comparison genotypes
CONCORD_GENO_PAIR       counts of concordant genotypes for reference and comparison samples
DISCORD_GENO_PAIR       counts of discordant genotypes for reference and comparison samples
ALLELE_MAP_CATEGORY     name of allele map applied, if any
ALLELE_MAPS             allele map applied, if any
======================= ===================================================================

Sample output
-------------

The optional sample output file, obtained from the -o option, contains
marginal counts and statistics for each sample across all loci that appear
in both genotype sets.

======================= ===================================================================
Column                  Description
======================= ===================================================================
REFKEY                  sample name of reference genotype
COMPKEY                 sample name of comparison genotype
CONCORD                 count of concordant genotypes
DISCORD_HET_HET         count of discordant genotypes where both reference and comparison are heterozygous
DISCORD_HET_HOM         count of discordant genotypes where the reference is heterozygous and comparison is homozygous
DISCORD_HOM_HET         count of discordant genotypes where the reference is homozygous and comparison is heterozygous
DISCORD_HOM_HOM         count of discordant genotypes where both reference and comparison are homozygous
CONCORDANCE_RATE        genotype concordance rate between reference and comparison for loci
======================= ===================================================================

Example
=======

Run::

    glu qc.concordance ref_sample.lbat comp_sample.lbat -o sample_discord.out -O locus_discord.out

.. seealso::

  :mod:`qc.dupcheck`
    Find expected and unexpected duplicates samples

  :mod:`qc.summary`
    Genotype summary statistics
