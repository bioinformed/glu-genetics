====================================================================
:mod:`qc.concordance` --- Check concordance of two sets of genotypes
====================================================================

.. module:: qc.concordance
   :synopsis: Check concordance of two sets of genotypes

.. module:: concordance
   :synopsis: Check concordance of two sets of genotypes

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
