===========================================================================
:mod:`qc.mendel_check` --- Check family data for non-Mendelian transmission
===========================================================================

.. module:: qc.mendel_check
   :synopsis: Check family data for non-Mendelian transmission

.. module:: mendel_check
   :synopsis: Check family data for non-Mendelian transmission

The :mod:`qc.mendel_check` module checks genotype data from related
individuals and detects patterns of non-Mendelian transmission.

Usage::

  glu qc.mendel_check [options] file

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
  -o FILE, --output=FILE
                        Output Mendelian concordance by sample
  -O FILE, --locout=FILE
                        Output Mendelian concordance by locus
