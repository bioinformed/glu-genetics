====================================================================
:mod:`struct.pca` --- Principle Components Analysis on genotype data
====================================================================

.. module:: struct.pca
   :synopsis: Principle Components Analysis on genotype data

.. module:: pca
   :synopsis: Principle Components Analysis on genotype data

Principle Components Analysis (PCA) to find large-scale correlations among
samples based on genotype data

Usage::

  glu struct.pca [options] genotypes

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
                        Output principle components (eigenvectors) to FILE
                        (default is "-" for standard out)
  -O FILE, --vecout=FILE
                        Output eigenvalues and statistics to FILE
  --vectors=N           Output the top N eigenvectors.  Default=10


.. seealso::

  :mod:`struct.admix`
    Estimate genetic admixture proportions
