=======================================================================================
:mod:`assoc.logit1` --- Association models between a logistic outcome and SNP genotypes
=======================================================================================

.. module:: assoc.logit1
   :synopsis: Association models between a logistic outcome and SNP genotypes

.. module:: logit1
   :synopsis: Association models between a logistic outcome and SNP genotypes

Single-locus association tests of dichotomous and unordered polytomous
outcomes with optional covariates and interactions.

Usage::

  glu assoc.logit1 [options] phenotypes genotypes

Options:
  -h, --help            show this help message and exit

  Input options:

    -f NAME, --format=NAME
                        Input genotype format
    -g REP, --genorepr=REP
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
    --renamesamples=FILE
                        Rename samples from a file containing rows of original
                        name, tab, new name
    --renameloci=FILE   Rename loci from a file containing rows of original
                        name, tab, new name
    --fixedloci=FILE    Genotypes for fixed loci that can be included in the
                        model
    --minmaf=N          Minimum minor allele frequency filter
    --mingenos=N        Minimum number of observed genotype filter.
                        default=10

  Analysis options:

    --model=F           General formula for model to fit
    --test=F            Formula terms to test.  Default is to test all
                        genotype effects if a model is specified, otherwise a
                        2df genotype test (GENO(locus)).
    --stats=T           Comma separated list of test statistics to apply to
                        each model.  Supported tests include score, Wald, and
                        likelihood ratio statistics.  Values: score, wald,
                        lrt.
    --scan=NAME         Name of locus over which to scan, used in --model,
                        --test and --display (default=locus)
    --pid=NAME          Subject column name or number in the phenotype file
                        (default=1)
    --pheno=NAME        Phenotype column name or number in the phenotype file
                        (default=2)
    --refalleles=FILE   Mapping of locus name to the corresponding reference
                        allele
    --allowdups         Allow duplicate individuals in the data (e.g., to
                        accommodate weighting or incidence density sampling)

  Output options:

    -o FILE, --output=FILE
                        Output summary results to FILE
    -O FILE, --details=FILE
                        Output detailed results to FILE
    --display=F         Formula terms to display in the summary output table.
                        Defaults to all test terms.
    --detailsmaxp=P     Output detailed results for only pvalues below P
                        threshold
    -v LEVEL, --verbose=LEVEL
                        Verbosity level of diagnostic output.  O for none, 1
                        for some (default), 2 for exhaustive.
    --ci=N              Show confidence interval around each estimate of width
                        N.  Set to zero to inhibit output.  Default=0
