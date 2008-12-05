================================================================
:mod:`convert.from_lbd` --- Import an Illumina LBD file into GLU
================================================================

.. module:: convert.from_lbd
   :synopsis: Import an Illumina LBD file into GLU

.. module:: from_lbd
   :synopsis: Import an Illumina LBD file into GLU

Command line
------------

Usage::

  glu convert.from_lbd [options] lbdfile...

Options:

  -h, --help            show this help message and exit
  -F NAME, --outformat=NAME
                        Output genotype format
  -G REP, --outgenorepr=REP
                        Output genotype representation
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
  --renamesamples=FILE  Rename samples from a file containing rows of original
                        name, tab, new name
  --renameloci=FILE     Rename loci from a file containing rows of original
                        name, tab, new name
  --renamealleles=FILE  Rename alleles based on file of locus name, tab, old
                        alleles (comma separated), tab, new alleles (comma
                        separated)
  --ordersamples=FILE   Order samples based on the order of names in FILE
  --orderloci=FILE      Order loci based on the order of names in FILE
  -o FILE, --output=FILE
                        Output genotype file name
  -m FILE, --abmap=FILE
                        Mappings from A and B probes to other allele codings
  -M FILE, --manifest=FILE
                        Illumina manifest file
  -s T, --targetstrand=T
                        Target strand based on Illumina manifest file: top,
                        bottom, forward, reverse, customer (default),
                        anticustomer, design, antidesign
  -t N, --gcthreshold=N
                        Genotypes with GC score less than N set to missing
  -w, --warnings        Emit warnings and A/B calls for SNPs with invalid
                        manifest data
  --samplestats=FILE    Output per sample average GC statistics to FILE
  --locusstats=FILE     Output per locus average GC statistics to FILE
