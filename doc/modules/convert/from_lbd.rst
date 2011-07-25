================================================================
:mod:`convert.from_lbd` --- Import an Illumina LBD file into GLU
================================================================

.. module:: convert.from_lbd
   :synopsis: Import an Illumina LBD file into GLU

.. module:: from_lbd
   :synopsis: Import an Illumina LBD file into GLU

The :mod:`convert.from_lbd` module reads an Illumina "Locus by DNA" report
file and converts it to a GLU genotype file.  Locus by DNA reports generated
by Illumina's BeadStudio/GenomeStudio software and are in a matrix format
that lists each sample's genotypes in terms of A and B probes and their
associated quality scores.

With the specification of an Illumina manifest file (-M), the A/B genotype
calls can be converted via the '-s/--targetstrand' option to genomic strand,
genomic strand, "top" strand, the customer specified strand, or the assay
design strand.  The complement of each of these strand designations are also
available.  Manifest files are supplied by Illumina and may be in the binary
BPM format or their CSV export format.

Support strand designations:

============  ================================================================
strand        Genotype calls relative to
designation
============  ================================================================
ab            A/B probe calls (default if no manifest or mapping is specified)
top           top genomic strand
bottom        bottom genomic strand
forward       forward genomic strand
reverse       reverse genomic strand
customer      customer strand (default when a manifest file is specified)
anticustomer  complement of the customer strand
design        assay design strand
antidesign    complement of the assay design strand
============  ================================================================

Alternatively, A/B probes to alleles can be specified explicitly via a
user-supplied mapping file (-m/--abmap).

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
  --filterfounders      Excludes founders
  --filternonfounders   Excludes non-founders
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
                        Illumina manifest file (BPM or CSV)
  -s T, --targetstrand=T
                        Target strand based on Illumina manifest file: ab,
                        top, bottom, forward, reverse, customer (default),
                        anticustomer, design, antidesign
  -t N, --gcthreshold=N
                        Genotypes with GC score less than N set to missing
  -w, --warnings        Emit warnings and A/B calls for SNPs with invalid
                        manifest data
  --samplestats=FILE    Output per sample average GC statistics to FILE
  --locusstats=FILE     Output per locus average GC statistics to FILE
