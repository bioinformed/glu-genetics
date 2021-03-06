:mod:`transform` --- Genotype data transformations
==================================================

.. module:: transform
   :synopsis: Genotype data transformations

The :mod:`transform` module is the primary data management component in GLU.
It allows users to read one or more genotype data files in any supported
format, apply very general data transformations, and then output a new
genotype data file in any supported format.  These transformations are
extremely general and include:

* converting among the various GLU file formats
* extracting a subset of samples and/or loci
* renaming samples and/or loci
* merging duplicate samples and loci
* merging multiple input files
* altering genotype coding and alleles
* altering the order of samples and/or loci
* all of the above in virtually any combination

:mod:`transform` is limited to outputting a single genotype file.  If
multiple output files are desired, see the :mod:`split` module.

Usage::

  glu transform [options] genotypes...

Options:

  -h, --help            show this help message and exit

  Input/Output Options:

    -o FILE, --output=FILE
                        Output transformed data to FILE (default is "-" for
                        standard out)
    -f NAME, --informat=NAME
                        Input genotype format
    -g REP, --ingenorepr=REP
                        Input genotype representation
    -F NAME, --outformat=NAME
                        Output genotype format
    -G REP, --outgenorepr=REP
                        Output genotype representation
    -l FILE, --loci=FILE
                        Locus description file and options
    -p FILE, --pedigree=FILE
                        Pedigree description file and options

  Genotype Merging and Reporting:

    --merge=METHOD      Genotype merge algorithm and optional consensus
                        threshold used to form a consensus genotypes.
                        Values=unique,unanimous,vote,ordered.  Value may be
                        optionally followed by a colon and a threshold.
                        Default=unanimous
    --samplemerge=FILE  Sample concordance statistics output to FILE
                        (optional)
    --locusmerge=FILE   Locus concordance statistics output to FILE (optional)

  Filtering:

    --filtermissing     Filters out the samples or loci with missing genotypes
    --includesamples=FILE
                        List of samples to include
    --includeloci=FILE  List of loci to include
    --excludesamples=FILE
                        List of samples to exclude
    --excludeloci=FILE  List of loci to exclude

  Transformations:

    --renamesamples=FILE
                        Rename samples from a file containing rows of original
                        name, tab, new name
    --renameloci=FILE   Rename loci from a file containing rows of original
                        name, tab, new name
    --renamealleles=FILE
                        Rename alleles based on file of locus name, tab, old
                        alleles (comma separated), tab, new alleles (comma
                        separated)
    --ordersamples=FILE
                        Order samples based on the order of names in FILE
    --orderloci=FILE    Order loci based on the order of names in FILE

Examples
--------

Convert an LDAT file to an SDAT file, including only those samples listed in
the "controls" file::

  > glu transform samples.ldat --includesamples=controls -o controls.sdat

Convert samples.ldat file to subjects.ldat, renaming samples according
to the mapping in the sampleid2subjectid file, sending sample concordance
statistics to sample_merge_report.txt, and sending locus concordance
statistics to locus_merge_report.txt::

  > glu transform  samples.ldat  --renamesamples=sampleid2subjectid -o subjects.ldat \
                                 --samplemerge=sample_merge_report.txt               \
                                 --locusmerge=locus_merge_report.txt
