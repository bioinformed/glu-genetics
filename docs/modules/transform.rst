:mod:`transform` --- Genotype data transformations
==================================================

.. module:: transform
   :synopsis: Genotype data transformations

.. index:: module: split

The Transform module performs various transformations on genotype data
files.  These transformations are extremely general and include:

* converting between the various GLU file formats
* extracting subsets of samples and/or loci
* renaming samples and/or loci
* merging duplicate samples and loci
* merging multiple input files
* altering genotype coding and alleles
* altering the order of samples and/or loci
* all of the above in virtually any combination

Transform is limited to outputting a single genotype file.  If multiple
output files are desired, see the :mod:`split` module.

Input/Output Options
--------------------

==================== ================================================
Option               Description
-------------------- ------------------------------------------------
\-o   --output       Output transformed data to FILE(default is "-" for standard out)
\-f   --informat     Input format for genotype data. Values=hapmap, ldat, sdat, trip or genotriple
\-F   --outformat    Output format for genotype data. Default is informat, cannot be hapmap
\-g   --ingenorepr   Input genotype representation.  Values=snp (default), hapmap, marker
\-G   --outgenorepr  Output genotype representation (see -g/--ingenorepr).  Default is ingenorepr
\-l   --limit        Limit the number of rows of data to N for testing purposes
==================== ================================================

Genotype Merging and Reporting Options
--------------------------------------

================== ================================================
Option             Description
------------------ ------------------------------------------------
--merge            Genotype merge algorithm and optional consensus threshold
                   used to form a consensus genotypes.
                   * Values=vote,ordered,unique
                   * Value may be optionally followed by a colon and a threshold.  Default=vote:1
--samplemerge      Per sample concordance statistics output to FILE (optional)
--locusmerge       Per locus concordance statistics output to FILE (optional)
================== ================================================

Filtering and Renaming Options
------------------------------

====================== ================================================
Option                 Description
---------------------- ------------------------------------------------
\-c   --filtermissing  Filter out the samples or loci with missing genotypes
\-n   --includesamples Include list for only samples to use in the transformation and output
\-u   --includeloci    Include list for only loci to use in the transformation and output
\-x   --excludesamples Exclude a list of samples from the transformation and output
\-e   --excludeloci    Exclude a list of loci from the transformation and output
\-m   --renamesamples  Rename samples from a file containing rows of original name, tab, new name
\-r   --renameloci     Rename loci from a file containing rows of original name, tab, new name
====================== ================================================

Transformations Options
-----------------------

====================  ================================================
Option                Description
--------------------  ------------------------------------------------
\-d  --ordersamples   Order samples based on the order of names in FILE
\-D  --orderloci      Order loci based on the order of names in FILE
\-a  --renamealleles  Rename alleles based on file of locus name, tab,
                      old alleles (comma separated), tab, new alleles (comma separated)
====================  ================================================

Examples
--------

Change an LDAT file to an SDAT file, including only those samples in the controls file::

   glu transform  samples.ldat  --includesamples=controls -o controls.ldat

Change the samples.ldat file to subjects.ldat, renaming samples according
to the mapping in the sampleid2subjectid file, sending sample concordance
statistics to sample_merge_report.txt, and sending locus concordance
statistics to locus_merge_report.txt::

  glu transform  samples.ldat  --renamesamples=sampleid2subjectid -o subjects.ldat \
                               --samplemerge=sample_merge_report.txt               \
                               --locusmerge=locus_merge_report.txt
