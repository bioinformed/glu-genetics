+++++++++++++++++++++++++++++++++++++++++++++++++
The Genotype Library and Utilities - Split module
+++++++++++++++++++++++++++++++++++++++++++++++++

The split module, predictably, is a utility to split a genotype file into subsets
based on a mapping from sample and/or locus to a partition name.  In
addition, large data sets can be split into smaller files, based on a
maximum number of data rows per file.

Typical applications include:

  1. splitting loci by chromosome or region of interest

  2. splitting large datasets (whole genome or large choromsomes) into more
     manageable chunks (e.g., no more than 500 loci per file)

  3. splitting samples by population or phenotype

  4. splitting by a combination of 1-3

Input/Output Options
====================

-f  --format
  Input format for genotype data. Values=hapmap, ldat, sdat, trip, or genotriple

-g  --genorepr
  genotype representation.  Values=snp (default), hapmap, marker

-d  --destdir
  Destination directory for output files.  Write to input file directory by default.

--maxrows
  Split matrix output so that each contains at most N rows of data

--locusgroups
  map from locus name to locus group

--samplegroups
  map from samples name to sample group

--defaultsamplegroup
  Default group for any unmapped sample

--defaultlocusgroup
  Default group for any unmapped sample

--template
  Template for names of the output files
