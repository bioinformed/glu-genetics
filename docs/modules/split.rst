==========================================
:mod:`split` --- Split genotype data files
==========================================

.. module:: split
   :synopsis: Split genotype data files

The split module, predictably, is a utility to split a genotype file into
subsets based on a mapping from sample and/or locus to a partition name.  In
addition, large data sets can be split into smaller files, based on a
maximum number of data rows per file.  Also, see the :mod:`transform` module
to perform other transformations than simple splitting.

Typical applications include:

  1. splitting loci by chromosome or region of interest

  2. splitting large datasets (whole genome or large chromosomes) into more
     manageable chunks (e.g., no more than 500 loci per file)

  3. splitting samples by population or phenotype

  4. splitting by a combination of 1-3

Usage::

  glu split [options] genotypes

Options:

  -h, --help            show this help message and exit
  -f NAME, --format=NAME
                        Input genotype format
  -g REP, --genorepr=REP
                        Input genotype representation
  -l FILE, --loci=FILE  Locus description file and options
  -p FILE, --pedigree=FILE
                        Pedigree description file and options
  -d DESTDIR, --destdir=DESTDIR
                        Destination directory for output files.  Write to
                        input file directory by default.
  --maxrows=N           Split matrix output so that each contains at most N
                        rows of data
  --locusgroups=FILE    map from locus name to locus group
  --samplegroups=FILE   map from samples name to sample group
  --defaultsamplegroup=NAME
                        Default group for any unmapped sample
  --defaultlocusgroup=NAME
                        Default group for any unmapped sample
  --template=NAME       The template for names of the output files
