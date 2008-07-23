==========================================================
:mod:`tagzilla` --- Robust and fast SNP tagging
==========================================================

.. module:: tagzilla
   :synopsis: Robust and fast SNP tagging program

Introduction
============

This user's guide explains the use of the :mod:`tagzilla`. TagZilla chooses
tag SNPs for any given set of SNPs with genotypes. The TagZilla SNP
selection algorithm estimates pair-wise r\ :sup:`2` and d' (d-prime) linkage
disequilibrium (LD) statistics on genotype data from unrelated
individuals.It then estimates bins using a greedy maximal approach similar
to that of Carlson et al. (2004), evaluates all the tags for each bin based
upon user-specified criteria, and recommends an optimal tag for each bin if
possible.

The program provides many options, such as minimum MAF (minor allele
frequency), d' threshold and r\ :sup:`2` threshold, include/exclude/subset,
optimization of numbers of bins for fixed sized panels, or thresholds to
reach a desired level of coverage. It also provides options for different
input formats such as HapMap, Linkage, and FESTA. Data in assay design
scores and weighting criteria can be incorporated into the analysis in order
to choosing optimal tags for each bin.

TagZilla options
================

TagZilla provides many options for controlling its execution. It can read in
multiple genotype files that contain data from disjoint genomic regions.
TagZilla allows each input genotype file to utilize different minimum MAF,
r\ :sup:`2`, completion rate and other parameters.

Usage::

  glu tagzilla [options] genotype_file [options] genotype_file


Example::

  glu tagzilla -p pedinfo2sample_CEU.txt -b summary -o pairs.out -O loci.out -D designscore.txt:0.6 \
               -f hapmap genotypes_chr21_CEU.txt.gz


TagZilla supports both short options and long options.

    * Short options: a token that starts with a dash followed by a letter. For example, -p

    * Long options: a token that starts with two dashes followed by a word. For example, --pedfile

The --version option is used to check the program version and --help/-h
option to print out help messages. The options are grouped into four
categories:

    * Genotype and LD estimate options
    * Binning options
    * Input options
    * Output options

The table for each category serves as a quick reference to the options. More
detailed discussions about the purpose and usage of the options are in the
notes below each table.

Genotype and LD estimation options
----------------------------------

    -a FREQ, --minmaf=FREQ
                        Minimum minor allele frequency (MAF) (default=0.05)
    -A FREQ, --minobmaf=FREQ
                        Minimum minor allele frequency (MAF) for obligate tags
                        (defaults to -a/--minmaf)
    -c N, --mincompletion=N
                        Drop loci with less than N valid genotypes (default=0)
    --mincompletionrate=N
                        Drop loci with completion rate less than N percent (0-100)
                        (default=0)
    -m D, --maxdist=D   Maximum inter-marker distance in kb for LD comparison
                        (default=200)
    -P p, --hwp=p       Filter out loci that fail to meet a minimum
                        signficance level (pvalue) for a test Hardy-Weignberg
                        proportion (no default)
	
Notes:

-a/--minmaf and -A/--minobmaf:

  Both options specify the MAF threshold to filter out loci with low MAF
  from the analysis. You can set different thresholds for obligates versus
  any other loci, but the default for -a/--minmaf is 0.05, and the default
  for -A/--minobmaf will take the value set for -a/--minmaf option.

-c/--mincompletion and --mincompletionrate:

  These options are used to drop loci with low number or rate of valid
  genotypes among all the genotyped samples.

-m/--maxdist:

  The linkage between loci usually decreases as the distance between loci
  increases. We won't consider the linkage disequilibrium between two loci
  if the distance between them is greater than the number specified in this
  option.

-P/--hwp:

  This option is used to specify the threshold of P value for the
  Hardy-Weinberg Equilibrium test. If the count of the minor alleles in the
  set of genotypes is less than 1000, TagZilla applies the exact test based
  on Wigginton JE et al. (2005), otherwise it simply uses the standard
  Chi-square test. Loci that fail to meet this threshold are filtered from
  the analysis.

Binning options
---------------

    -d DPRIME, --dthreshold=DPRIME
                        Minimum d-prime threshold to output (default=0)
    -M POPS, --multipopulation=POPS
                        Multipopulation tagging where every N input files
                        represent a group of populations. May be specified as
                        an integer N or a comma separated list of population
                        labels.
    -r N, --rthreshold=N
                        Minimum r-squared threshold to output (default=0.8)
    -t N, --targetbins=N
                        Stop when N bins have been selected (default=0 for
                        unlimited)
    -T N, --targetloci=N
                        Stop when N loci have been tagged (default=0 for
                        unlimited)
    -C crit, --tagcriteria=crit
                        Use the specified criteria to choose the optimal tag
                        for each bin

                        Currently supported tag selection criteria:

                          maxtag: choose the tag having largest minimum-r\ :sup:`2` with any tag snps in the bin

                          maxsnp: choose the tag having largest minimum-r\ :sup:`2` with all snps in the bin

                          avgtag: choose the tag having maximum average-r\ :sup:`2` with non-tag snps in the bin

                          avgsnp: choose the tag having maximum average-r\ :sup:`2` with all snps in the bin

    -z N, --locipertag=N
                        Ensure that bins contain more than one tag per N loci.
                        Bins with an insufficient number of tags will be
                        reduced.
    -Z B, --loglocipertag=B
                        Ensure that bins contains more than one tag per
                        log_B(loci).  Bins with an insufficient number of tags
                        will be reduced.
    --skipbinning       Skip binning step.  Typically used in conjunction with
                        -u/--saveldpairs

Notes:

-C/--tagcriteria:

  Example: -C maxsnp:2

  gives half the weight the each tag that does not meet the maxsnp criteria.

  This option can be used together with the -D/--designscores option to
  specify how the optimal tag should be selected for each bin.
  -C/--tagcriteria provides the weights, and -D/--designscores provides the
  designscores. TagZilla will compute a weighted score and thus determine
  which tag is recommended to the user.

-d/--dthreshold and -r/--rthreshold:

  Both are used as cut-off criteria so that only locus pairs satisfying
  these thresholds are considered in the binning process.

-t/--targetbins and -T/--targetloic:

  Both options are used as stopping criteria. In either case, once the
  criteria are met, Tagzilla produces residual bins instead of maximal bins.

-M/--multipopulation:

  You can specify the number of populations via -M/--multipopulation option.
  Tagzilla uses minLD method if -multimerge hasn't been set to bin the loci
  with genotypes from different populations and thus generate a set of tags
  applicable for all the populations.

-z/--loicpertag and -Z/--loglocipertag:

  Both options control the ratio between the tags and loci.If the size of
  the bin is too large and thus the number of loci per tag is too big, the
  genotype failure on the tag will lead to losing information on lots of
  loci surrogated only by that tag.Instead of picking another candidate tag
  from large bin as in a post-process, TagZilla incorporates this user
  requirement into the binning process and generates bins only satisfying
  these requirements.

Input options
-------------

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
    -e FILE, --excludetag=FILE
                        File containing loci that are excluded from being a
                        tag
    -i FILE, --includeuntyped=FILE
                        File containing loci that are obligatorily tags and
                        untyped (may not cover another obligate locus)
    -I FILE, --includetyped=FILE
                        File containing loci that are obligatorily tags but
                        have been typed (may cover another typed locus)
    -s FILE, --subset=FILE
                        File containing loci to be used in analysis
    -S FILE, --ldsubset=FILE
                        File containing loci within the region these loci LD
                        will be analyzed (see -d/--maxdist)
    -R RANGE, --range=RANGE
                        Ranges of genomic locations to analyze, specified as a
                        comma seperated list of start and end coordinates
                        "S-E".  If either S or E is not specified, then the
                        ranges are assumed to be open.  The end coordinate is
                        exclusive and not included in the range.

                        Example: -R 10000-20000,30000-80000

    -D FILE, --designscores=FILE
                        Read in design scores or other weights to use as
                        criteria to choose the optimal tag for each bin

                        Example: -D designscore1.txt:0.5:1

                        0.5 is minimum threshold for what is considered
                        designable, 1 is the scale. Both are optional for
                        this option entry, the default value for threshold
                        is 0 and the default value for scale is 1.

                        This option can be specified multiple times on the
                        command line.

    --designdefault=N   Default design score for any locus not found in a
                        design file
    -L N, --limit=N     Limit the number of loci considered to N for testing
                        purposes (default=0 for unlimited)



Notes:

-p/--pedfile:

  This option specifies the pedigree file for those genotypes provided in
  the format of Hapmap, Prettybase or raw. It is not meaningful to specify a
  pedigree file when reading genotype or LD data in linkage or FESTA format.
  The genotypes for the non-founders as found in the pedigree file won't be
  considered in the binning process. Note that if the pedigree file is
  incomplete, we assume all the individuals not contained in the pedigree
  file are founders.

-s/--subset:

  Besides providing a file containing the subset of loci to be analyzed, the
  user can also specify a comma separated list of loci from the command
  line. The string value for this option has to start with a colon. For
  example:

  -s :rs12355,rs12365,rs12488


-l/--loci:

  The locus description file for genotype input in linkage format. TagZilla
  reads in the location for each locus from this file.

-i/--includetag and -e/--excludetag:

  Similar to -s/--subset option you can also specify a list of loci as the
  string value for both options in addition to a file name. The specified
  list of loci are either forced in as tags or excluded from being chosen as
  tags for non-excludes.

-D/--designscores:

  This option can be used alone or together with -C/--tagcriteria to choose
  the optimal tag among all the valid tags for a bin.

-L/--limit:

  This option is useful for testing purposes. If the genotype data are too
  big to complete a run quickly, you can limit the number of loci by
  specifying a value for the option.

Output options
--------------

    -b FILE, --summary=FILE
                        Output summary tables FILE (default='-' for standard
                        out)
    -B FILE, --bininfo=FILE
                        Output summary information about each bin to FILE
    -H N, --histomax=N  Largest bin size output in summary histogram output
                        (default=10)
    -k, --skip          Skip output of untagged or excluded loci
    -o FILE, --output=FILE
                        Output tabular LD information for bins to FILE ('-'
                        for standard out)
    -O FILE, --locusinfo=FILE
                        Output locus information to FILE
    -u FILE, --saveldpairs=FILE
                        Output pairwise LD estimates to FILE
    -x, --extra         Output inter-bin LD statistics

Notes:

-o/--output, -x/--extra and -k/--skip:

  -o/--output specifies the name of the output file containing LD
  information for the bins,-x/--extra triggers appending the inter-bin LD
  statistics to the same file, and -k/--skip skips output of the pair of
  loci if the disposition of the bin is either obligate-exclude or residual,
  or either one of the pair is in the exclude set.

-b/--summary and -H/--histomax:

  -b/--summary specifies the name of the output file containing the
  histogram table summaries for all the bins, and -H/--histomax is the
  largest bin size that is included in the table.

-B/--bininfo:

  This option specifies the name of the output file containing the summary
  information including tags, non-tags, bin size, location and spacing about
  each bin.

-O/--locusinfo:

  This option specifies the name of the output file containing the locus
  information such as location, MAF, bin number and disposition for each
  locus.

File formats
============

This section describes the file format for both input and output files.
Following are the allowable input files:

    * Genotype data in standard GLU formats
    * Pedigree file in standard GLU format
    * Pre-computed pair-wise LD values in FESTA format
    * SNP list files
    * Design score files

Following are the output files:

    * Bin info file
    * LD data file
    * Locus info file
    * Bin summary statistics file

FESTA formatted linkage disequilibrium files
--------------------------------------------

TagZilla can read in these files containing the pre-computed pair-wise LD
parameter between the SNPs in certain region.For details about the format of
these files, user can refer to this link:
http://www.sph.umich.edu/csg/qin/FESTA/sample_files/

SNP list files
--------------

These files contain lists of loci for the purpose of sub-setting, specifying
loci that must be included as tags, or excluding loci from being tags during
the analysis process.

    * '-s' is used to tell TagZilla to read in a subset of all genotyped loci
    * '-i' is used for the list of obligatorily included loci
    * '-e' is used for the list of obligatorily excluded loci

    These set of files have the same simple format, no headers, with one
    locus name on each line, and the locus name is case-sensitive. For
    example::

      rs150379
      rs469673
      rs212121
      rs210499
      rs469536

  However, if the first character of the argument on any of these options is
  a colon ':', then the remainder of the argument is processed as a
  comma-delimited list of loci.For example, -i :rs512331,rs1221. This method
  is sometimes convenient when running TagZilla iteratively from the
  command-line.

Design score files
------------------

These files contain the design score information for SNPs. Each line of the
file must contain the name of the SNP and its design score. TagZilla allows
multiple design score files to be specified from the command line, and
information in all files will considered during tag selection stage.

If the design score for a SNP is 0 or below the given threshold, that SNP
will be forced into the exclude set. If this SNP also happens to be in the
include set, then the disposition of the bin containing this SNP will be
obligate-include, the SNP will be reported as obligate_tag (because it is in
the include set) and also as one of the excluded_as_tags (because it is
forced into the exclude set). Therefore, include will take priority over
exclude in our program.

Following are some sample lines of a design score file::

  rs150379  0.8
  rs469673  0.9
  rs212121  0.7
  rs210499  0.6
  rs469536  0.5

There are four different output files, and only one of these files can be
directed to standard output, others must be output to the files with names
specified in the command line options. The output will contain the following
information about each bin chosen by TagZilla:

    * all possible tags
    * one recommend tag
    * the total number of loci in the bin
    * the summary statistics in tabular format for all the bins
    * Pair-wise LD statistics for each bin (with an option to also include the inter-bin LD statistics)
    * Locus information including the MAF, disposition and bin number for each locus.

Bin info file
-------------

The name and location of this file are specified in the '-B' option. The bin
number will appear multiple times as we output all the information for that
bin.This format is an expanded version of the output produced by the program
ldSelect (Carlson et al., 2004).The following table describes the
information produced for each bin:

=== ========================================================================
Row Description
=== ========================================================================
 1  summary line: contains the total number of sites for the bin, the number
    of tags, the number of non-tags, the number of required tags, the width,
    and the average MAF for the bin
 2  detailed location information: minimum, median, average and maximum
    location of all the loci in the bin
 3  detailed spacing information: minimum, median, average and maximum
    spacing among all the loci in the bin
 4  tag SNPs
 5  recommend tag SNP
 6  Other SNPs
 7  excluded tag SNPs
 8  bin disposition (four possible values: 'obligate-include',
    'maximal-bin', 'residual', 'obligate-exclude')
 9  Number of loci that would have been covered by the bin, note that for
    obligate include bins only the obligatory tags are considered.
=== ========================================================================

Here is an example bin info file generated by TagZilla::

  Bin 1    sites: 9, tags 3, other 6, tags required 1, width 40229, avg. MAF 49.0%
  Bin 1    Location: min 119730461, median 119769349, average 119760254, max 119770690
  Bin 1    Spacing: min 214, median 1855, average 5028, max 16495
  Bin 1    TagSnps: G11-SN-3PS10 G11-SN-3PS11 rs6204
  Bin 1    RecommendedTags: rs6204
  Bin 1    other_snps: G11-SN-3PS3 G11-SN-3PS7 rs1998182 rs2064902 rs6200 rs6686779
  Bin 1    Excluded_as_tags: rs6200
  Bin 1    Bin_disposition: maximal-bin
  Bin 1    Loci_covered: 9

  Bin 2    sites: 9, tags 5, other 4, tags required 1, width 9729, avg. MAF 48.5%
  Bin 2    Location: min 35922900, median 35926631, average 35926790, max 35932629
  Bin 2    Spacing: min 299, median 919, average 1216, max 4179
  Bin 2    TagSnps: G22-SN-E2S28 G22-SN-E2S33 rs69264 rs86582 rs9622573
  Bin 2    RecommendedTags: rs69264
  Bin 2    other_snps: rs229559 rs229566 rs6413537 rs739040
  Bin 2    Excluded_as_tags: rs6413537
  Bin 2    Bin_disposition: maximal-bin
  Bin 2    Loci_covered: 9

  Bin 3    sites: 5, tags 1, other 4, tags required 1, width 39002, avg. MAF 33.7%
  Bin 3    Location: min 119734713, median 119770024, average 119757461, max 119773715
  Bin 3    Spacing: min 939, median 2966, average 9750, max 32130
  Bin 3    TagSnps: rs1812256
  Bin 3    RecommendedTags: rs1812256
  Bin 3    other_snps: rs10754400 rs4659182 rs6667572 rs7535128
  Bin 3    Bin_disposition: maximal-bin
  Bin 3    Loci_covered: 5

  Bin 4    sites: 3, tags 3, other 0, tags required 1, width 2251, avg. MAF 32.0%
  Bin 4    Location: min 35924206, median 35924854, average 35925172, max 35926457
  Bin 4    Spacing: min 648, median 1125, average 1125, max 1603
  Bin 4    TagSnps: G22-SN-E2S24 rs1861945 rs229565
  Bin 4    RecommendedTags: rs1861945
  Bin 4    other_snps:
  Bin 4    Bin_disposition: maximal-bin
  Bin 4    Loci_covered: 3

  Bin 5    sites: 3, tags 3, other 0, tags required 1, width 4452, avg. MAF 12.8%
  Bin 5    Location: min 35926946, median 35927521, average 35928621, max 35931398
  Bin 5    Spacing: min 575, median 2226, average 2226, max 3877
  Bin 5    TagSnps: G22-SN-E1S1 rs2071710 rs229567
  Bin 5    RecommendedTags: rs229567
  Bin 5    other_snps:
  Bin 5    Bin_disposition: maximal-bin
  Bin 5    Loci_covered: 3

LD data output file
-------------------

The name and location of this file are specified in '-o' option. The first
line is the header line. The following table describes each column in the LD
data output file:

====== ======================================================================
Column Description
====== ======================================================================
   1   sequence number for identifying each bin
   2   the first locus name of the pair
   3   the second locus name of the pair
   4   the rsquared value for the pair
   5   Disposition (see the table below for details)
====== ======================================================================

All possible values for the disposition of each LD pair are summarized in
the two tables below. The first table describes different dispositions for
the tags paired with themselves, and the second table is for the rest of the
LD pairs within each bin.

Tags paired with themselves:

================ ==================================================================
Disposition      Description
================ ==================================================================
obligate-tag     An obligate tag
alternate-tag    tag in an obligate-include bin, but not the obligate tag
excluded-tag     tag for a bin that contains all obligatorily excluded loci
candidate-tag    tag for a non obligate bin with more than one possible tags
necessary-tag    tag for a bin that has only one possible tag
lonely-tag       A tag for bin with no other loci, but originally covered
                 more loci. These additional loci were removed by previous
                 iterations of the binning algorithm. This disposition is
                 primarily to distinguish these bins from singletons, which
                 intrinsically are in insufficient LD with any other locus.
singleton-tag    A tag that is not in significant LD with any other locus
                 based upon specified LD threshold.
================ ==================================================================

Note: 'recommended' will be appended to the above disposition to indicate
that it is also an optimal tag chosen among all the possible tags for a bin
by comparing the score and checking certain criteria provided that these
options are set from the command line.

Other LD pairs in the bin:

================ ==================================================================
Disposition      Description
================ ==================================================================
tag-tag          LD between tags within a bin
other-tag        LD between a non-tag and a tag
tag-other        LD between a tag and non-tag
other-other      LD between non-tags within a bin
================ ==================================================================

Note that for residual bins, the dispositions for all LD pairs within each
bin will have a 'residual' qualifier appended to them, and for obligate
exclude bins, the dispositions for all LD pairs will have an 'excluded'
qualifier appended to them.Also if the user specifies the '-x' option, the
'interbin' qualifier will appear in the disposition column for all residual
LD pairs that sit in the bottom part of this output file.The LD pairs are
formed based on each individual genotype input file, i.e., TagZilla doesn't
look for significant LD among loci in multiple input files. The LD data is
presorted by rsquared, and then locus1, and then locus2 for each bin.

Following is an example of an LD data output file::

  BIN	LNAME1	LNAME2	RSQUARED	DPRIME	DISPOSITION
  1	rs150379	rs150379	1	1	obligate-tag
  1	rs469673	rs469673	1	1	alternate-tag
  1	rs212121	rs212121	1	1	alternate-tag,recommended
  1	rs212111	rs212111	1	1	alternate-tag
  1	rs210499	rs210499	1	1	alternate-tag
  1	rs469536	rs469536	1	1	alternate-tag
  1	rs469667	rs469667	1	1	alternate-tag
  1	rs210499	rs469667	1	1	tag-tag
  1	rs150379	rs469673	1	1	tag-tag
  1	rs210534	rs150379	1	1	other-tag
  1	rs210534	rs469673	1	1	other-tag
  1	rs150379	rs212121	1	1	tag-tag
  1	rs212121	rs469673	1	1	tag-tag
  1	rs150379	rs469536	1	1	tag-tag
  1	rs210534	rs469536	1	1	other-tag
  1	rs469673	rs469536	1	1	tag-tag
  1	rs210534	rs212121	1	1	other-tag
  1	rs212121	rs469536	1	1	tag-tag
  1	rs212111	rs212121	0.932	1	tag-tag
  1	rs210499	rs212121	0.928	1	tag-tag
  1	rs212121	rs469667	0.928	1	tag-tag
  1	rs210534	rs469667	0.919	1	other-tag
  1	rs210499	rs210534	0.919	1	tag-other
  1	rs150379	rs212111	0.919	1	tag-tag
  1	rs212111	rs469673	0.919	1	tag-tag
  1	rs210499	rs150379	0.914	1	tag-tag
  1	rs210499	rs469673	0.914	1	tag-tag
  1	rs150379	rs469667	0.914	1	tag-tag
  1	rs469673	rs469667	0.914	1	tag-tag
  1	rs212111	rs469536	0.914	1	tag-tag
  1	rs210499	rs469536	0.907	1	tag-tag
  1	rs469536	rs469667	0.907	1	tag-tag
  1	rs212111	rs469667	0.865	1	tag-tag
  1	rs210499	rs212111	0.865	1	tag-tag
  1	rs210534	rs212111	0.859	0.927	other-tag
  2	rs4913553	rs4913553	1	1	candidate-tag,recommended
  2	rs1827997	rs1827997	1	1	candidate-tag
  2	rs4913553	rs1827997	0.897	1	tag-tag
  3	rs169757	rs169757	1	1	candidate-tag,recommended
  3	rs456706	rs456706	1	1	candidate-tag
  3	rs169757	rs456706	1	1	tag-tag
  4	rs240446	rs240446	1	1	singleton-tag,recommended
  5	rs10439884	rs10439884	1	1	singleton-tag,recommended
  6	rs11088417	rs11088417	1	1	singleton-tag,recommended,residual
  7	rs12172917	rs12172917	1	1	excluded-tag,recommended,excluded

Locus info data output file
---------------------------

The file name and location can be specified by the '-O' option on the
command line. The first line is the header line. The following table
describes each column in the locus info data output file:

====== ===================================================================
Column Description
====== ===================================================================
  1    Locus name
  2    Location of the locus
  3    MAF(Minor Allele Frequency)
  4    Bin number
  5    Disposition
====== ===================================================================

There are two possible disposition categories for each locus.

    * tag category, refer to the self pair-wise LD disposition table
      for details.

    * non tag category, there are two possible dispositions: 'exclude' or
      'other'.Note that for a residual bin, the disposition for all loci in
      that bin will have a 'residual' qualifier, and for an obligate exclude
      bin, the disposition for all loci in that bin will have an 'excluded'
      qualifier.

The Following is an example of a locus data info output file. The contents
in the file are sorted by bin number, and within each bin sorted by tags
first and then non tags::

  LNAME         LOCATION    MAF           BINNUM      DISPOSITION
  rs150379      9978594     0.116         1           obligate-tag
  rs469673      10011786    0.116         1           alternate-tag
  rs212121      9986010     0.136         1           alternate-tag,recommended
  rs212111      9981677     0.15          1           alternate-tag
  rs210499      9929079     0.127         1           alternate-tag
  rs469536      10016358    0.109         1           alternate-tag
  rs469667      10018800    0.127         1           alternate-tag
  rs210534      9972502     0.138         1           exclude
  rs4913553     9941912     0.375         2           candidate-tag,recommended
  rs1827997     9947160     0.35          2           candidate-tag
  rs169757      9928594     0.067         3           candidate-tag,recommended
  rs456706      10022975    0.067         3           candidate-tag
  rs240446      10000969    0.092         4           singleton-tag,recommended
  rs10439884    9993822     0.083         5           singleton-tag,recommended
  rs11088417    13262512    0.067         6           singleton-tag,recommended,residual
  rs12172917    13279705    0.417         7           excluded-tag,recommended,excluded

Bin summary statistics output file
----------------------------------

There are four types of bins: obligate-include, maximal-bin, residualand obligate-exclude.

  * Obligate-include bin is a bin with an obligatorily included locus.
  * Obligate-exclude bin is a bin with an obligatorily excluded locus chosen as a tag for the bin.
  * The rest of the bins will fall into the maximal-bin category unless the targeted number of loci or targeted number of bins have been met, in which case they will be called residual bins.

The output file includes a table summarizing the bin statistics by bin size
(a histogram) for each type of bin. The maximum size of bin shown as one row
in the table can be configured with the '-H' option, and these tables in
this output file share the common set of columns. The following table
describes the content of each column:

============ =======================================================================
Column       Description
============ =======================================================================
Bin size     number of loci in the bin
Bins         number of bins with specified bin size
%            percent of bins with specified bin size
Loci         number of loci contained in all the bins with specified bin size
%            percent of loci contained in all the bins with specified bin size
Tags         number of tags for all the bins with specified bin size
Non tags     number of non tags for all the bins with specified bin size
Avg tags     average number of tags per bin for all the bins with specified bin size
Avg width    average width for all the bins with specified bin size. The
             width for a bin is the difference between the maximum location
             and minimum location of the loci in a bin.
============ =======================================================================

At the end of the file there is a final summary table showing the bin
statistics further summarized by bin disposition. Following is an example of
the bin summary file::

  Bin statistics by bin size for obligate-include:

   bin                                     non-    avg    avg
   size  bins     %   loci      %   tags   tags   tags  width
   ----- ----- ------ ------ ------ ------ ------ ---- ------
     8       1 100.00      8 100.00      7      1  7.0  89721
   Total     1 100.00      8 100.00      7      1  7.0  89721


  Bin statistics by bin size for maximal-bin:

   bin                                     non-    avg    avg
   size  bins     %   loci      %   tags   tags   tags  width
   ----- ----- ------ ------ ------ ------ ------ ---- ------
   singl     2  50.00      2  33.33      2      0  1.0      0
     1       0   0.00      0   0.00      0      0  0.0      0
     2       2  50.00      4  66.67      4      0  2.0  49814
   Total     4 100.00      6 100.00      6      0  1.5  24907


  Bin statistics by bin size for residual:

   bin                                     non-    avg    avg
   size  bins     %   loci      %   tags   tags   tags  width
   ----- ----- ------ ------ ------ ------ ------ ---- ------
   singl     1 100.00      1 100.00      1      0  1.0      0
   Total     1 100.00      1 100.00      1      0  1.0      0


  Bin statistics by bin size for obligate-exclude:

   bin                                     non-    avg    avg
   size  bins     %   loci      %   tags   tags   tags  width
   ----- ----- ------ ------ ------ ------ ------ ---- ------
   singl     1 100.00      1 100.00      1      0  1.0      0
   Total     1 100.00      1 100.00      1      0  1.0      0


  Bin statistics by disposition:
                                                          non-    avg    avg
   disposition          bins     %   loci      %   tags   tags   tags  width
   -------------------- ----- ------ ------ ------ ------ ------ ---- ------
   obligate-include         1  14.29      8  50.00      7      1  7.0  89721
   maximal-bin              4  57.14      6  37.50      6      0  1.5  24907
   residual                 1  14.29      1   6.25      1      0  1.0      0
   obligate-exclude         1  14.29      1   6.25      1      0  1.0      0
                 Total      7 100.00     16 100.00     15      1  2.1  27050

References
==========

   1. Carlson C.S. et al. (2004) Selecting a maximally informative set of
      single-nucleotide polymorphisms for association analysis using linkage
      disequilibrium. Am. J. Hum. Genet. 74, 106-120

   2. Wigginton J.E. et al. (2005) A Note on Exact Tests of Hardy-Weinberg
      Equilibrium. Am. J. Hum. Genet. 76, 887-93

   3. Zhaohui S. Qin et al. (2006) An efficient comprehensive search
      algorithm for tagSNP selection using linkage disequilibrium criteria.
      Bioinformatics. 22(2):220-5.
