*********************
What's New in GLU 1.0
*********************

Version 1.0a6 (2009-01-06)
==========================

* Major re-write of genotype model encoding.  This corrects a major design
  flaw which caused excessive amounts of memory to be used to process
  monomorphic SNPs or other instances of incomplete genotype models.  The
  details are fairly low-level and technical, but the net result is that GLU
  is much smarter about allocating new model objects, performs faster for many
  operations, and requires less memory.

  Although known in principle, this issue was first reported in the wild
  when Jun Lu was utilizing HapMap build 23, which includes 125k monomorphic
  SNPs (incomplete models).  Over 4.7 GB of RAM and 2m22s were needed to
  subset the data using GLU 1.0a5 with the old model management strategy,
  but now only 315 MB of RAM and 5.7s are needed to perform the same
  operations.  A pleasant side-effect is that runtime performance is greatly
  improved for this and many other operations.  This 15x reduction in the
  amount of memory and a 25x reduction in time required is a substantive
  start on optimizing GLU for operation on more modest desktop hardware,
  though clearly more work is needed.

  Special thanks to Jun Lu for his help in testing this fairly significant
  set of changes.

* GLU's genotype file format support is now fully "pluggable", in that new
  formats can be added by placing code in a plug-ins directory and will be
  automatically made available to all programs.  The API is not yet
  documented, but this feature removes a major barrier for adding custom and
  user-defined file formats to GLU.  This feature also fixes a number of
  internal limitations and bugs.

  Some formats can no longer be specified by file extension.  e.g.,::

    glu transform mydata.lbat -o mydata.structure

  is now invalid.  Really, there are no files with the .structure extension in
  the wild (nor should there be).  What any sensible user wants is::

    glu transform mydata.lbat -o mydata.dat:format=structure

  or::

    glu transform mydata.lbat -F structure -o mydata.dat

* :mod:`tagzilla`'s founder filter is now based on an exclude filter, since phenome
  information may not include all individuals.  This ensures that those with
  unknown descent are assumed to be founders, rather than non-founders.

* Enhance association testing output to include standard errors
  (logit1/linear1), genotype counts by category (logit1), maf by category
  (logit1), and align degenerate categories (logit1).

* Renaming alleles and recoding models is now done before applying sample
  and locus renaming.  The original behavior had identifiability issues.

* Added support for Illumina's genotype matrix format by adding a new
  genotype representation (missing genotype='--').  To use, specify
  representation "-g isnp".

* Enabled genotype filter command line parameters in ginfo.

* Standardize command-line help output to use the standard error output stream.

* Documentation updates based on contributions from Dennis Maeder, Jun
  Lu, Zhaoming Wang and Dan Eisenberg.

Version 1.0a5 (2008-10-01)
==========================

* Added support for raw WTCCC genotype files

* Fix bugs in locus and pedigree file readers

* Performance optimizations to binary file reader

* Add support for dbSNP 129 to genedb

* Support new Illumina manifest format strangeness

* Renamed --outfile options to --output in the following modules:

  * :mod:`genedb.annotate`
  * :mod:`genedb.find_snps`
  * :mod:`genedb.find_regions`
  * :mod:`tagzilla.coverage`
  * :mod:`tagzilla.surrogates`
  * :mod:`tagzilla.ldmatrix`
  * :mod:`tagzilla.tagzilla`

* Rename fileutils internal APIs, with backwardly compatible aliases

* Add checks to prevent pathological weights in GLM iterations.  This
  prevents infinite loops in the LAPACK fitting code for extremely sparse or
  ill-conditioned data.  Improve :mod:`assoc.logit1` and
  :mod:`assoc.linear1` to be robust to these new failure conditions.

* Add non-founder filter to :mod:`tagzilla` and related modules.

Version 1.0a4 (2008-08-11)
==========================

* Fix multiple programs that broke due to overhasty standardization just
  before 1.0a3 was released.

* Correct :mod:`tagzilla` -u/--saveldpairs not respecting region boundaries
  and metadata lifetimes.  Reported by Nick Orr.

* Specification of genotypes is now optional in :mod:`assoc.logit1` and
  :mod:`assoc.linear1` to to allow for fitting pure null models

* Fix output routines in :mod:`qc.summary` to work when no valid samples or loci
  are observed.

* Minor internal tweaks and documentation improvements

Version 1.0a3 (2008-07-31)
==========================

* Modified standard options to use -f/--informat and -g/--ingenorepr in all
  cases to be more consistent.  Similarly, moved --renameloci and
  --renamesamples to the transformation section, as they are somewhat out of
  place in the filter section.

* Add support for transparent bzip2 (.bz2) stream compression and decompression

* Fix to logistic regression due to a change in Numpy 1.1

* Added concordance rate to :mod:`qc.dupcheck` output and an option to check only
  expected duplicates.

* Several documentation updates

Version 1.0a2 (2008-07-27)
==========================

* Rate parameters for :mod:`tagzilla` and :mod:`qc.dupcheck` now take
  decimal rates and not integer percentages.

* Fixed a missing import that prevented :mod:`qc.dupcheck` from running.

* Corrected a metadata sequencing bug when recoding or merging genotriple
  streams.

* Corrected a bug in code that selects the optimal genotype merge
  algorithm that affected merging genotriple files (tdat/tbat/PrettyBase).

Version 1.0a1 (2008-07-23)
==========================

* Update version numbers and tag release

* :mod:`assoc.logit1` and :mod:`assoc.linear1` are smarter about dropping
  records with missing data.  Only columns used in the model are checked for
  missing values, which allows use of phenotype files with many more
  variables than will be used in a given analysis.  In addition, the subject
  ID and phenotype columns are now configurable.

* Refactored genedb and related code to search for database files based on
  an optional database name and search path. If not specified, a series of
  standard database names and paths will be explored.

  The following modules no longer take the database name as the first argument:

    * :mod:`genedb.find_snps`
    * :mod:`genedb.find_regions`
    * :mod:`genedb.annotate`

  Instead, a '-g/--genedb' option is provided.  E.g.::

    > glu genedb.annotate -g genome36.3 assoc.txt -o assoc_annotated.txt

  This will look for the genome36.3.db file in the standard GLU genedb paths
  (places like /usr/local/share/genedb/).  Absolute paths are also allowed::

    > glu genedb.annotate -g /path/to/genome36.3.db assoc.txt -o assoc_annotated.txt

* Many documentation improvements

* Minor bug fixes, including an internal issue with the genotype counts in
  :mod:`qc.summary` (r725) and to the PLINK genotype writers (r724,r741).
