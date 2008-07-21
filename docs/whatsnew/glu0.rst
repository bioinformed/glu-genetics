*********************
What's New in GLU 0.x
*********************

Version 0.98
========================

* The following modules have been deprecated and moved to the cgfqc package:

    a. :mod:`qc.maf`, :mod:`qc.hets`, :mod:`qc.completion`, :mod:`qc.completion2`

    b. :mod:`qc.dupcheck`, :mod:`qc.dupcheck2`, :mod:`qc.concordance`

  Modules in (a) have been deprecated by :mod:`qc.summary`. Modules in (b) still
  exist in the qc package but have been consolidated, converted to use
  purely tabular output, and streamlined in potentially non-backward
  compatible ways.

* Many genotype input, output and filter command line arguments have been
  standardized.  The biggest change is genotype filters no longer support
  short argument forms.  E.g., '--excludesamples' is no longer abbreviated
  as '-x'.  e.g., see 'glu transform'.

  This change was telegraphed months ago and affects :mod:`assoc.logit1`,
  :mod:`assoc.linear1`, :mod:`qc.summary`, :mod:`transform` and :mod:`struct.pca`.  Modules in cgfqc
  are _not_ affected.

* Major updates to genedb:

  - Updated genedb to import dbSNP builds from GoldenPath, mirs from
    miRBase, and more assay panels from Affy and Illumina.  This allows us
    to create very comprehensive genedb files.

  - Lookup of genes has been expanded so that can be found by Entrez
    GeneID (e.g. 672 for BRCA1) or miRBase name (e.g., hsa-mir-311)

  - Renamed genedb.preprocess to genedb.find_regions

  - genedb.find_snps no longer needs its input piped through
    genedb.preprocess (or genedb.find_snps)


* 'glu list' now dynamically discovers modules and their descriptions

* TagZilla and its utility programs have been fully ported to use GLU
  libraries and infrastructure.  This allows tagzilla to read all
  supported GLU file formats and utilize the efficient in-memory genotype
  representations.  In addition, multiple-chromosomes and
  non-communicating regions are supported.

  These changes result in significant speedups over classical TagZilla.
  To tag all of HapMap build 22 CEU:

  ======== ============ ============= =======================
    Time     Version       Input      Main reason for speedup
  ======== ============ ============= =======================
  1h22m39s TagZilla 1.1 HapMap (text)
    49m60s GLU 0.98     HapMap (text) Improved LD calculation
    36m56s GLU 0.98     lbat (binary) Faster file parsing
  ======== ============ ============= =======================

* Hagzilla has been streamlined and moved to the tagzilla module

* Added :mod:`qc.mendel_check` to detect non-Mendelian transmission among parents
  and their offspring.

* The usual minor bug fixes and performance improvements.

Version 0.97
========================

* Major performance optimizations in many areas, with emphasis on those
  involved in genotype merging.  The net result is a 15-25x speed up when
  merging data sets with non-identical samples or loci.  Memory usage is
  also reduced dramatically (18GB -> 1GB for one specific data set).
  Notably, matrix transpose operations are 2-4x faster.

* Add support for fixed locus terms in :mod:`assoc.linear1` and :mod:`assoc.logit1`.
  Removed several older assoc scripts that have been superseded by
  :mod:`assoc.logit1`.

* :mod:`qc.summary` has been polished and proper completion logic for attempted,
  observed, and non-empty samples and loci has been added.  Attempted
  samples and loci are communicated via include options.  hwp computation
  is now optional, since it is only desired in fairly specialized
  situations and is currently very inefficient to compute.

* convert.from_lbd now normalizes chromosome 'Mt' to 'M' since Illumina
  decided it was a good idea to change conventions.

* Added options to ignore genotype or phenotype metadata to the GLU binary
  formats (:ignoreloci=y, :ignorephenos=y).  The use-case that inspired
  this change was merging data typed on the Illumina HumanHap 1M and 610Q
  chips.

Version 0.96
========================

* NEW FEATURE: Major reworking of association models, centered around the
  addition of a model formula parser (r635)::

  > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI"

  By default all terms containing genotype effects are tested, so the above
  will result in a 3df test of two genotype main effects and a single trend
  by BMI interaction term.

  To explicitly choose terms to test, in this case only the interaction::

    > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI" \
                                             --test="TREND(locus)*BMI"

  By default summary output only includes terms that are tested.  To explicitly choose terms to display::

    > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI" \
                                             --test="TREND(locus)*BMI"                 \
                                          --display="GENO(locus)+BMI+TREND(locus)*BMI"

  If a test is specified but not a model, then the model will be formed by
  taking the test terms, plus all phenotype marginal effects from the phenotype
  file::

    > glu assoc.logit1 pheno.def genos.lbat --test="GENO(locus)+TREND(locus)*BMI"

  If the phenotype file includes BMI and SMOKING, the resulting analysis
  will be the same as specifying::

    --model="GENO(locus)+TREND(locus)*BMI+BMI+SMOKING"

  To specify a model with no intercept term, e.g.::

    --model="0+GENO(locus)+EV1+EV2+EV3"

  To force an explicit intercept term (the default, anyhow), e.g.::

    --model="1+GENO(locus)+EV1+EV2+EV3"

  N.B. BACKWARD INCOMPATIBLE CHANGE AT BOTH COMMAND LINE AND API LEVELS.
  THIS INTERFACE IS STILL UNDER REVIEW AND THE FORMULA SYNTAX MAY CHANGE.

* NEW FEATURE: New module :mod:`qc.summary` that streamlines and will eventually
  superceed :mod:`qc.completion`, :mod:`qc.completion2`, :mod:`qc.hets`, :mod:`qc.maf`, :mod:`qc.hwp`. (r639)

* NEW DEV FEATURE: Fast binary iterators in C (r639)

* BUG FIX: :mod:`transform` now merges triple streams by default, in line with
  all other formats.  This can be disabled by specifying --merge=none
  (r633)

* The usual minor bug fixes and tweaks

Version 0.95
========================

* Major improvements to delimited file readers and writers, including
  better/working column range selection, output column ordering and
  filter, and support for recognizing 'csv' file extensions as csv files
  (and setting the dialect appropriately).

  Examples::

    $ glu util.table foo.csv -o foo.tsv  # Now does CSV->TSV

    # Output only the two desired columns
    $ glu assoc.logit1 phenos.def genos.lbat -o "out.txt:c=Locus,Score P-value"

* Support reading and writing Microsoft Excel files in all places the
  delimited file readers and writers are used.

  The worksheet to read or write may be specified as an augmented parameter.
  The reader accepts sheets by name or by index (strings are taken to be
  1-based, ints are taken to be 0-based, like in other fileutils).  When
  writing, the sheet name must be specified.

  Examples:

    Translate a tab to comma delimited file::

      > glu util.table foo.txt -o foo.csv

    Translate a comma delimited file to XLS::

      > glu util.table foo.csv -o foo.xls:sheet="Sheet 1"

    Translate XLS back to tab delimited::

      > glu util.table foo.xls                 -o foo.txt
      > glu util.table foo.xls:sheet=1         -o foo.txt
      > glu util.table foo.xls:sheet="Sheet 1" -o foo.txt

* The usual minor bug fixes and tweaks

Version 0.94
========================

* License and copyright have been updated as agreed upon in SAIC-Frederick
  Subcontract S07-041 modification 7, effective 2/5/2008.

* Improved handling and reporting of genotype model errors.

  When updating existing trees, this requires rebuilding the
  genoarray C module.  Usually, this is done by running::

    > python setup.py build_ext -i

  ** FAILURE TO DO SO WILL RESULT IN GLU FALLING BACK TO THE PURE PYTHON
  GENOARRAY AND RUNNING VERY VERY SLOW. DON'T SAY I DIDN'T WARN YOU. **

* Genotype loader APIs now accept all transformations as augmented
  filename parameters.  This is convenient for executing GLU commands that
  do not make all transformations available and avoids the creation of
  temporary files.  e.g.::

    > glu qc.completion data.lbat:excludeloci=naughty.lst

* Module convert.from_lbd now understands Illumina's latest manifest file
  formats (Infinium Super HD).

* All extended filename parameters for column numbers are now use 1-based
  indexing.  API level values entered as integers still assume 0-based.

  ** WARNING: BACKWARD INCOMPATIBLE CHANGE **

  NOTICE: This change has the potential to break old _shell_ scripts that
  passed columns by index.  Python API users should be immune, since it
  is silly to specify indexes as strings within Python code (right?).

* load_list, load_map, and load_table all accept c/cols/columns to specify
  columns. The previous syntax (i/index and k/key,v/value) is retained for
  backward-compatibility.

* Improved support for augmented paramters to stdin/stdout, file extension
  detection, and other corner-cases.

* Add a new module util.table to expose the flexibility of GLU's
  load_table and table_writer APIs at the commandline level.

* Rewrite of util.join to support inner and left joins, equi-join
  semantics, compound keys, and non-key header prefixes.  Docstrings and
  doctests were added.  Command-line options have changed.  This utility
  is now considered near-production grade pending only feedback from usage
  in the wild.

* Major optimizations when merging genotype matricies with identical
  columns and disjoint rows.

* Association output from :mod:`assoc.logit1` and :mod:`assoc.linear1` now includes
  genotype counts.

* Add integer 'skip' parameter to the GLU triple file reader to skip
  header row(s).

Version 0.93
========================

* New PCA analysis module for population structure in :mod:`struct.pca`.
  Functionality is rather complete, except for proper Tracy-Wisdom pvalues
  and outlier removal.

* Support for EIGENSOFT's SmartPCA genotype, locus, and subject file
  formats (via :format=smartpca or :format=eigensoft)

* New utils module for helper modules.  The first is a utility to join two
  files on a common lookup key.  This is something like the Unix 'join(1)'
  command-line application, except it uses GLU's spiffy table readers.
  This avoids the need for Excel VLOOKUPs and other ad hoc nonsense.  More
  cool utilities to come.

* Teach triple file format parsers (GLU trip and PrettyBase) how to read
  stream order, loci, and samples.  These are needed to efficiently transform
  appropriately sorted triple streams into matrix form.

* File format guessing now uses the :format=xxx to infer file formats.
  This helps, e.g., when running :mod:`transform`.  Previously one had to write::

    > glu transform -f smartpca foo.genos -o bar.xyz

  and now::

    > glu transform foo.genos:format=smartpca -o bar.xzy

* Other minor bug fixes, as usual.

Version 0.92
========================

* Support for writing STRUCTURE parameter files, since they are highly
  intertwined with the genotype data format and contents.

* Fix categorization of allele remapping in :mod:`qc.concordance`

* Major improvements to convert.from_lbd to handle new manifest
  permutations, indels, and incorrect assay alleles

* Fix completion denominator in :mod:`qc.completion2` for empty loci

* Other minor fixes, as usual

Version 0.91.2
==========================

* Allow generation of section metadata when the user environment is wonky

* Add GLU version number to standard module information output

Version 0.91.1
==========================

* Tweaks and optimizations to genotype model management

Version 0.91
========================

* Restore previous semantics for unambiguous genotypes when merging
  streams

* Fix HapMap loader to accept tri- and quadalleleic data gracefully.

Version 0.90
========================

* Preliminary phenome support added

* Added support for new file formats:
  - PLINK (ped, tped, bed, tbed) readers and writers
  - Merlin/MACH reader and writer
  - STRUCTURE writer
  - PHASE writer
  - WTCCC writer

* Performance improvements when loading SDAT-type files

* Improved error reporting

* The usual minor fixes and tweaks

Version 0.81
========================

* Fix preliminary PLINK support

* Update documentation version

Version 0.80
========================

* Added filters on minimum genotype counts, completion, and maf to
  :mod:`qc.hwp`, :mod:`qc.hets`, :mod:`qc.maf` and rename --mincount to --mingenos in
  :mod:`qc.dupcheck` (r438)

* Minor allele for monomorphic loci is now undefined (versus arbitrary)

* Fix association sample counts for non-canonical categories

* Small API changes (r432)

Version 0.69
========================

* Chromosomes are normalized from chrXX to XX when reading HapMap files.

* Add version 2 binary format that reads and writes chromosome, location,
  and strand.  Also add backward and forward compatibility metadata to all
  binary files.

* :mod:`convert.from_lbd` how to parse chromosome, location, and strand from
  Illumina manifest files.

* Many improvements and tweaks to locus metadata management.

* Refactored genotype input and output formats to a sub-package

* Bug fix for reading models with no registered alleles/genotypes from
  binary formats.

* Added new utility izip_exact to glu.lib.utils to detect errors when
  merging sequences (useful for developers)

* Add tdat file extension as a new alias for triple files
* The usual array of minor fixes and cleanups

Version 0.69
========================

* fix :mod:`assoc.linear1` output for Harvard collaborators

Version 0.61
========================

* Binary file support has been added.  All modules can read and write highly
  compressed binary files using the lbat, sbat, and tbat format names and
  file extensions.  These correspond to the text file formats ldat, sdat,
  and trip, respectively.

  The advantages of these binary files include:

    * loading genotype data is up to 20x faster
    * saving genotype data is up to 10x faster
    * genotype file size is significantly reduced, often by more than 2x over
      a compressed text version

* The in-memory representation of genotypes

* A lot of other good things behind th scenes.