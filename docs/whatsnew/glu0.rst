******************************
What's New in GLU 0.50 -> 0.99
******************************

Version 0.98 (2008-07-15)
=========================

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
  :mod:`assoc.linear1`, :mod:`qc.summary`, :mod:`transform` and
  :mod:`struct.pca`.  Modules in cgfqc are _not_ affected.

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

Version 0.97 (2008-06-24)
=========================

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

Version 0.96 (2008-06-09)
=========================

* NEW FEATURE: Major reworking of association models, centered around the
  addition of a model formula parser (r635)::

  > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI"

  By default all terms containing genotype effects are tested, so the above
  will result in a 3df test of two genotype main effects and a single trend
  by BMI interaction term.

  To explicitly choose terms to test, in this case only the interaction::

    > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI" \
                                             --test="TREND(locus)*BMI"

  By default summary output only includes terms that are tested.  To
  explicitly choose terms to display::

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

* NEW FEATURE: New module :mod:`qc.summary` that streamlines and will
  eventually superceed :mod:`qc.completion`, :mod:`qc.completion2`,
  :mod:`qc.hets`, :mod:`qc.maf`, :mod:`qc.hwp`. (r639)

* NEW DEV FEATURE: Fast binary iterators in C (r639)

* BUG FIX: :mod:`transform` now merges triple streams by default, in line with
  all other formats.  This can be disabled by specifying --merge=none
  (r633)

* The usual minor bug fixes and tweaks

Version 0.95 (2008-05-14)
=========================

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

Version 0.94 (2008-04-28)
=========================

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

Version 0.93 (2008-04-02)
=========================

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

Version 0.92 (2008-03-03)
=========================

* Support for writing STRUCTURE parameter files, since they are highly
  intertwined with the genotype data format and contents.

* Fix categorization of allele remapping in :mod:`qc.concordance`

* Major improvements to convert.from_lbd to handle new manifest
  permutations, indels, and incorrect assay alleles

* Fix completion denominator in :mod:`qc.completion2` for empty loci

* Other minor fixes, as usual

Version 0.91.2 (2008-02-11)
===========================

* Allow generation of section metadata when the user environment is wonky

* Add GLU version number to standard module information output

Version 0.91.1 (2008-02-11)
===========================

* Tweaks and optimizations to genotype model management

Version 0.91 (2008-02-08)
=========================

* Restore previous semantics for unambiguous genotypes when merging
  streams

* Fix HapMap loader to accept tri- and quadalleleic data gracefully.

Version 0.90 (2008-02-06)
=========================

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

Version 0.81 (2007-12-20)
=========================

* Fix preliminary PLINK support

* Update documentation version

Version 0.80 (2007-12-19)
=========================

* Added filters on minimum genotype counts, completion, and maf to
  :mod:`qc.hwp`, :mod:`qc.hets`, :mod:`qc.maf` and rename --mincount to --mingenos in
  :mod:`qc.dupcheck` (r438)

* Minor allele for monomorphic loci is now undefined (versus arbitrary)

* Fix association sample counts for non-canonical categories

* Small API changes (r432)

Version 0.69 (2007-11-28)
=========================

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

Version 0.68 (2007-11-19)
=========================

* Several modules.assoc bug fixes and output improvements

* Improved modules.genedb handling of alises and official gene symbols.

* Rework of locus model and medata management:

  * Track chromosome, location, and strand per locus
  * Streams derived from HapMap files will now include locus metadata
  * Fixes to allow full use of user-specified models
  * Support added to all modules for user-specified models

* Update of modules.qc.completion as a precursor to the eventual
  (soon!) merge with completion2.

* Updates to modules.qc.{maf,hwp,hets}, including parameter renamings

* lib.fileutils.load_table now pads data rows to match the length of file
  headers

* Reference counting bug fix to _genoarray.c (r406)

* Improved profiler output to include caller summaries (r405)

* Many other minor fixes and tweaks

Version 0.67 (2007-10-23)
=========================

* fix :mod:`assoc.linear1` output for Harvard collaborators

Version 0.66 (2007-10-23)
=========================

* Missing support for allowdups in assoc.linear1

* Numerical fix for condition number bound

* Improved logit1 output

Version 0.65 (2007-10-15)
=========================

* Incremental documentation work

* Start building a locus file parser and modelmap constructor API

* Support default and non-zero design scores in TagZilla

* Add distutils package data options to allow glu/VERSION to be included.

Version 0.64 (2007-10-09)
=========================

* Optimize block matrix calls in lib.glm

* Made all GLU imports absolute in order to avoid conflicts with multi-version
  installations

* version number handling to use a common embedded version number

Version 0.63 (2007-10-05)
=========================

* Update allele renaming format

* Output and parsing improvements to hagzilla

* Allow ginfo to materialized if necessary to produce the requested output

* Add option to allow duplicate subjects in association models.  This allows
  incidence density sampling data sets once again.

* Preliminary support for reading LINKAGE format genotypes

* Added a preminary GLU coding style guide to the doc tree

* Read and write genotriples in prettybase (SeattleSNP) format

Version 0.62 (2007-09-18)
=========================

* Added parallel tagging capability to hagzilla

* Add allele1 and allele2 attributes to genolib genotype objects

* Update documentation tree to use Sphynx, which is what the Python folks
  are now using to generate their documentation.  It uses reST
  (restructured text) format.

* Initial version of ginfo script to extract metadata from gentoype files

Version 0.61 (2007-09-04)
=========================

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

* A lot of other good things behind the scenes.

Version 0.60 (2007-09-03)
=========================

* Remove support for GLU_PATH, since the built-in PYTHONPATH works nicely
  now.

* Improved genotype stream checking and debugging mode added.  Many bugs
  were detectected and fixed.

* Allow parsing of slightly mangled Illumina manifests, like those run
  through Excel, which adds extra comma delimiters in the header.

* :mod:`split` now supports reading and writing binary files

* Taught :mod:`from_lbd` how to create locus models and other minor updates
  to support the new stream APIs.

* dupcheck and dupcheck2 were updated to worki efficiently with the new
  bit-packed genotype representations, including optimized genoarray
  comparison functions for 2-bit and 4-bit homogeneous binary arrays for
  speedups of 100-1000x.

* Added repr unparser from alleles to strings.  This allows the genomatrix
  encoders to more aggressively populate model-specific string->genotype
  caches, increase the fast-path hit rates (typically from 50% -> 95%), and
  ultimately double the speed of reading large sdat files without default
  locus models.

* load_genosteam/save_genostream can now load and save binary genotype
  matrix (lbat and sbat) and binary triple files (tbat).  Performance is
  impressive

* Major polishing to binary representations

Version 0.54 (2007-08-07)
=========================

* This is the first bold step towards four major refactorings needed to
  bring in-memory and on-disk heterogeneous binary representations to GLU.

To recap, the advantages:

  * 2-bits of storage per SNP genotype in core,
    versus 8 or 64 or more required now

  * 2-bits of storage per SNP genotype on disk uncompressed,
    <2-bits of storage per SNP genotype on disk compressed,
    versus 24 bits uncompressed or 2-4 bits compressed

  * No more trinity of string/genotuple/repr

     * string representations are input and output options
     * geno tuples and reprs are now one and the same.
     * missing genotypes always evaluate to false
     * genotype objects have useful methods and attributes
     * genotypes know their model and will complain loudly if
       they used with an incompatible model

  * Due to some judicious use of C code, most operations on in core
    bit-packed genotypes are 40% faster.

  * File loading is 10x faster, saving is 5x faster.

  * Support for non-SNP genotypes, and soon non-diploid organisms

Version 0.50 (2007-06-12)
=========================

* First step in creating the "new" GLU structure: create a branches/ and trunk/ repository structure

* Major rearrangement to allow SVN to support additional tool and analysis trees

* Embed ez-setup bootstrapper to ensure we can build on systems without setuptools installed

* Move non-GLU modules out of the "core" to the tools directory
