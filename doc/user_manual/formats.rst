.. _user_manual-formats:

++++++++++++
File Formats
++++++++++++

GLU recognizes several file formats for specifying
genotypes and other data used in GLU operations.

Summary
=======

- Genotype data (ASCII text formats):

  tdat:
    tab-delimited ASCII file containing a single sample, locus, and genotype per line
  ldat:
    tab-delimited ASCII file containing a matrix of genotypes with a locus per row and a sample per column
  sdat:
    tab-delimited ASCII file containing a matrix of genotypes with a sample per row and a locus per column
  hapmap:
    tab-delimited ASCII file of HapMap data, format is same as LDAT

- Genotype data (compressed binary format):

  tbat:
    binary file containing a single sample, locus, and genotype per line
  lbat:
    binary file containing a matrix of genotypes with a locus per row and a sample per column
  sbat:
    binary file containing a matrix of genotypes with a sample per row and a locus per column

- Other data:

  list file:
    file containing a list of values in a single column
  map file:
    file containing two columns of values, used to indicate a mapping from the
    value in one column to the value in the other column
  counts file:
    file containing genotype counts (homozygous major allele, heterozygous,
    homozygous minor allele), by locus or by sample - used in the hwp, maf,
    and hets modules.


Genotype data file format details
=================================

tdat format
-----------

The GLU genotype triple (tdat) format is a tab-delimited ASCII format with
exactly three columns of data and one row (line) of data per genotype to be
stored.  Rows may appear in any order and no header line is allowed.

Columns:

  1. Sample name:
       User specified sample name with no length limit or format restriction.

  2. Locus name:
       User specified locus name with no length limit or format restriction.

  3. Genotype:
       Two character genotype with alleles ' ACGT'.  Genotypes containing one
       blank are treated as hemizygous, while genotypes composed of an empty
       string or two blanks are treated as missing.

Example::

  S1	rs12345	AA
  S2	rs54321	GG
  S3	rs12345	AT	
  S1	rs54321	AG

*Note: TBAT format is the compressed binary version of TDAT format.*


LDAT format
-----------

The GLU LDAT (locus text) format is a tab-delimited ASCII format that
contains a matrix of genotypes with a locus per row and a sample per
column:

Row 1: 'ldat' in the first column, each sample name in subsequent columns
Row 2: locus name in the first column, genotypes for that locus for each sample in the subsequent columns

Column 1: 'ldat' in the first row, each locus name in subsequent rows
Column 2: sample name in the first row, genotypes for that sample for each locus in subsequent rows

Genotypes are coded as a two character string with alleles ' ACGT'.
Genotypes containing one blank are treated as hemizygous, while genotypes
composed of an empty string or two blanks are treated as missing.

Example::

 ldat 	S1	S2	S3
 rs123	AA	AT	
 rs321	AG	GG	AA
 rs555	CC	CC	CC

*Note: LBAT format is the compressed binary version of LDAT format.*


SDAT format
-----------

The GLU SDAT (sample text) format is a tab-delimited ASCII format that
contains a matrix of genotypes with a sample per row and a locus per
column.

Row 1: 'sdat' in the first column, each locus name in subsequent columns
Row 2: sample name in the first column, genotypes for that sample for each locus in the subsequent columns

Column 1: 'sdat' in the first row, each sample name in subsequent rows
Column 2: locus name in the first row, genotypes for that locus for each sample in subsequent rows

Genotypes are coded as a two character string with alleles ' ACGT'.
Genotypes containing one blank are treated as hemizygous, while genotypes
composed of an empty string or two blanks are treated as missing.

Example::

 sdat	rs123	rs321	rs555
 S1	AA	AG	CC
 S2	AT	GG	CC
 S3	  	AA	CC


*Note: SBAT format is the compressed binary version of SDAT format.*


*See also:*

  Transform module: The transform module converts between GLU file formats,
  filters items, renames samples and loci, merges data files, and performs
  many other useful functions.
