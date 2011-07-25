# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'Overview of GLU file formats'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

formats = '''\
GLU uses several tab-delimited ASCII file formats for specifying genotypes
and other data.

SUMMARY:

Genotype data:
  triples: single sample, locus, and genotype per line
  ldat   : matrix of genotypes with a locus per rows and a sample per column
  sdat   : matrix of genotypes with a sample per rows and a locus per column

Other data:
  list file: file containing a list of values to be used in a GLU operation
             in a single column
  map file : file containing two columns of values used to indicate a
             mapping from one value to another

--------------------------------------------------------------------------------

DETAILS:

TRIPLE format:

The GLU triple file format is a tab-delimited ASCII file with exactly three
columns of data and one row/line of data per genotype to be stored.  Rows
may appear in any order and no header line is allowed.

Columns:

  1. Sample name: User specified sample name with no length limit or format
                  restriction.
  2. Locus name : User specified locus name with no length limit or format
                  restriction.
  3. Genotype   : A two character genotype with alleles ' ACGT'.  Genotypes
                  containing one blank are treated as hemizygous, while
                  genotypes composed of an empty string or two blanks are
                  treated as missing.

Example:

S1	rs12345	AA
S2	rs54321	GG
S3	rs12345	AT
S1	rs54321	AG

--------------------------------------------------------------------------------

LDAT format:

The GLU LDAT format is a tab-delimited ASCII file that contains a matrix of
genotypes with a locus per row and a sample per column.  The format is as
follows:

Row    1  : 'ldat' followed by each sample name
Row    2..: The name and genotypes from each locus

Column 1  : 'ldat' in the first row, the name of each locus on each
             subsequent row
Column 2..: subject name in the first row, each of the subject's genotypes
            in each subsequent row

Genotypes are coded as a two character string with alleles ' ACGT'.
Genotypes containing one blank are treated as hemizygous, while genotypes
composed of an empty string or two blanks are treated as missing.

Example:

ldat	S1	S2	S3
rs123	AA	AT	
rs321	AG	GG	AA
rs555	CC	CC	CC

--------------------------------------------------------------------------------

SDAT format:

The GLU SDAT format is a tab-delimited ASCII file that contains a matrix of
genotypes with a sample per row and a locus per column.  The format is as
follows:

Row    1  : 'sdat' followed by each locus name
Row    2..: The name and genotypes from each sample

Column 1  : 'sdat' in the first row, the name of each sample on each
            subsequent row
Column 2..: locus name in the first row, each of the genotypes at that locus
            in each subsequent row

Genotypes are coded as a two character string with alleles ' ACGT'.
Genotypes containing one blank are treated as hemizygous, while genotypes
composed of an empty string or two blanks are treated as missing.

Example:

sdat	rs123	rs321	rs555
S1	AA	AG	CC
S2	AT	GG	CC
S3	  	AA	CC

--------------------------------------------------------------------------------

See also:

  Transform module: The transform module is used to convert between these
  formats, filter items, rename samples and loci, merge data files, and many
  other useful functions.
'''


def main():
  import pydoc
  pydoc.pager(formats)
