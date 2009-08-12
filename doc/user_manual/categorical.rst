.. _user_manual-categorical:

+++++++++++++++++++++
Categorical Variables
+++++++++++++++++++++

Many datasets feature categorical variables and GLU provides tools to filter
and aid in forming models using these variables.  A categorical variable
(sometimes called a nominal variable) is one that has two or more categories
or values, but there is no intrinsic ordering to the categories.  For
example, gender is a categorical variable having at least two categories
(male and female) and there is no intrinsic ordering to the categories.
Hair color is also a categorical variable having a number of categories
(blond, brown, brunette, red, etc.) and again, there is no agreed way to
order these from highest to lowest.

GLU provides methods to filter data by inclusion or exclusion of specific
category values.  In addition, GLU allows categorical variables to be
encoded as one or more binary (0/1) variables, which is useful when fitting
association models (as in :mod:`assoc.logit1` or :mod:`assoc.linear1`).  It
should be noted that genotype effects, although categorical, are encoded via
a different mechanism :ref:`user_manual-formulae`.

The following options are available in :mod:`util.table`, :mod:`util.join`,
:mod:`assoc.logit1`, and :mod:`assoc.linear1`:

  -c VAR, --categorical=VAR  Create indicator variables based on values of VAR
  --includevar=VARVAL        Include only records with variable VAR equal to VAL
  --excludevar=VARVAL        Exclude all records with variable VAR equal to VAL

The `-c` or `--categorical` option also accepts several additional
parameters that may be appended to VAR:

+------------------+------------------------------------------------------------+
| Parameter        | Description                                                |
+==================+============================================================+
| :prefix=VALUE    | Prefix of indicator variable names.  Defaults to           |
|                  | <VAR>_, where VAR is the column heading in the input       |
|                  | file.                                                      |
+------------------+------------------------------------------------------------+
| :ref=VALUES      | Comma separated list of referent value(s).  No             |
|                  | indicators will be added for these values and all values   |
|                  | will be set to 0 (or the :no value).  This option is most  |
|                  | useful when constructing variables for  association tests. |
+------------------+------------------------------------------------------------+
| :include=VALUES  | Create indicators for only the listed values (comma        |
|                  | separated).  All other values will be treated as missing.  |
+------------------+------------------------------------------------------------+
| :exclude=VALUES  | Do not create indicators for the listed values (comma      |
|                  | separated) and treat these values as missing.              |
+------------------+------------------------------------------------------------+
| :yes=VALUE       | Output value for fields that belong to a given category.   |
|                  | Default value is "1".                                      |
+------------------+------------------------------------------------------------+
| :no=VALUE        | Output value for fields that do not belong to a given      |
|                  | category, default value is "0".                            |
+------------------+------------------------------------------------------------+
| :missing=VALUE   | Output value for categories with missing or excluded       |
|                  | values.  Defaults to blank.                                |
+------------------+------------------------------------------------------------+

Example Usage
-------------

For the follow examples, assume the following tab-delimited input file,
named example.txt:

===== ======= =======
ID    SEX     REGION
===== ======= =======
001   MALE    ARIZONA
002   MALE    ALASKA
003   FEMALE
004   FEMALE  ALASKA
005   MALE    ELBONIA
006   FEMALE  ARIZONA
007   MALE    ALASKA
008   ?
009   MALE    ELBONIA
010           ALASKA
===== ======= =======

To create indicator (0/1) variables for the SEX column, the command::

  > glu util.table example.txt -c SEX

or::

  > glu util.table example.txt --categorical=SEX

which result in:

===== ======= ========= ====== ========== ========
ID    SEX     REGION    SEX\_? SEX_FEMALE SEX_MALE
===== ======= ========= ====== ========== ========
001   MALE    ARIZONA     0       0         1
002   MALE    ALASKA      0       0         1
003   FEMALE              0       1         0
004   FEMALE  ALASKA      0       1         0
005   MALE    ELBONIA     0       0         1
006   FEMALE  ARIZONA     0       1         0
007   MALE    ALASKA      0       0         1
008   ?                   1       0         0
009   MALE    ELBONIA     0       0         1
010           ALASKA
===== ======= ========= ====== ========== ========

which now includes variables for for each non-missing value in the input
file SEX column.  The default prefix (SEX\_) is used based on the name of the
categorical variable.

To avoid the creation of an indicator variable for "?", the command::

  > glu util.table example.txt -c SEX:exclude=?

produces:
	
===== ======= ========= ========== ========
ID    SEX     REGION    SEX_FEMALE SEX_MALE
===== ======= ========= ========== ========
001   MALE    ARIZONA       0         1
002   MALE    ALASKA        0         1
003   FEMALE                1         0
004   FEMALE  ALASKA        1         0
005   MALE    ELBONIA       0         1
006   FEMALE  ARIZONA       1         0
007   MALE    ALASKA        0         1
008   ?                     0         0
009   MALE    ELBONIA       0         1
010           ALASKA
===== ======= ========= ========== ========

In the previous examples, missing values were set to blank.  An alternate
value may be specified::

  > glu util.table example.txt -c SEX:exclude=?:missing=MISSING

===== ======= ========= ========== ========
ID    SEX     REGION    SEX_FEMALE SEX_MALE
===== ======= ========= ========== ========
001   MALE    ARIZONA       0         1
002   MALE    ALASKA        0         1
003   FEMALE                1         0
004   FEMALE  ALASKA        1         0
005   MALE    ELBONIA       0         1
006   FEMALE  ARIZONA       1         0
007   MALE    ALASKA        0         1
008   ?                     0         0
009   MALE    ELBONIA       0         1
010           ALASKA     MISSING   MISSING
===== ======= ========= ========== ========

These options can be combined to generate sophisticated results.  For
example, we can exclude SEX=? as a missing value and create indicator
variables for regions::

  > glu util.table example.txt -c SEX:exclude=? -c REGION:ref=ELBONIA:exclude=ALASKA:prefix=

ELBONIA is taken as the reference region, thus no column is created for it
and rows with the value ELBONIA are coded with all "0" indicators.  Also,
the ALASKA region is excluded and rows with the value ALASKA are coded with
blank indicators:

=== ====== ======= ========== ======== =======
ID  SEX    REGION  SEX_FEMALE SEX_MALE ARIZONA
=== ====== ======= ========== ======== =======
001 MALE   ARIZONA     0         1        1
002 MALE   ALASKA      0         1
003 FEMALE             1         0
004 FEMALE ALASKA      1         0
005 MALE   ELBONIA     0         1        0
006 FEMALE ARIZONA     1         0        1
007 MALE   ALASKA      0         1
008 ?
009 MALE   ELBONIA     0         1        0
010        ALASKA
=== ====== ======= ========== ======== =======

Rows can be filtered by the value of categorical variables using the
--includevar and --excludevar options.

For example::

  > glu util.table example.txt --includevar=SEX=FEMALE --includevar=REGION

returns:

=== ====== =======
ID  SEX    REGION
=== ====== =======
004 FEMALE ALASKA
006 FEMALE ARIZONA
=== ====== =======
