.. _user_manual-formulae:

++++++++++++
GLU Formulae
++++++++++++

Formulae are GLU expressions that state the structural form of a model in
terms of the variables involved (similar to `S+`_, R_, SAS_, STATA_, and
other commonly used statistical packages).

..    _S+: http://www.insightful.com/
..     _R: http://www.r-project.org/
..   _SAS: http://www.sas.com/
.. _STATA: http://www.stata.com/

In GLU, these variables are derived from numerical data read in from files,
genetic effects based on a specified model, indicator variables based on
categorical variables, or arbitrary multiplicative interactions among any of
these terms.

For example, the formula::

  CANCER = BMI+SEX

indicates that the response variable, ``CANCER``, is to be modeled using an
additive model in terms of an intercept term (grand mean) and two
predictors, ``BMI`` and ``SEX``.  In this example, ``CANCER`` can be either
a categorical or continuous response variable, ``BMI``, short for *Body Mass
Index*, is usually a continuous variable and ``SEX`` is an indicator
variable coded as 0 or 1 depending if the subject is male or female.  In
your own models virtually any variable names and coding may be used, though
we shall refer back to these specific definitions in some of the examples
below.

This model is the same as::

  CANCER = 1+BMI+SEX

where ``1`` indicates that an intercept term is to be included in the model.
To exclude the intercept::

  CANCER = 0+BMI+SEX

Interactions are specified with the ``*`` operator::

  CANCER = BMI + SEX + BMI*SEX

which adds a multiplicative interaction between BMI and SEX to the model for
a total of three independent variates in the model.  Interactions do not
automatically include main effect parameters in models.  For example::

  CANCER = BMI*SEX

specifies a model with only one independent variate.

Genotype effects can be incorporated into models by including terms that
specify how genotype states are encoded as numerical parameters (see Table 1).


Table 1: Genotype encoding models
---------------------------------

+----------+----------+----+----+--------+---------------------------------------------+
| TERM     | Genotype | β1 | β2 | Effect | Description                                 |
+==========+==========+====+====+========+=============================================+
| GENO     |    AA    | 0  |  0 |   0    | Fit two independent/unconstrained genotype  |
|          +----------+----+----+--------+ effects: one for heterozygotes and one for  |
|          |    AB    | 1  |  0 |   β1   | rare homozygotes.                           |
|          +----------+----+----+--------+                                             |
|          |    BB    | 0  |  1 |   β2   |                                             |
+----------+----------+----+----+--------+---------------------------------------------+
| TREND    |    AA    | 0  |    |   0    | Fit an additive model based on the number   |
|          +----------+----+----+--------+ of rare alleles present for each genotype.  |
|          |    AB    | 1  |    |   β1   | This model is a more general form of the    |
|          +----------+----+----+--------+ model used by the Cochran-Armitage trend    |
|          |    BB    | 2  |    | 2xβ1   | test.                                       |
+----------+----------+----+----+--------+---------------------------------------------+
| DOM      |    AA    | 0  |    |   0    | Fit a dominant genetic model                |
|          +----------+----+----+--------+                                             |
|          |    AB    | 1  |    |   β1   |                                             |
|          +----------+----+----+--------+                                             |
|          |    BB    | 1  |    |   β1   |                                             |
+----------+----------+----+----+--------+---------------------------------------------+
| REC      |    AA    | 0  |    |   0    | Fit a recessive genetic model.              |
|          +----------+----+----+--------+                                             |
|          |    AB    | 0  |    |   0    |                                             |
|          +----------+----+----+--------+                                             |
|          |    BB    | 1  |    |   β1   |                                             |
+----------+----------+----+----+--------+---------------------------------------------+

where ``A`` is the reference allele and ``B`` is the other allele at the
given locus.  Unless otherwise specified, ``A`` is the more common or major
allele and ``B`` is the minor allele.

e.g. fitting two unconstrained genotype effects::

  CANCER = GENO(locus)

Covariate may be included::

  CANCER = GENO(locus) + BMI + SEX

as well as sophisticated interactions::

  CANCER = GENO(locus) + BMI + SEX + SEX*TREND(locus)

Test formulae
=============

Formulae can also be used to specify which terms in a model are to be tested
for significance.  Test formulae must include a subset of terms in the model
or else an error will be generated.

Correct::

  model: CANCER = GENO(locus) + SEX + BMI
  test:  GENO(locus) + SEX

Incorrect::

  model: CANCER = GENO(locus) + SEX*BMI
  test:  TREND(locus) + SEX

due to mismatched genotype effects and because no main-effect term for
``SEX`` is present.

**N.B.** If a test formula is specified and a model is not, then the model
will default to the terms in the test plus any available covariates.
Conversely, if a model is specified and a test is not, then the default test
includes all non-fixed genotype effect terms.

Display formulae
================

Formulae can also be used to specify terms to include in summary result
output.  Like text formulae, all terms must be present in the model or else
an error will be generated.

-----------------------------------------------------------------------

Examples taken from release notes
=================================

Simple models, like in older versions of assoc.logit1 and assoc.linear1, can
be fit with:

 > glu assoc.logit1 pheno.def genos.lbat --test="GENO(locus)"

The "locus" is the "scan" variable or placeholder for the current SNP.  GENO
adds the usual two genotype effect terms, TREND, DOM, and REC are also
supported for trend, dominant, and recessive models.  In addition, terms for
missing or non-missing genotype status can be included with MISSING(locus),
NOTMISSING(locus).

The "--test" parameter acts like the --genomodel in previous versions.  By
default, all covariates in the phenotype file are used as covariates.  In
contrast:

> glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)"

specifies a model with an intercept term and only two genotype effect terms.
In general when no model is specified, the default model includes an
intercept term, all terms listed in "--test" and all covariates in the
phenotype file.  In all cases, test must include a subset of the terms that
appear in the model.

Models (and tests) can contain more than single terms.  These can be
specified by more complex formulae that include arbitrary genotype effect
and covariate interactions:

> glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI"

By default, the summary output file has columns for parameter estimates for
all terms in test.  This can be altered by specifying a "--display"
parameter, which may list any subset of terms from the model.

Major reworking of association models, centered around the addition of a
model formula parser (r635)::

> glu assoc.logit1 pheno.def genos.lbat
> --model="GENO(locus)+BMI+TREND(locus)*BMI"

By default all terms containing genotype effects are tested, so the above
will result in a 3df test of two genotype main effects and a single trend by
BMI interaction term.

To explicitly choose terms to test, in this case only the interaction::

  > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI"
                                           \ --test="TREND(locus)*BMI"

By default summary output only includes terms that are tested.  To
explicitly choose terms to display::

  > glu assoc.logit1 pheno.def genos.lbat --model="GENO(locus)+BMI+TREND(locus)*BMI" \
                                           --test="TREND(locus)*BMI"                 \
                                        --display="GENO(locus)+BMI+TREND(locus)*BMI"

If a test is specified but not a model, then the model will be formed by
taking the test terms, plus all phenotype marginal effects from the
phenotype file::

  > glu assoc.logit1 pheno.def genos.lbat --test="GENO(locus)+TREND(locus)*BMI"

If the phenotype file includes BMI and SMOKING, the resulting analysis will
be the same as specifying::

  --model="GENO(locus)+TREND(locus)*BMI+BMI+SMOKING"

To specify a model with no intercept term, e.g.::

  --model="0+GENO(locus)+EV1+EV2+EV3"

To force an explicit intercept term (the default, anyhow), e.g.::

  --model="1+GENO(locus)+EV1+EV2+EV3"
