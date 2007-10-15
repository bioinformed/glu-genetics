==========================================================
:mod:`qc.completion` --- Estimate genotype completion rate
==========================================================

.. module:: qc.completion
   :synopsis: Estimate genotype completion rate

.. module:: completion
   :synopsis: Estimate genotype completion rate

The completion module looks at genotype data by locus and by sample, and
determines the percentage of missing values. The user can indicate the
number of rows and columns previously dropped from the dataset, so that
overall completion can be computed accurately.
