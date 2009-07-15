==========================================================
:mod:`ld.surrogates` --- Find LD surrogates for SNPs
==========================================================

.. module:: ld.surrogates
   :synopsis: Find LD surrogates for SNPs

.. module:: surrogates
   :synopsis: Find LD surrogates for SNPs

Introduction
============

This user's guide explains the use of the :mod:`ld.surrogates`, a part of
the :mod:`ld.tagzilla` module.  The options for this module are very similar
to that of TagZilla, except that the best LD surrogate is detected for all
excludes (-e) and design entries (-D) that do not meet the required minimum
threshold.  As no binning is performed, obligate includes are not accepted.
