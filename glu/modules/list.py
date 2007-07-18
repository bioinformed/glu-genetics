# -*- coding: utf-8 -*-
'''
File:          list.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-06-13

Abstract:      List of "published" GLU modules

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

modlist = '''\
Welcome to GLU: The Genotype Library and Utilities.

The following general modules are available:

  intro            : Basic introductory information
  intro.overview   : More detailed overview of GLU
  intro.formats    : Basic information on GLU file formats
  intro.quickstart : Just the basics for how to get started
  list             : List major modules

Data management:

  transform        : Genotype data file manipulation
  split            : Split a genotype file into subsets

Genotype quality control:

  qc.completion    : Assay and sample completion
  qc.dupcheck      : Duplicate sample analysis
  qc.concordance   : Genotype concordance between data sets
  qc.hwp           : Test for deviations from Hardy-Weinberg Proportions

Genotype-phenotype association testing:

  assoc.logit1     : Single-locus association tests of dichotomous and
                     unordered polytomous outcomes with optional covariates
  assoc.linear1    : Single-locus association tests of continuous outcomes
                     with optional covariates
'''

def main():
  import pydoc
  pydoc.pager(modlist)
