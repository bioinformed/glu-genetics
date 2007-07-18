# -*- coding: utf-8 -*-
'''
File:          quickstart.py

Authors:       Kevin Jacobs (jacobske@bioinformed.com)

Created:       2007-06-13

Abstract:      Quickstart help of GLU

Requires:      Python 2.5, glu

Revision:      $Id$
'''

__copyright__ = 'Copyright (c) 2007 Science Applications International Corporation ("SAIC")'
__license__   = 'See GLU license for terms by running: glu license'

quickstart = '''\
Welcome to GLU: The Genotype Library and Utilities.

For those impatient to get started, it is best to begin with a brief review
of the available modules, using the command:

  glu list

and then learning about the data file formats used by GLU wuth the command:

  glu intro.formats

Here is a sample quality control analysis of genotype data using data found
in the glu/examples/qc1 directory:

  glu qc.completion samples.ldat  -o completion_report.txt
  glu qc.dupcheck   samples.ldat  --duplicates=sampleid2subjectid -o dupcheck_report.txt
  glu transform     samples.ldat  --includesamples=controls -o controls.ldat
  glu qc.hwp        controls.ldat -o hwp_report.txt
  glu transform     samples.ldat  --renamesamples=sampleid2subjectid -o subjects.ldat \\
                                  --samplemerge=sample_merge_report.txt               \\
                                  --locusmerge=locus_merge_report.txt
  glu split         subjects.ldat --locusgroups=locus.map:k=SNP:v=CHROMOSOME
'''


def main():
  import pydoc
  pydoc.pager(quickstart)
