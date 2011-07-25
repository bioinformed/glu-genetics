# -*- coding: utf-8 -*-

__gluindex__  = True
__abstract__  = 'New user quick start guide'
__copyright__ = 'Copyright (c) 2007-2009, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

quickstart = '''\
Welcome to GLU: The Genotype Library and Utilities.

For those impatient to get started, it is best to begin with a brief review
of the available modules, using the command:

  glu list

and then learning about the data file formats used by GLU with the command:

  glu intro.formats

Here is a sample quality control analysis of genotype data using data found
in the glu/examples/qc1 directory:

  glu qc.summary    samples.ldat  -o locus_summary.txt -O sample_summary.txt
  glu qc.dupcheck   samples.ldat  --duplicates=sampleid2subjectid -o dupcheck_report.txt
  glu transform     samples.ldat  --includesamples=controls -o controls.ldat
  glu transform     samples.ldat  --renamesamples=sampleid2subjectid -o subjects.ldat \\
                                  --samplemerge=sample_merge_report.txt               \\
                                  --locusmerge=locus_merge_report.txt
  glu split         subjects.ldat --locusgroups=locus.map:k=SNP:v=CHROMOSOME
'''


def main():
  import pydoc
  pydoc.pager(quickstart)
