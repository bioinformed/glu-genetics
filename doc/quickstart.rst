+++++++++++++++++
Quick Start Guide
+++++++++++++++++

For those impatient to get started, it is best to begin with a brief review
of the available modules, using the command::

  glu list

** WARNING: THIS GUIDE IS TERRIBLY OUT OF DATE AND MUST BE REPLACED WITH SOMETHING USEFUL ASAP **

and then learn about the file formats used by GLU using the command::

  glu intro.formats

Here are the steps involved in a sample quality control analysis of genotype
data, using data found in the **glu/examples/qc1 directory**: ::

  glu qc.completion samples.ldat  -o completion_report.txt

  glu qc.dupcheck   samples.ldat  --duplicates=sampleid2subjectid -o dupcheck_report.txt

  glu transform     samples.ldat  --includesamples=controls -o controls.ldat

  glu qc.summary    controls.ldat -o summary_locus.txt -O summary_sample.txt

  glu transform     samples.ldat  --renamesamples=sampleid2subjectid -o subjects.ldat \
                                  --samplemerge=sample_merge_report.txt               \
                                  --locusmerge=locus_merge_report.txt

  glu split         subjects.ldat --locusgroups=locus.map:k=SNP:v=CHROMOSOME
