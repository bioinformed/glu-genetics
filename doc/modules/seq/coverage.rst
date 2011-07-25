============================================================================
:mod:`seq.coverage` --- Compute summary statistics on SAM/BAM files
============================================================================

.. module:: seq.coverage
   :synopsis: Compute summary statistics on SAM/BAM files

.. module:: coverage
   :synopsis: Compute summary statistics on SAM/BAM files

Usage:

  glu seq.coverage [options] in.bam

Options:

  -h, --help            show this help message and exit
  --targets=BED         Single track BED file containing all targetd intervals
  --region=REGION       Region over which to compute as "", "contig", or
                        "contig:start-stop".  Default="" (all aligned reads)
  --maxcoverage=N:M     Maximum coverage depth N to track in intervals of
                        width M
  -o FILE, --output=FILE
                        Overall coverage statistics
  -O FILE, --targetout=FILE
                        Per-target coverage statistics
