============================================================================
:mod:`seq.filter` --- Filter SAM/BAM files by read length and target overlap
============================================================================

.. module:: seq.filter
   :synopsis: Filter SAM/BAM files by read length and target overlap

.. module:: filter
   :synopsis: Filter SAM/BAM files by read length and target overlap

Usage::

  glu seq.filter [options] in.bam

Options:

  -h, --help            show this help message and exit
  --minreadlen=N        Minimum read length filter (default=85)
  --targets=BED         Single track BED file containing all targeted
                        intervals
  --minoverlap=N        Minimum read overlap with any target (default=20)
  --action=X            Action to perform on failing reads (filter or color,
                        default=filter)
  -o FILE, --output=FILE
                        Output BAM file
