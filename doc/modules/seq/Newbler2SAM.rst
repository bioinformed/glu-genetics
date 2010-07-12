============================================================================
:mod:`seq.Newbler2SAM` --- Convert Newbler 454PairAlign.txt and the corresponding SFF files into SAM/BAM format
============================================================================

.. module:: seq.Newbler2SAM
   :synopsis: Convert Newbler 454PairAlign.txt and the corresponding SFF files into SAM/BAM format

.. module:: Newbler2SAM
   :synopsis: Convert Newbler 454PairAlign.txt and the corresponding SFF files into SAM/BAM format

Usage::

  glu seq.Newbler2SAM [options] 454PairAlign.txt[.gz|.bz2] [SFFfiles.sff..]

Options:
  -h, --help            show this help message and exit
  --reflist=FILE        Reference genome contig list
  --remapcontig=FILE    Contig remapping
  --trim=ACTION         Trim feature(s) of reads.  Comma separated list of:
                        flowkey, adapter, quality, all.  Default=all
  --maligned=ACTION     Action to perform for multiply aligned reads: keep-
                        primary, keep-all, unalign, drop.  Default=keep-all
  --mpick=METHOD        Method of selecting primary alignment when keeping
                        multiply aligned reads: best, random.  Default=best
  --unaligned=ACTION    Action to perform for unaligned reads: keep, drop.
                        Default=keep
  -o FILE, --output=FILE
                        Output SAM file
