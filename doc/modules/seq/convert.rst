===================================================================
:mod:`seq.convert` --- Convert between sequence file formats
===================================================================

.. module:: seq.convert
   :synopsis: Convert between sequence file formats

.. module:: convert
   :synopsis: Convert between sequence file formats

This module converts among various sequence data file formats and can take
as input::

  * ace
  * clustal
  * fasta/fa
  * fastq/fq
  * fastq-sanger
  * fastq-solexa
  * fastq-illumina
  * genbank/gb
  * ig
  * nexus
  * phd
  * phylip
  * pir
  * stockholm
  * sff
  * sff-trim
  * swiss
  * qual

and can output::

  * clustal
  * fasta/fa
  * fastq/fq
  * fastq-sanger
  * fastq-solexa
  * fastq-illumina
  * genbank/gb
  * nexus
  * phd
  * phylip
  * stockholm
  * sff
  * qual


Usage::

  glu seq.convert [options] [input files..]

Options:

  -h, --help            show this help message and exit
  -f FORMAT, --informat=FORMAT
                        Input sequence format (default=sff).  Formats include:
                        ace, clustal, embl, fasta, fastq/fastq-sanger, fastq-
                        solexa, fastq-illumina, genbank/gb, ig
                        (IntelliGenetics), nexus, phd, phylip, pir, stockholm,
                        sff, sff-trim, swiss (SwissProt), tab (Agilent
                        eArray), qual
  -F FORMAT, --outformat=FORMAT
                        Output sequence format (default=fasta).  As above,
                        except ace, ig, pir, sff-trim, swiss.
  -o FILE, --output=FILE
                        Output file


  :mod:`seq.Newbler2SAM`
    Conver Roche/454 Newbler output to SAM/BAM format
