.. _user_manual-preferred_formats:

++++++++++++++++++++++++++++
Preferred Input File Formats
++++++++++++++++++++++++++++

Most GLU modules will accept any of the supported input file formats,
however some methods naturally process genotypes individually, by sample, or
by locus.  For example, :mod:`qc.dupcheck` performs operations over samples,
while :mod:`assoc.logit1` iterates over loci.  Both modules will transform
the data provided into the appropriate format, but such transformations can
be extremely memory- and time-intensive.  It may be better to use the
:mod:`transform` module to generate data in the preferred formats of the
modules you wish to use to avoid unnecessary overhead.

+----------------------------+------------------------------------------------------+
|Module                      |                     File Format                      |
+----------------------------+----------------------+----------------+--------------+
|                            |   LDAT/LBAT/HAPMAP   |     SDAT/SBAT  |   TRIP/TBAT  |
+============================+======================+================+==============+
|Data management             |                      |                |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`split`            |         x            |        x       |       x      |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`transform`        |         x            |        x       |       x      |
+----------------------------+----------------------+----------------+--------------+
|Quality control             |                      |                |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`qc.completion`    |         x            |        x       |       x      |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`qc.concordance`   |         x            |                |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`qc.dupcheck`      |                      |        x       |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`qc.summary`       |         x            |        x       |              |
+----------------------------+----------------------+----------------+--------------+
|Association tests           |                      |                |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`assoc.linear1`    |         x            |                |              |
+----------------------------+----------------------+----------------+--------------+
|  * :mod:`assoc.logit1`     |         x            |                |              |
+----------------------------+----------------------+----------------+--------------+
