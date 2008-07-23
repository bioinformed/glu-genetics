++++++++++++++++++
Installation Guide
++++++++++++++++++

Installation the easy way (binary installation)
===============================================

*Instructions for binary installation will be available soon...*

Major GLU releases will be available in the form of binary installation
packages for the following operating systems:

  ========================== ====== =====
     Operating system        bits   CPU
  ========================== ====== =====
  Microsoft Windows XP/Vista 32 bit Intel
  Microsoft Windows XP/Vista 32 bit AMD
  Fedora Linux 7             64 bit Intel
  Fedora Linux 7             64 bit AMD
  Apple OS X 10.4.6 (Tiger)         Intel
  ========================== ====== =====

Additional platforms can be added by request, although limited by my access
to suitable systems.  Contributed binary versions for additional platforms
are gratefully accepted, though only from contributors who are willing to
track and provide builds for major releases in a timely fashion.

Installation from source
========================

For most users, a binary installation of GLU will be the preferred method of
installation.  However, power users and software developers who wish to
extend or modify GLU will want to consider installing from source. GLU is
easily built from source on Unix-like systems, but source installation
requires several prerequisite packages, and access to Fortran and C
compilers.

Prerequisites for using GLU:

GLU requires the following free and open source packages:

 * Python   2.5.2 or newer from http://www.python.org/
 * NumPy    1.1.0 or newer from http://numpy.scipy.org/
 * SciPy    0.6.0 or newer from http://www.scipy.org/
 * PyTables 2.0.4 or newer from http://www.pytables.org/
 * SQLite   3.5.9 or newer from http://www.sqlite.org/
 * Ply      2.5   or newer from http://www.dabeaz.com/ply/

Prerequisites of the prerequisites:

 * NumPy and SciPy require a full version of BLAS and LAPACK libraries for
   many numerical functions.  As these are performance-critical, optimized
   versions have been developed, including:

   - MKL     by Intel from http://www3.intel.com/cd/software/products/asmo-na/eng/307757.htm
   - ACML    by AMD   from http://www.amd.com/acml
   - PerfLib by Sun   from http://developers.sun.com/sunstudio/overview/topics/perflib_index.html
   - ATLAS+LAPACK     from http://math-atlas.sourceforge.net/

   Otherwise, unoptimized versions are available from:

   - BLAS from http://www.netlib.org/blas/
   - LAPACK from http://www.netlib.org/lapack/

   For more information on installing NumPy and SciPy see http://www.scipy.org/Installing_SciPy

 * In addition, PyTables is based on the HDF5 standard and requires:

   - HDF5 1.8.1 or newer from http://hdfgroup.org/HDF5/

To build and install GLU, run::

        python setup.py install

To run GLU::

        glu -h
        glu <module>

For information on getting started::

        glu intro

Please contact me if you have any questions, suggestions, or would like to
help improve GLU.

| *Kevin Jacobs*
| *BioInformed LLC*
| *jacobs@bioinformed.com*
