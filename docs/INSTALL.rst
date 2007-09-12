++++++++++++++++++
Installation Guide
++++++++++++++++++

GLU provides command-line tools for the management of large amounts of SNP
genotype data, and programs to check data quality and association between
SNP markers with continuous or discrete trait phenotypes.

Installation the easy way (binary installation)
===============================================

*Instructions for binary installation will be available soon...*

Major GLU releases will be available in the form of binary installation
packages for the following operating systems:

 ========================== ===== ===== =====================
    Operating system        bits  CPU     Installer Types
 ========================== ===== ===== =====================
 Microsoft Windows XP/Vista 32bit Intel MSI
 Microsoft Windows XP/Vista 32bit AMD   MSI
 Fedora Linux 7             64bit Intel RPM, egg
 Fedora Linux 7             64bit AMD   RPM, egg
 Apple OS X 10.4.6 (Tiger)  Intel       egg
 ========================== ===== ===== =====================

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

 * Python   2.5.1 or newer from http://www.python.org/
 * NumPy    1.0.3 or newer from http://numpy.scipy.org/
 * SciPy    0.5.3 or newer from http://www.scipy.org/
 * PyTables 2.0   or newer from http://www.pytables.org/
 * SQLite   2.4.1 or newer from http://www.sqlite.org/

Prerequisites of the prerequisites:

 * NumPy and SciPy require a full version of BLAS and LAPACK libraries for
   many numerical functions.  As these are performance-critical, optimized
   versions have been developed, including:

   - MKL  by Intel Corp from http://www3.intel.com/cd/software/products/asmo-na/eng/307757.htm
   - ACML by   AMD Corp from http://developer.amd.com/acml.jsp
   - ATLAS+LAPACK       from http://math-atlas.sourceforge.net/

 * In addition, PyTables is based on the HDF5 standard and requires:

   - HDF5 1.6.5 or newer from http://hdfgroup.org/HDF5/

Then to build and install GLU, run::

        python setup.py install

To run GLU::

        glu -h
        glu <module>

To get started::

        glu intro

Please contact me if you have any questions, suggestions, or would like to contribute better documentation.

| *Kevin Jacobs*
| *BioInformed LLC*
| *jacobs@bioinformed.com*