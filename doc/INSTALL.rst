++++++++++++++++++
Installation Guide
++++++++++++++++++

Installation the easy way (binary installation)
===============================================

Major GLU releases will be available in the form of binary installation
packages for the following operating systems:

  ============================ ====== =====
     Operating system          bits   CPU
  ============================ ====== =====
  Microsoft Windows XP/Vista   32 bit Intel
  Microsoft Windows XP/Vista   32 bit AMD
  Fedora Linux 7               64 bit Intel
  Fedora Linux 7               64 bit AMD
  Apple OS X 10.4.6 (Tiger)    32 bit Intel
  ============================ ====== =====

Additional platforms can be added by request, although limited by my access
to suitable systems.  Contributed binary versions for additional platforms
are gratefully accepted, though only from contributors who are willing to
track and provide builds for major releases in a timely fashion.

All-in-one binary installer for Windows XP and Vista (32 bit)
-------------------------------------------------------------

1. Download zip archive for your CPU:

   * `Intel CPU`__

     __ http://code.google.com/p/glu-genetics/downloads/detail?name=glu-1.0a6-win32_Intel.zip

   * AMD CPU *(coming soon)*

2. Double-clock on the ZIP archive and extract the contents an accessible location.
   Here we use ``C:\glu-1.0a6-win32_Intel``.

3. Add GLU to your executable path (recommended)

   * Right click on your ``My Computer`` icon on your desktop or in your ``Start Menu``.

   * Select ``Properties`` | ``Advanced`` or ``System Settings`` | ``Environment Variables``

   * Under ``System variables`` select ``Path`` and click ``Edit...``

   * **Append** the directory name where you installed GLU, prefixed by a semi-colon ';'::

       ;C:\glu-1.0a6-win32_Intel

     (remember, enter the directory where GLU was extracted)

4. Create a short-cut for GLU (optional)

   * On your desktop or in your data directory, create a shortcut (right click | New | Shortcut)

   * Set the location to ``cmd`` and name to ``Command prompt for GLU``

   * Click ``Next`` and then ``Finish``

   * Right click the new ``cmd.exe`` icon

   * Select ``Properties`` and then on the Shortcut tab

   * Replace the contents of the ``Start in`` box with a default directory, if you have one

   * Click ``OK``

5. Test your settings

   * Open a command shell

     - Use the icon you created in step 4 to open a Windows command prompt
     - Or click on ``Start`` | ``Run...`` and enter ``cmd``

   * Type::

      glu

     You should see::

      Usage: glu [options] [module] [args...]

      Options:
        --version          show program's version number and exit
        -h, --help         show this help message, then exit
        -s, --stats        display program runtime statistics
        -v, --verbose      verbose error output

        Options for software developers & power users:
          --path           Display GLU package installation path
          -p, --profile    Profile GLU code to find performance bottlenecks
          --profiler=P     Set the profiler to use when -p is specified
          --gcstats        Generate statistics from the runtime object garbage
                           collector
          --gcthreshold=N  Set the threshold for triggering collection of
                           generation-0 objects

      For information on how to get started run the "intro" module,
      usually as "glu intro".  For a list of available modules run
      "glu list".

  If you see the above help message, then you are in business!


All-in-one binary archive for Linux and OS X
--------------------------------------------

1. Download binary archive from the main project site for your operating system and CPU:

   * `Linux 64 bit with an Intel CPU`__

     __ http://code.google.com/p/glu-genetics/downloads/detail?name=glu-1.0a6-Linux_Intel_EM64T.tar.gz

   * `Linux 64 bit with an AMD CPU`__

     __ http://code.google.com/p/glu-genetics/downloads/detail?name=glu-1.0a6-Linux_AMD_EM64T.tar.gz

   * `Apple OS X with an Intel CPU`__

     __ http://code.google.com/p/glu-genetics/downloads/detail?name=glu-1.0a6-OSX_Intel.tar.gz


   Below we will assume the Linux 64 bit version for Intel CPUs.

2. Extract the GLU binary archive to an an accessible location.  Here we use ``/opt/glu`` as a base directory::

     cd /opt
     mkdir glu
     tar xzf /path/to/glu-1.0a6-Linux_Intel_EM64T.tar.gz

   GLU will be extracted to a new directory named ``glu-1.0a6-Linux_Intel_EM64T``.

3. Add GLU to your executable path (recommended)

   * In **bash** or similar shells::

       export PATH=$PATH:/opt/glu/glu-1.0a6-Linux_Intel_EM64T

   * In **tcsh** or similar shells::

       setenv PATH $PATH:/opt/glu/glu-1.0a6-Linux_Intel_EM64T

   NOTE: The exact GLU directory name will depend on the version of GLU, so use
   the one appropriate for your installation.  Also, these settings are only
   temporary.  Please consult your local system administrator on how to save
   these settings into your login scripts.

4. Test your settings

   * Open a command shell

   * Type::

      glu

     You should see::

      Usage: glu [options] [module] [args...]

      Options:
        --version          show program's version number and exit
        -h, --help         show this help message, then exit
        -s, --stats        display program runtime statistics
        -v, --verbose      verbose error output

        Options for software developers & power users:
          --path           Display GLU package installation path
          -p, --profile    Profile GLU code to find performance bottlenecks
          --profiler=P     Set the profiler to use when -p is specified
          --gcstats        Generate statistics from the runtime object garbage
                           collector
          --gcthreshold=N  Set the threshold for triggering collection of
                           generation-0 objects

      For information on how to get started run the "intro" module,
      usually as "glu intro".  For a list of available modules run
      "glu list".

  If you see the above help message, then you are in business!


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

 * Python     2.7.2           from http://www.python.org/
 * setuptools 0.6c11 or newer http://pypi.python.org/pypi/setuptools
 * NumPy      1.1.0  or newer from http://numpy.scipy.org/
 * SciPy      0.9.0  or newer from http://www.scipy.org/
 * numexpr    1.4.2  or newer from http://code.google.com/p/numexpr/
 * PyTables   2.0.4  or newer from http://www.pytables.org/
 * SQLite     3.5.9  or newer from http://www.sqlite.org/
 * Ply        2.5    or newer from http://www.dabeaz.com/ply/
 * h5py       2.0    or newer from http://code.google.com/p/h5py/
 * pysam      0.5    or newer from http://code.google.com/p/pysam/
 * BioPython  1.58   or newer from http://biopython.org/wiki/Main_Page
 * matplotlib 1.0.1  or newer from http://matplotlib.sourceforge.net/
 * Cython     0.15   or newer from http://cython.org/
 * xlrd       0.7.1  or newer from http://pypi.python.org/pypi/xlrd
 * xlwt       0.7.2  or newer from http://pypi.python.org/pypi/xlwt
 * cvxopt     1.1.1  or newer from http://abel.ee.ucla.edu/cvxopt/

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

 * In addition, PyTables and h5py allows access to files based on the HDF5 standard and requires:

   - HDF5 1.8.3 or newer from http://hdfgroup.org/HDF5/

To build and install GLU, run::

        python setup.py install

Need help?
==========

Please use the GLU User Group mailing list if you have any questions,
suggestions, or would like to help improve GLU.  See the accompanying README
file for information on how to do so.
