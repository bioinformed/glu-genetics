===================================================================
:mod:`util.join` --- Merge rows with matching values from two files
===================================================================

.. module:: util.join
   :synopsis: Merge rows with matching values from two files

.. module:: join
   :synopsis: Merge rows with matching values from two files

Utility to merge rows with matching values from two files

Usage::

  glu util.join [options] table1 table2

Options:

  -h, --help            show this help message and exit
  -o FILE, --output=FILE
                        Output results (default is "-" for standard out)
  -1 KEY1, --key1=KEY1  Table 1 key column numbers or names.  Comma separated,
                        with default to common headers
  -2 KEY2, --key2=KEY2  Table 2 key column numbers or names
  --prefix1=P1          prefix to prepend to each non-key header name in
                        table1
  --prefix2=P2          prefix to prepend to each non-key header name in
                        table2
  -u, --unique          Require that rows with valid keys from table2 are
                        unique
  -j JOIN, --join=JOIN  Join type: 'left' (default) or 'inner'
