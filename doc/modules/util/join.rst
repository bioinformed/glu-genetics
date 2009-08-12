===================================================================
:mod:`util.join` --- Merge rows with matching values from two files
===================================================================

.. module:: util.join
   :synopsis: Merge rows with matching values from two files

.. module:: join
   :synopsis: Merge rows with matching values from two files

Utility to merge rows with matching values from two files.  Like
:mod:`util.table`, the results may be optionally manipulated by creating
columns based on categorical values, filtering rows based on values, sorting
the resulting rows, and only outputting unique rows.

See the section on :ref:`user_manual-categorical` for usage of :option:`-c`,
:option:`--categorical`, :option:`--includevar`, and :option:`--excludevar`.

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
  -U, --uniquejoin      Require that rows with valid keys from table2 are
                        unique
  -j JOIN, --join=JOIN  Join type: 'left' (default) or 'inner'
  -c VAR, --categorical=VAR
                        Create indicator variables based on values of VAR
  --includevar=VARVAL  Include only records with variable VAR equal to VAL
  --excludevar=VARVAL  Exclude all records with variable VAR equal to VAL
  -s VAR, --sort=VAR    Sort rows based on values in column VAR
  -u, --uniq            Produce only unique rows by collapsing consecutive
                        duplicate rows

.. seealso::

  :mod:`util.table`
    Convert and manipulate delimited files
