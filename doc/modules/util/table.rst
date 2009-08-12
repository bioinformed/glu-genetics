===================================================================
:mod:`util.table` --- Convert and manipulate delimited files
===================================================================

.. module:: util.table
   :synopsis: Convert and manipulate delimited files

.. module:: table
   :synopsis: Convert and manipulate delimited files

Convert and manipulate delimited files, optionally creating columns based on
categorical values, filtering rows based on values, sorting the resulting
rows, and only outputting unique rows.

See the section on :ref:`user_manual-categorical` for usage of :option:`-c`,
:option:`--categorical`, :option:`--includevar`, and :option:`--excludevar`.

Usage::

  glu util.table [options] table

Options:

  -h, --help            show this help message and exit
  -c VAR, --categorical=VAR
                        Create indicator variables based on values of VAR
  --includevar=VARVAL   Include only records with variable VAR equal to VAL
  --excludevar=VARVAL   Exclude all records with variable VAR equal to VAL
  -s VAR, --sort=VAR    Sort rows based on values in column VAR
  -u, --uniq            Produce only unique rows by collapsing consecutive
                        duplicate rows
  -o FILE, --output=FILE
                        Output results (default is "-" for standard out)

.. seealso::

  :mod:`util.join`
    Merge rows with matching values from two files
