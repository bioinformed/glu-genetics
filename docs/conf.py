# -*- coding: utf-8 -*-
#
# GLU documentation build configuration file
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#

import os

# The default replacements for |version| and |release|.
# If 'auto', Sphinx looks for the Include/patchlevel.h file in the current Python
# source tree and replaces the values accordingly.
#
# The short X.Y version.
# version = '2.6'
version = '0.81'
# The full version, including alpha/beta/rc tags.
# release = '2.6a0'
release = '0.81'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# Set templates to the local path
templates_path=os.path.join(os.getcwd(), 'templates')

# The base URL for download links.
html_download_base_url = 'http://cgf.nci.nih.gov/glu/'

# List of files that shouldn't be included in the build.
unused_files = [
]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
last_updated_format = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
use_smartypants = True

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True
