# -*- coding: utf-8 -*-
#
# Python documentation build configuration file
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).

import sys, os, time
sys.path.append('tools/sphinxext')

# General configuration
# ---------------------

extensions = ['sphinx.addons.refcounting']

# General substitutions.
project = 'GLU'
copyright = '2008 BioInformed LLC'

# The default replacements for |version| and |release|.
#
# The short X.Y version.
version = '1.0'
# The full version, including alpha/beta/rc tags.
release = '1.0'

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# Set templates to the local path
templates_path=[os.path.join(os.getcwd(), 'templates')]

# List of files that shouldn't be included in the build.
unused_files = [
]

# Relative filename of the reference count data file.
refcount_file = 'data/refcounts.dat'

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True


# Options for HTML output
# -----------------------

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# Content template for the index page, filename relative to this file.
#html_index = 'tools/sphinxext/indexcontent.html'

# Custom sidebar templates, filenames relative to this file.
html_sidebars = {
#    'index': 'tools/sphinxext/indexsidebar.html',
}

# Additional templates that should be rendered to pages.
html_additional_pages = {
#    'download': 'tools/sphinxext/download.html',
}

# Output file base name for HTML help builder.
htmlhelp_basename = 'pydoc'


# Options for LaTeX output
# ------------------------

# The paper size ('letter' or 'a4').
latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
_stdauthor = r'Guido van Rossum\\Fred L. Drake, Jr., editor'
latex_documents = []

# Collect all HOWTOs individually

# Additional stuff for the LaTeX preamble.
latex_preamble = r'''
\authoraddress{
  \strong{GLU}\\
  Email: \email{GLU@gmail.com}
}
'''

# Documents to append as an appendix to all manuals.
latex_appendices = []
