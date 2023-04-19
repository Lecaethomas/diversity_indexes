# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from pathlib import Path

code_dir = Path(__file__).resolve().parent.parent / 'code'
code_dir = sys.path.insert(0, str(code_dir))

# sys.path.insert(0, os.path.abspath('../code'))


# -- Project information -----------------------------------------------------

project = 'Diversity Indexes'
copyright = '2023, Thomas Lecae'
author = 'Thomas Lecae'

# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc'
    ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
                    '_build',
                    'Thumbs.db',
                    '.DS_Store'
                    ]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


master_doc= 'index'

latex_documents = [
  (master_doc, 'Diversity_indexes.tex', 'Diversity Indexes',
   'Lecaé Thomas - INRAE - ODR', 'article', False),
]

# latex_logo = '..\supports\logo.png'
latex_elements = {
    'preamble': r'''
        \usepackage{graphicx}
    
    ''',
    'maketitle': r'''
        \begin{titlepage}
            \centering
            \includegraphics[width = 60mm]{logo.png}\\[8ex]
            \vspace{1cm}
            {\Huge\bfseries Diversity Indexes \par}
            \vspace{1cm}
            {\Large User Manual \par}
            \vspace{1cm}
            {\Large Author : LECAE Thomas\\[4ex]}
            \vfill
            {\large \today\par}
        \end{titlepage}
    '''
}
