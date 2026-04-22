.. _dev-documenting:

Documentation
=============

Kynema-SGF comes with two different types of documentation:

- The manual, i.e., the document you are reading now, that is
  written using `Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_, and

- Inline documentation within C++ source code that are written in a format that can be
  processed automatically by `Doxygen <http://www.doxygen.nl/manual/index.html>`_

Manual
------

The Kynema-SGF manual is written using a special format called
ReStructured Text (ReST) and is converted into HTML and PDF formats
using a python package Sphinx. Since the manuals are written in simple
text files, they can be version controlled alongside the source
code. Documentation is automatically generated with new updates to the
GitHub repository and deployed at `Kynema-SGF documentation site
<https://exawind.github.io/amr-wind>`_.

Writing documentation
`````````````````````

As mentioned previously, documentation is written using a special text format
called reStructuredText. Sphinx user manual provides a `reST Primer
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_ that
provides an overview of this format and how to write documentation using this format.


Source code documentation
-------------------------

Source code (C++ files) are commented using a special format that
allows Doxygen to extract the annotated comments and create source
code documentation as well as inheritance diagrams. The
:doc:`../doxygen/html/index` documentation for the latest snapshot of
the codebase can be browsed in this manual. The `Doxygen manual
<http://www.doxygen.nl/manual/index.html>`_ provides an overview of
the syntax that must be used. Please follow the Doxygen style of
commenting code when commenting Kynema-SGF sources.

When commenting code, try to use self-documenting code, i.e., descriptive names
for variables and functions that eliminate the need to describe what is going on
in comments. In general, comments should address *why* something is being coded
in a particular way, rather than how the code does things. Try to write the code
in a clear manner so that it is obvious from reading the code instead of having
to rely on comments to follow the code structure.

Building documentation
----------------------

Dependencies
````````````

To generate the Kynema-SGF manual on a local machine, `doxygen`,
`graphviz`, `doxysphinx`, `enchant`, and `sphinx` are required.

To install the required packages on Linux:

.. code-block:: console

    sudo apt-get install -y --no-install-recommends graphviz libenchant-2-dev

and on OSX with homebrew:

.. code-block:: console

    brew install doxygen graphviz enchant

To install the required python packages:

.. code-block:: console

    pip install sphinx sphinx_rtd_theme sphinx_toolbox sphinx_copybutton pyenchant sphinxcontrib-spelling doxysphinx

Build
`````

Run the following command to build the documentation:

.. code-block:: console

    cd build && cmake -DAMR_WIND_ENABLE_DOCUMENTATION:BOOL=ON .. && cmake --build . -t docs

.. tip::

   On OSX, before running the cmake, you may need to :code:`export PYENCHANT_LIBRARY_PATH=/opt/homebrew/lib/libenchant-2.dylib`

The resulting documentation is in `docs/spinx/html`
directory. Documentation can also be generated in other formats,
consult `Sphinx docs
<https://www.sphinx-doc.org/en/master/usage/builders/index.html>`_ for
available formats and their usage.
