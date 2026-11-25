.. _dev-documenting:

Documentation
=============

Kynema-SGF provides a manual written using
`Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ and comments in
the C++ source code.

Manual
------

The Kynema-SGF manual is written using a special format called
ReStructured Text (ReST) and is converted into HTML and PDF formats
using a python package Sphinx. Since the manuals are written in simple
text files, they can be version controlled alongside the source
code. Documentation is automatically generated with new updates to the
GitHub repository and deployed at `Kynema-SGF documentation site
<https://kynema.github.io/kynema-sgf>`_.

Writing documentation
`````````````````````

As mentioned previously, documentation is written using a special text format
called reStructuredText. Sphinx user manual provides a `reST Primer
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_ that
provides an overview of this format and how to write documentation using this format.


Source code documentation
-------------------------

Use concise, ordinary C++ comments to describe the purpose of functions and
other code where useful. Avoid repeating function arguments in lists.

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

To generate and spell-check the Kynema-SGF manual locally, install
`enchant` and the Sphinx packages below.

To install the required packages on Linux:

.. code-block:: console

    sudo apt-get install -y --no-install-recommends libenchant-2-dev

and on OSX with homebrew:

.. code-block:: console

    brew install enchant

To install the required python packages:

.. code-block:: console

    python3 -m pip install "sphinx<9" sphinx_rtd_theme sphinx_toolbox sphinx_copybutton pyenchant sphinxcontrib-spelling sphinxcontrib-bibtex

Build
`````

Run the following command to build the documentation:

.. code-block:: console

    sphinx-build -b spelling -W --keep-going -n docs/sphinx build-docs/spelling
    sphinx-build -b html -W --keep-going -n docs/sphinx build-docs/html

.. tip::

   On OSX, you may need to :code:`export PYENCHANT_LIBRARY_PATH=/opt/homebrew/lib/libenchant-2.dylib`

The resulting documentation is in ``build-docs/html``
directory. Documentation can also be generated in other formats,
consult `Sphinx docs
<https://www.sphinx-doc.org/en/master/usage/builders/index.html>`_ for
available formats and their usage.
