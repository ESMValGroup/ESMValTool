.. _new-diagnostic:

*********************************
Making a new diagnostic or recipe
*********************************

Getting started
===============

Please discuss your idea for a new diagnostic or recipe with the development team before getting started,
to avoid disappointment later. A good way to do this is to open an
`issue on GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`_.
This is also a good way to get help.

Creating a recipe and diagnostic script(s)
==========================================
First create a recipe in esmvaltool/recipes to define the input data your analysis script needs
and optionally preprocessing and other settings. Also create a script in the esmvaltool/diag_scripts directory
and make sure it is referenced from your recipe. The easiest way to do this is probably to copy the example recipe
and diagnostic script and adjust those to your needs.

If you have no preferred programming language yet, Python 3 is highly recommended, because it is most well supported.
However, NCL, R, and Julia scripts are also supported.

Good example recipes for the different languages are:

-  python: esmvaltool/recipes/examples/recipe_python.yml
-  R: esmvaltool/recipes/examples/recipe_r.yml
-  julia: esmvaltool/recipes/examples/recipe_julia.yml
-  ncl: esmvaltool/recipes/examples/recipe_ncl.yml

Good example diagnostics are:

-  python: esmvaltool/diag_scripts/examples/diagnostic.py
-  R: esmvaltool/diag_scripts/examples/diagnostic.R
-  julia: esmvaltool/diag_scripts/examples/diagnostic.jl
-  ncl: esmvaltool/diag_scripts/examples/diagnostic.ncl

Unfortunately not much documentation is available at this stage,
so have a look at the other recipes and diagnostics for further inspiration.

Re-using existing code
======================
Always make sure your code is or can be released under a license that is compatible with the Apache 2 license.

If you have existing code in a supported scripting language, you have two options for re-using it. If it is fairly
mature and a large amount of code, the preferred way is to package and publish it on the
official package repository for that language and add it as a dependency of esmvaltool.
If it is just a few simple scripts or packaging is not possible (i.e. for NCL) you can simply copy
and paste the source code into the esmvaltool/diag_scripts directory.

If you have existing code in a compiled language like
C, C++, or Fortran that you want to re-use, the recommended way to proceed is to add Python bindings and publish
the package on PyPI so it can be installed as a Python dependency. You can then call the functions it provides
using a Python diagnostic.

Recording provenance
====================
When ESMValCore (the ``esmvaltool`` command) runs a recipe, it will first find all data and run the default preprocessor steps plus any
additional preprocessing steps defined in the recipe. Next it will run the diagnostic script defined in the recipe
and finally it will store provenance information. Provenance information is stored in the
`W3C PROV XML format <https://www.w3.org/TR/prov-xml/>`_
and also plotted in an SVG file for human inspection. In addition to provenance information, a caption is also added
to the plots.

In order to communicate with the diagnostic script, two interfaces have been defined, which are described below.
Note that for Python and NCL diagnostics much more convenient methods are available than
directly reading and writing the interface files. For other languages these are not implemented (yet).

Using the interfaces from Python
--------------------------------
Always use :meth:`esmvaltool.diag_scripts.shared.run_diagnostic` to start your script and make use of a
:class:`esmvaltool.diag_scripts.shared.ProvenanceLogger` to log provenance. Have a look at the example
Python diagnostic in esmvaltool/diag_scripts/examples/diagnostic.py for a complete example.

Using the interfaces from NCL
-----------------------------
Always call the ``log_provenance`` procedure after plotting from your NCL diag_script. You could find available shortcuts for
statistics, domain, plottype, authors and references in the ``config-references.yml`` file.

.. code-block:: console

  log_provenance(nc-file,plot_file,caption,statistics,domain,plottype,authors,references,input-files)

Have a look at the example NCL diagnostic in ``esmvaltool/diag_scripts/examples/diagnostic.ncl`` for a complete example.

Adding references
=================
Recipe and diagnostic script can include references.
When ESMValCore (the ``esmvaltool`` command) runs a recipe, it will store citation information in `BibTeX <https://en.wikipedia.org/wiki/BibTeX>`__ format.
Follow the steps below to add a reference to a recipe (or a diagnostic):

-  make a ``tag`` that is representative of the reference entry.
   For example, ``righi15gmd`` shows the last name of the first author, year and journal abbreviation.
-  add the ``tag`` to the ``references`` section in the recipe (or the diagnostic).
-  make a BibTeX file for the reference entry. There are some online tools to convert a doi to BibTeX format like https://doi2bib.org/
-  rename the file to the ``tag``.
-  add the file to the folder ``esmvaltool/references``.

Note: the ``references`` section in ``config-references.yaml`` has been replaced by the folder ``esmvaltool/references``.
