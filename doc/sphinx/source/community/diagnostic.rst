.. _new-diagnostic:

Making a new diagnostic or recipe
*********************************

Getting started
===============

Please discuss your idea for a new diagnostic or recipe with the development team before getting started,
to avoid disappointment later. A good way to do this is to open an
`issue on GitHub <https://github.com/ESMValGroup/ESMValTool/issues>`_.
This is also a good way to get help.

.. _diagnostic_from_example:

Creating a recipe and diagnostic script(s)
==========================================
First create a recipe in esmvaltool/recipes to define the input data your analysis script needs
and optionally preprocessing and other settings.
Also create a script in the
`esmvaltool/diag_scripts <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/diag_scripts>`_
directory and make sure it is referenced from your recipe.
The easiest way to do this is probably to copy the example recipe and diagnostic
script and adjust those to your needs.

If you have no preferred programming language yet, Python 3 is highly recommended, because it is most well supported.
However, NCL, R, and Julia scripts are also supported.

Good example recipes for the different languages are:

-  python: `esmvaltool/recipes/examples/recipe_python.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_python.yml>`_
-  R: `esmvaltool/recipes/examples/recipe_r.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_r.yml>`_
-  julia: `esmvaltool/recipes/examples/recipe_julia.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_julia.yml>`_
-  ncl: `esmvaltool/recipes/examples/recipe_ncl.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_ncl.yml>`_

Good example diagnostics are:

-  python: `esmvaltool/diag_scripts/examples/diagnostic.py <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.py>`_
-  R: `esmvaltool/diag_scripts/examples/diagnostic.R <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.R>`_
-  julia: `esmvaltool/diag_scripts/examples/diagnostic.jl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.jl>`_
-  ncl: `esmvaltool/diag_scripts/examples/diagnostic.ncl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.ncl>`_

For an explanation of the recipe format, you might want to read about the
:ref:`ESMValTool recipe <esmvalcore:recipe_overview>` and have a look at the
available :ref:`preprocessor functions <esmvalcore:preprocessor>`.
For further inspiration, check out the already
:ref:`available recipes and diagnostics <recipes>`.

There is a directory
`esmvaltool/diag_scripts/shared <https://github.com/ESMValGroup/ESMValTool/tree/master/esmvaltool/diag_scripts/shared>`_
for code that is shared by many diagnostics.
This directory contains code for creating common plot types, generating output
file names, selecting input data, and other commonly needed functions.
See :ref:`api_shared` for the documentation of the shared Python code.

Re-using existing code
======================
Always make sure your code is or can be released under a license that is compatible with the Apache 2.0 license.

If you have existing code in a supported scripting language, you have two options for re-using it. If it is fairly
mature and a large amount of code, the preferred way is to package and publish it on the
official package repository for that language and add it as a dependency of ESMValTool.
If it is just a few simple scripts or packaging is not possible (i.e. for NCL) you can simply copy
and paste the source code into the ``esmvaltool/diag_scripts`` directory.

If you have existing code in a compiled language like
C, C++, or Fortran that you want to re-use, the recommended way to proceed is to add Python bindings and publish
the package on PyPI so it can be installed as a Python dependency. You can then call the functions it provides
using a Python diagnostic.

.. _recipe_documentation:

Recipe and diagnostic documentation
===================================

This section describes how to document a recipe.
For more general information on writing documentation, see :ref:`documentation`.

On readthedocs
--------------

Recipes should have a page in the :ref:`recipes` chapter which describes what
the recipe/diagnostic calculates.

When adding a completely new recipe, please start by copying
`doc/sphinx/source/recipes/recipe_template.rst.template <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/recipe_template.rst.template>`_
to a new file ``doc/sphinx/source/recipes/recipe_<name of diagnostic>.rst``
and do not forget to add your recipe to the
`index <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/index.rst>`_.

Fill all sections from the template:

- Add a brief description of the method
- Add references
- Document recipe options for the diagnostic scripts
- Fill in the list of variables required to run the recipe
- Add example images

An example image for each type of plot produced by the recipe should be added
to the documentation page to show the kind of output the recipe produces.
The '.png' files can be stored in a subdirectory specific for the recipe under
`doc/sphinx/source/recipes/figures <https://github.com/ESMValGroup/ESMValTool/blob/master/doc/sphinx/source/recipes/figures>`_
and linked from the recipe documentation page.
A resolution of 150 `dpi <https://en.wikipedia.org/wiki/Dots_per_inch>`_ is
recommended for these image files, as this is high enough for the images to look
good on the documentation webpage, but not so high that the files become large.

In the recipe
-------------
Fill in the ``documentation`` section of the recipe as described in
:ref:`esmvalcore:recipe_documentation` and add a ``description`` to each
diagnostic entry.
When reviewing a recipe, check that these entries have been filled with
descriptive content.

In the diagnostic scripts
-------------------------
Functions implementing scientific formula should contain comments with
references to the source paper(s) and formula number(s).

When reviewing diagnostic code, check that formulas are implemented according
to the referenced paper(s) and/or other resources and that the computed numbers
look as expected from literature.

.. _diagnostic_output:

Diagnostic output
=================

Typically, diagnostic scripts create plots, but any other output such as e.g.
text files or tables is also possible.
Figures should be saved in the ``plot_dir``, either in both ``.pdf`` and
``.png`` format (preferred), or
respect the ``output_file_type`` specified in the
:ref:`esmvalcore:user configuration file`.
Data should be saved in the ``work_dir``, preferably as a ``.nc``
(`NetCDF <https://www.unidata.ucar.edu/software/netcdf/>`__) file, following the
`CF-Conventions <https://cfconventions.org/>`__ as much as possible.

Have a look at the :ref:`example scripts <diagnostic_from_example>` for how to
access the value of ``work_dir``, ``plot_dir``, and ``output_file_type`` from
the diagnostic script code.
More information on the interface between ESMValCore and the diagnostic script
is available :ref:`here <esmvalcore:interface_esmvalcore_diagnostic>` and
the description of the :ref:`outputdata` may also help to understand this.

If a diagnostic script creates plots, it should save the data used to create
those plots also to a NetCDF file.
If at all possible, there will be one NetCDF file for each plot the diagnostic
script creates.
There are several reasons why it is useful to have the plotted data available
in a NetCDF file:

- for interactive visualization of the recipe on a website
- for automated regression tests, e.g. checking that the numbers are still the
  same with newer versions of libraries

If the output data is prohibitively large, diagnostics authors can choose to
implement a ``write_netcdf: false`` diagnostic script option, so writing the
NetCDF files can be disabled from the recipe.

When doing a scientific review, please check that the figures and data look as
expected from the literature and that appropriate references have been added.

.. _recording-provenance:

Recording provenance
====================

When ESMValCore (the ``esmvaltool`` command) runs a recipe,
it will first find all data and run the default preprocessor steps plus any
additional preprocessing steps defined in the recipe. Next it will run the diagnostic script defined in the recipe
and finally it will store provenance information. Provenance information is stored in the
`W3C PROV XML format <https://www.w3.org/TR/prov-xml/>`_
and provided that the provenance tree is small, also plotted in an SVG file for
human inspection.
In addition to provenance information, a caption is also added to the plots.
When contributing a diagnostic, please make sure it records the provenance,
and that no warnings related to provenance are generated when running the recipe.
To allow the ESMValCore to keep track of provenance (e.g. which input files
were used to create what plots by the diagnostic script), it needs the
:ref:`esmvalcore:interface_diagnostic_esmvalcore`.

Provenance items provided by the recipe
---------------------------------------
Provenance tags can be added in several places in the recipe.
The :ref:`esmvalcore:recipe_documentation` section provides information about
the entire recipe.

For each diagnostic in the recipe, ESMValCore supports the following additional information:

- :code:`realms` a list of high-level modeling components
- :code:`themes` a list of themes

Please see the (installed version of the) file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/config-references.yml>`_
for all available information on each item.

Provenance items provided by the diagnostic script
--------------------------------------------------
For each output file produced by the diagnostic script, ESMValCore supports the following additional information:

- :code:`ancestors` a list of input files used to create the plot.
- :code:`caption` a caption text for the plot

Note that the level of detail is limited, the only valid choices for ``ancestors`` are files produced by
:ref:`ancestor tasks<esmvalcore:ancestor-tasks>`.

It is also possible to add more information for the implemented diagnostics using the following items:

- :code:`authors` a list of authors
- :code:`references` a list of references, see :ref:`adding_references` below
- :code:`projects` a list of projects
- :code:`domains` a list of spatial coverage of the dataset
- :code:`plot_types` a list of plot types if the diagnostic created a plot, e.g. error bar
- :code:`statistics` a list of types of the statistic, e.g. anomaly

Arbitrarily named other items are also supported.

Please see the (installed version of the) file
`esmvaltool/config-references.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/config-references.yml>`_
for all available information on each item, see :ref:`esmvalcore:config-ref` for
an introduction.
In this file, the information is written in the form of ``key: value``.
Note that we add the keys to the diagnostics.
The keys will automatically be replaced by their values in the final provenance records.
For example, in the ``config-references.yml`` there is a category for types of the plots:

.. code-block:: console

  plot_types:
    errorbar: error bar plot

In the diagnostics, we add the key as:
:code:`plot_types: [errorbar]`
It is also possible to add custom provenance information by adding items to each category in this file.

In order to communicate with the diagnostic script, two interfaces have been defined,
which are described in the `ESMValCore documentation <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/interfaces.html>`_.
Note that for Python and NCL diagnostics much more convenient methods are available than
directly reading and writing the interface files. For other languages these are not implemented (yet).

Depending on your preferred programming language for developing a diagnostic,
see the instructions and examples below on how to add provenance information:

Recording provenance in a Python diagnostic script
--------------------------------------------------
Always use :meth:`esmvaltool.diag_scripts.shared.run_diagnostic` at the end of your script:

.. code-block:: python

  if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

And make use of a :class:`esmvaltool.diag_scripts.shared.ProvenanceLogger` to log provenance:

.. code-block:: python

  with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

The ``diagnostic_file`` can be obtained using :class:`esmvaltool.diag_scripts.shared.get_diagnostic_filename`.

The ``provenance_record`` is a dictionary of provenance items, for example:

.. code-block:: python

  provenance_record = {
        'ancestors': ancestor_files,
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'caption': caption,
        'domains': ['global'],
        'plot_types': ['zonal'],
        'references': [
            'acknow_project',
        ],
        'statistics': ['mean'],
      }

Have a look at the example Python diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.py <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.py>`_
for a complete example.

Recording provenance in an NCL diagnostic script
------------------------------------------------
Always call the ``log_provenance`` procedure after plotting from your NCL diag_script:

.. code-block:: console

  log_provenance(nc-file,plot_file,caption,statistics,domain,plottype,authors,references,input-files)

For example:

.. code-block:: console

  log_provenance(ncdf_outfile, \
                 map@outfile, \
                 "Mean of variable: " + var0, \
                 "mean", \
                 "global", \
                 "geo", \
                 (/"righi_mattia", "gottschaldt_klaus-dirk"/), \
                 (/"acknow_author"/), \
                 metadata_att_as_array(info0, "filename"))

Have a look at the example NCL diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.ncl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.ncl>`_
for a complete example.

Recording provenance in a Julia diagnostic script
-------------------------------------------------
The provenance information is written in a ``diagnostic_provenance.yml`` that will be located in ``run_dir``.
For example a ``provenance_record`` can be stored in a yaml file as:

.. code-block:: julia

  provenance_file = string(run_dir, "/diagnostic_provenance.yml")

  open(provenance_file, "w") do io
      JSON.print(io, provenance_records, 4)
  end

The ``provenance_records`` can be defined as a dictionary of provenance items.
For example:

.. code-block:: julia

  provenance_records = Dict()

  provenance_record = Dict(
      "ancestors" => [input_file],
      "authors" => ["vonhardenberg_jost", "arnone_enrico"],
      "caption" => "Example diagnostic in Julia",
      "domains" => ["global"],
      "projects" => ["crescendo", "c3s-magic"],
      "references" => ["zhang11wcc"],
      "statistics" => ["other"],
  )

  provenance_records[output_file] = provenance_record

Have a look at the example Julia diagnostic in
`esmvaltool/diag_scripts/examples/diagnostic.jl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.jl>`_
for a complete example.

Recording provenance in an R diagnostic script
----------------------------------------------
The provenance information is written in a ``diagnostic_provenance.yml`` that will be located in ``run_dir``.
For example a ``provenance_record`` can be stored in a yaml file as:

.. code-block:: R

  provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
  write_yaml(provenance_records, provenance_file)

The ``provenance_records`` can be defined as a list of provenance items.
For example:

.. code-block:: R

  provenance_records <- list()

  provenance_record <- list(
    ancestors = input_filenames,
    authors = list("hunter_alasdair", "perez-zanon_nuria"),
    caption = title,
    projects = list("c3s-magic"),
    statistics = list("other"),
  )

  provenance_records[[output_file]] <- provenance_record

.. _adding_references:

Adding references
=================
Recipes and diagnostic scripts can include references.
When a recipe is run, citation information is stored in `BibTeX <https://en.wikipedia.org/wiki/BibTeX>`__ format.
Follow the steps below to add a reference to a recipe (or a diagnostic):

-  make a ``tag`` that is representative of the reference entry.
   For example, ``righi15gmd`` shows the last name of the first author, year and journal abbreviation.
-  add the ``tag`` to the ``references`` section in the recipe (or the diagnostic script provenance, see recording-provenance_).
-  make a BibTeX file for the reference entry. There are some online tools to convert a doi to BibTeX format like https://doi2bib.org/
-  rename the file to the ``tag``, keep the ``.bibtex`` extension.
-  add the file to the folder ``esmvaltool/references``.

Note: the ``references`` section in ``config-references.yaml`` has been replaced by the folder ``esmvaltool/references``.

.. _testing_recipes:

Testing recipes
===============

To test a recipe, you can run it yourself on your local infrastructure or you
can ask the `@esmvalbot <https://github.com/apps/esmvalbot>`_ to run it for you.
To request a run of ``recipe_xyz.yml``, write the following comment below a pull
request:

::

   @esmvalbot Please run recipe_xyz.yml

Note that only members of the `@ESMValGroup/esmvaltool-developmentteam`_
can request runs. The memory of the `@esmvalbot`_ is limited to 16 GB and it only
has access to data available at DKRZ.

When reviewing a pull request, at the very least check that a recipes runs
without any modifications.
For a more thorough check, you might want to try out different datasets or
changing some settings if the diagnostic scripts support those.
A simple :ref:`tool <recipe_test_tool>` is available for testing recipes
with various settings.

.. _diagnostic_checklist:

Detailed checklist for reviews
==============================

This (non-exhaustive) checklist provides ideas for things to check when reviewing
pull requests for new or updated recipes and/or diagnostic scripts.

Technical reviews
-----------------

Documentation
~~~~~~~~~~~~~

Check that the scientific documentation of the new diagnostic has been added to
the user’s guide:

* A file ``doc/sphinx/source/recipes/recipe_<diagnostic>.rst`` exists
* New documentation is included in ``doc/sphinx/source/recipes/index.rst``
* Documentation follows template `doc/sphinx/source/recipes/recipe_template.rst.template`_
* Description of configuration options
* Description of variables
* Valid image files
* Resolution of image files (~150 dpi is usually enough; file size should be
  kept small)

Recipe
~~~~~~

Check yaml syntax (with ``yamllint``) and that new recipe contains:

* Documentation: description, authors, maintainer, references, projects
* Provenance tags: themes, realms

Diagnostic script
~~~~~~~~~~~~~~~~~

Check that the new diagnostic script(s) meet(s) standards.
This includes the following items:

* In-code documentation (comments, docstrings)
* Code quality (e.g. no hardcoded pathnames)
* No Codacy errors reported
* Re-use of existing functions whenever possible
* Provenance implemented

Run recipe
~~~~~~~~~~

Make sure new diagnostic(s) is working by running the ESMValTool with the recipe.

Check output of diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~

After successfully running the new recipe, check that:

* NetCDF output has been written
* Output contains (some) valid values (e.g. not only nan or zeros)
* Provenance information has been written

Check automated tests
~~~~~~~~~~~~~~~~~~~~~

Check for errors reported by automated tests

* Codacy
* CircleCI
* Documentation build

Scientific reviews
------------------

Documentation added to user’s guide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check that the scientific documentation of the new diagnostic
in ``doc/sphinx/source/recipes/recipe_<diagnostic>.rst``:

* Meets scientific documentation standard and
* Contains brief description of method
* Contains references
* Check for typos / broken text
* Documentation is complete and written in an understandable language
* References are complete

Recipe
~~~~~~

Check that new recipe contains valid:

* Documentation: description, references
* Provenance tags: themes, realms

Diagnostic script
~~~~~~~~~~~~~~~~~

Check that the new diagnostic script(s) meet(s) scientific standards.
This can include the following items:

* Clear and understandable in-code documentation including brief description of
  diagnostic
* References
* Method / equations match reference(s) given

Run recipe
~~~~~~~~~~

Make sure new diagnostic(s) is working by running the ESMValTool.

Check output of diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~

After successfully running the new recipe, check that:

* Output contains (some) valid values (e.g. not only nan or zeros)
* If applicable, check plots and compare with corresponding plots in the
  paper(s) cited


.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
