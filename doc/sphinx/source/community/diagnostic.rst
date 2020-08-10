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

-  python: `esmvaltool/recipes/examples/recipe_python.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_python.yml>`_
-  R: `esmvaltool/recipes/examples/recipe_r.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_r.yml>`_
-  julia: `esmvaltool/recipes/examples/recipe_julia.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_julia.yml>`_
-  ncl: `esmvaltool/recipes/examples/recipe_ncl.yml <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_ncl.yml>`_

Good example diagnostics are:

-  python: `esmvaltool/diag_scripts/examples/diagnostic.py <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.py>`_
-  R: `esmvaltool/diag_scripts/examples/diagnostic.R <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.R>`_
-  julia: `esmvaltool/diag_scripts/examples/diagnostic.jl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.jl>`_
-  ncl: `esmvaltool/diag_scripts/examples/diagnostic.ncl <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/diag_scripts/examples/diagnostic.ncl>`_

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

Provenance items provided by the recipe
---------------------------------------
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
`ancestor tasks <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/recipe/overview.html#ancestor-tasks>`_.

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
for all available information on each item.
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
-  add the ``tag`` to the ``references`` section in the recipe (or the diagnostic).
-  make a BibTeX file for the reference entry. There are some online tools to convert a doi to BibTeX format like https://doi2bib.org/
-  rename the file to the ``tag``, keep the ``.bibtex`` extension.
-  add the file to the folder ``esmvaltool/references``.

Note: the ``references`` section in ``config-references.yaml`` has been replaced by the folder ``esmvaltool/references``.
