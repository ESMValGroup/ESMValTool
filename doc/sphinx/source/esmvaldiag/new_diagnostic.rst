.. _new-diagnostic:

***************************************
Contributing a new diagnostic or recipe
***************************************

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
A good example recipe is esmvaltool/recipes/examples/recipe_python.yml
and a good example diagnostic is esmvaltool/diag_scripts/examples/diagnostic.py.

If you have no preferred programming language yet, Python 3 is highly recommended, because it is most well supported.
However, NCL, R, and Julia scripts are also supported.

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

Interfaces and provenance
=========================
When ESMValTool runs a recipe, it will first find all data and run the default preprocessor steps plus any
additional preprocessing steps defined in the recipe. Next it will run the diagnostic script defined in the recipe
and finally it will store provenance information. Provenance information is stored in the
`W3C PROV XML format <https://www.w3.org/TR/prov-xml/>`_
and also plotted in an SVG file for human inspection. In addition to provenance information, a caption is also added
to the plots.

In order to communicate with the diagnostic script, two interfaces have been defined, which are described below.
Note that for Python and NCL diagnostics much more convenient methods are available than
directly reading and writing the interface files. For other languages these are not implemented yet.

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

Generic interface between backend and diagnostic
------------------------------------------------
To provide the diagnostic script with the information it needs to run (e.g. location of input data, various settings),
the backend creates a YAML file called settings.yml and provides the path to this file as the first command line
argument to the diagnostic script.

The most interesting settings provided in this file are

.. code-block:: yaml

  run_dir:  /path/to/recipe_output/run/diagnostic_name/script_name
  work_dir: /path/to/recipe_output/work/diagnostic_name/script_name
  plot_dir: /path/to/recipe_output/plots/diagnostic_name/script_name
  input_files:
    - /path/to/recipe_output/preproc/diagnostic_name/ta/metadata.yml
    - /path/to/recipe_output/preproc/diagnostic_name/pr/metadata.yml

Custom settings in the script section of the recipe will also be made available in this file.

There are three directories defined:

- :code:`run_dir` use this for storing temporary files
- :code:`work_dir` use this for storing NetCDF files containing the data used to make a plot
- :code:`plot_dir` use this for storing plots

Finally :code:`input_files` is a list of YAML files, containing a description of the preprocessed data. Each entry in these
YAML files is a path to a preprocessed file in NetCDF format, with a list of various attributes.
An example preprocessor metadata.yml file could look like this:

.. code-block:: yaml

  ? /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
  : alias: GFDL-ESM2G
    cmor_table: CMIP5
    dataset: GFDL-ESM2G
    diagnostic: diagnostic_name
    end_year: 2002
    ensemble: r1i1p1
    exp: historical
    filename: /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
    frequency: mon
    institute: [NOAA-GFDL]
    long_name: Precipitation
    mip: Amon
    modeling_realm: [atmos]
    preprocessor: preprocessor_name
    project: CMIP5
    recipe_dataset_index: 1
    reference_dataset: MPI-ESM-LR
    short_name: pr
    standard_name: precipitation_flux
    start_year: 2000
    units: kg m-2 s-1
    variable_group: pr
  ? /path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
  : alias: MPI-ESM-LR
    cmor_table: CMIP5
    dataset: MPI-ESM-LR
    diagnostic: diagnostic_name
    end_year: 2002
    ensemble: r1i1p1
    exp: historical
    filename: /path/to/recipe_output/preproc/diagnostic1/pr/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc
    frequency: mon
    institute: [MPI-M]
    long_name: Precipitation
    mip: Amon
    modeling_realm: [atmos]
    preprocessor: preprocessor_name
    project: CMIP5
    recipe_dataset_index: 2
    reference_dataset: MPI-ESM-LR
    short_name: pr
    standard_name: precipitation_flux
    start_year: 2000
    units: kg m-2 s-1
    variable_group: pr


Generic interface between diagnostic and backend
------------------------------------------------

After the diagnostic script has finished running, the backend will try to store provenance information. In order to
link the produced files to input data, the diagnostic script needs to store a YAML file called :code:`diagnostic_provenance.yml`
in it's :code:`run_dir`.

For output file produced by the diagnostic script, there should be an entry in the :code:`diagnostic_provenance.yml` file.
The name of each entry should be the path to the output file.
Each file entry should at least contain the following items

- :code:`ancestors` a list of input files used to create the plot
- :code:`caption` a caption text for the plot
- :code:`plot_file` if the diagnostic also created a plot file, e.g. in .png format.

Each file entry can also contain items from the categories defined in the file :code:`esmvaltool/config_references.yml`.
The short entries will automatically be replaced by their longer equivalent in the final provenance records.
It is possible to add custom provenance information by adding custom items to entries.

An example :code:`diagnostic_provenance.yml` file could look like this

.. code-block:: yaml

  ? /path/to/recipe_output/work/diagnostic_name/script_name/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002_mean.nc
  : ancestors:[/path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_GFDL-ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc]
    authors: [andela_bouwe, righi_mattia]
    caption: Average Precipitation between 2000 and 2002 according to GFDL-ESM2G.
    domains: [global]
    plot_file: /path/to/recipe_output/plots/diagnostic_name/script_name/CMIP5_GFDL ESM2G_Amon_historical_r1i1p1_T2Ms_pr_2000-2002_mean.png
    plot_type: zonal
    references: [acknow_project]
    statistics: [mean]

  ? /path/to/recipe_output/work/diagnostic_name/script_name/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002_mean.nc
  : ancestors: [/path/to/recipe_output/preproc/diagnostic_name/pr/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002.nc]
    authors: [andela_bouwe, righi_mattia]
    caption: Average Precipitation between 2000 and 2002 according to MPI-ESM-LR.
    domains: [global]
    plot_file: /path/to/recipe_output/plots/diagnostic_name/script_name/CMIP5_MPI-ESM-LR_Amon_historical_r1i1p1_T2Ms_pr_2000-2002_mean.png
    plot_type: zonal
    references: [acknow_project]
    statistics: [mean]

You can check whether your diagnostic script successfully provided the provenance information to the backend by
verifying that

- for each output file in the :code:`work_dir`, a file with the same name, but ending with _provenance.xml is created
- any NetCDF files created by your diagnostic script contain a 'provenance' global attribute
- any PNG plots created by your diagnostic script contain the provenance information in the 'Image History' attribute

Note that this is done automatically by the ESMValTool Core.

********************************************
How to prepare and run your first diagnostic
********************************************

Instructions for personal diagnostic
====================================

Anyone can run a personal diagnostic, no matter where the location of it;
there is no need to install esmvaltool in developer mode nor is it to
git push or for that matter, do any git operations; the example recipe

.. code-block:: console

  esmvaltool/recipes/recipe_my_personal_diagnostic.yml

shows the use of running a personal diagnostic; the example

.. code-block:: console

  esmvaltool/diag_scripts/examples/my_little_diagnostic.py

and any of its alterations may be used as training wheels for the future ESMValTool
diagnostic developer. The purpose of this example is to familiarize the user with
the framework of ESMValTool without the constraints of installing and running the
tool as developer.

Functionality
=============

`my_little_diagnostic` (or whatever the user will call their diagnostic) makes full use
of ESMValTool's preprocessor output (both phyisical files and run variables); this output
comes in form of a nested dictionary, or config dictionary, see an example below;
it also makes full use of the ability to call any of the preprocessor's functions,
note that relative imports of modules from the esmvaltool package are allowed and
work without altering the $PYTHONPATH.

The user may parse this dictionary so that they execute a number of operations on the
preprocessed data; for example the `my_little_diagnostic.plot_time_series` grabs the
preprocessed data output, computes global area averages for each model, then plots
a time-series for each model. Different manipulation functionalities for grouping,
sorting etc of the data in the config dictionary are available,
please consult ESMValTool User Manual.


Writing a basic recipe
======================
The user will need to write a basic recipe to be able to run their own personal diagnostic.
An example of such a recipe is found in `esmvaltool/recipes/recipe_my_personal_diagnostic.yml`.
For general guidelines with regards to ESMValTool recipes please consult the User Guide;
the specific parameters needed by a recipe that runs a personal diagnostic are:

.. code-block:: yaml

  scripts:
    my_diagnostic:
    script: /path/to/your/my_little_diagnostic.py

i.e. the full path to the personal diagnostic that the user needs to run.

Example of config dictionary
============================
To be added (use python-style code-block).
