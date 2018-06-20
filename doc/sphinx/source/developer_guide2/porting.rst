.. _porting:

**************************************************************
Porting namelists (recipes) and diagnostics to ESMValTool v2.0
**************************************************************

This guide summarizes the main steps to be taken in order to port an ESMValTool namelist (now called **recipe**) and the corresponding diagnostic(s) from v1.0 to v2.0, hereafter also referred as the *"old"* and the *"new version"*, respectively. The new ESMValTool version is being developed in the public git branch ``version2_development``. An identical version of this branch is maintained in the private repository as well and kept synchronized on an hourly basis.

In the following, it is assumed that the user has successfully installed ESMValTool v2 and has a rough overview of its structure (see `Technical Overview <http://www.esmvaltool.org/download/Righi_ESMValTool2-TechnicalOverview.pdf>`_).

Create a github issue
=====================

Create an issue in the public repository to keep track of your work and inform other developers. See an example `here <https://github.com/ESMValGroup/ESMValTool/issues/293>`_. Use the following title for the issue: "PORTING <recipe> into v2.0".
Do not forget to assign it to yourself.

Create your own branch
======================

A dedicated branch shall be created from ``version2_development`` for each namelist (recipe) to be ported.

.. code-block:: bash

    git checkout version2_development
    git checkout -b version2_<recipe>

``version2_development`` contains only v2.0 under the ``./esmvaltool/`` directory. 

Convert xml to yml
==================

In ESMValTool v2.0, the namelist (now recipe) is written in yaml format (`Yet Another Markup Language format <http://www.yaml.org/>`_). It may be useful to activate the yaml syntax highlighting for the editor in use. This improves the readability of the recipe file and facilitates the editing, especially concerning the indentations which are essential in this format (like in python). Instructions can be easily found online, for example for `emacs <https://www.emacswiki.org/emacs/YamlMode>`_ and `vim <http://www.vim.org/scripts/script.php?script_id=739>`_.


The yaml recipe can now be edited and tested, starting with a few models and one diagnostics and proceed gradually. The recipe file ``./esmvaltool/recipes/recipe_perfmetrics_CMIP5.yml`` can be used as an example, as it covers most of the common cases.

Create a copy of the diag script in v2.0
========================================

The diagnostic script to be ported shall first be created into ./esmvaltool/diag_script/. It is recommended to get a copy of the very last version of the script to be ported from the development branch (either in the public or in the private repository). Just create a local (offline) copy of this file from the repository and add it to ../esmvaltool/diag_script/ as a new file.
 
Note that (in general) this is not necessary for plot scripts and for the libraries in ``./esmvaltool/diag_script/ncl/lib/``, which have already been ported. Changes may however still be necessary, especially in the plot scripts which have not yet been fully tested with all diagnostics.

Check and apply renamings
=========================

The new ESMValTool version includes a completely revised interface, handling the communication between the python workflow and the (NCL) scripts. This required several variables and functions to be renamed or removed. These chagnes are listed in the following table and shall be applied to the diagnostic code before starting with testing.

.. tabularcolumns:: |p{6cm}|p{6cm}|p{3cm}|

+-------------------------------------------------+-----------------------------------------------------+------------------+
| Name in v1.0                                    | Name in v2.0                                        | Affected code    |
+=================================================+=====================================================+==================+
| ``getenv("ESMValTool_wrk_dir")``                | ``config_user_info@work_dir``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``getenv(ESMValTool_att)``                      | ``diag_script_info@att`` or                         | all .ncl scripts |
|                                                 | ``config_user_info@att``                            |                  |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``xml``                                         | ``yml``                                             | all scripts      |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``var_attr_ref(0)``                             | ``variable_info@reference_model``                   | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``var_attr_ref(1)``                             | ``variable_info@alternative_model``                 | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``models``                                      | ``model_info`` or ``input_file_info``               | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``models@name``                                 | ``model_info@model`` or                             | all .ncl scripts |
|                                                 | ``input_file_info@model``                           |                  |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``verbosity``                                   | ``config_user_info@log_level``                      | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``isfilepresent_esmval``                        | ``fileexists``                                      | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``messaging.ncl``                               | ``logging.ncl``                                     | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``               | ``log_info(arg1)`` if ``arg3=1``                    | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``               | ``log_debug(arg1)`` if ``arg3>1``                   | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``verbosity = config_user_info@verbosity``      | remove this statement                               | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``enter_msg(arg1, arg2, arg3)``                 | ``enter_msg(arg1, arg2)``                           | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``leave_msg(arg1, arg2, arg3)``                 | ``leave_msg(arg1, arg2)``                           | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``noop()``                                      | appropriate ``if-else`` statement                   | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``nooperation()``                               | appropriate ``if-else`` stsatement                  | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``fullpaths``                                   | ``input_file_info@filename``                        | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``get_output_dir(arg1, arg2)``                  | ``config_user_info@plot_dir``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``get_work_dir``                                | ``config_user_info@work_dir``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``inlist(arg1, arg2)``                          | ``any(arg1.eq.arg2)``                               | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load interface_scripts/*.ncl``                | ``load interface_scripts/interface.ncl``            | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``<varname>_info.tmp``                          | ``<varname>_info.ncl`` in ``preproc`` dir           | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``ncl.interface``                               | ``settings.ncl`` in ``run_dir`` and                 | all .ncl scripts |
|                                                 | ``interface_scripts/interface.ncl``                 |                  |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/``                  | ``load diag_scripts/shared/``                       | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load plot_scripts/ncl/``                      | ``load diag_scripts/shared/plot/``                  | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/rgb/``              | ``load diag_scripts/shared/plot/rgb/``              | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/styles/``           | ``load diag_scripts/shared/plot/styles``            | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``load diag_scripts/lib/ncl/misc_function.ncl`` | ``load diag_scripts/shared/plot/misc_function.ncl`` | all .ncl scripts |
+-------------------------------------------------+-----------------------------------------------------+------------------+
| ``LW_CRE``, ``SW_CRE``                          | ``lwcre``, ``swcre``                                | some yml nmls    |
+-------------------------------------------------+-----------------------------------------------------+------------------+

The following changes shall also be considered:

- ``run_dir`` (previous ``interface_data``), ``plot_dir``, ``work_dir`` are now unique to each diagnostic script, so it is no longer necessary to define specific paths in the diagnostic scripts to prevent file collision;
- the interface functions ``interface_get_*`` and ``get_figure_filename`` are no longer available: their functionalities can be easily reproduced using the ``model_info`` and ``input_file_info`` logicals and their attributes;
- there are now only 4 log levels (``debug``, ``info``, ``warning``, and ``error``) instead of (infinite) numerical values in ``verbosity``
- diagnostic scripts are now organized in subdirectories in ``esmvaltool/diag_scripts/``: all scripts belonging to the same diagnostics shall be collected in a single subdirectory (see ``esmvaltool/diag_scripts/perfmetrics/`` for example). This applies also to the ``aux_`` scripts, unless they are shared among multiple diagnostics (in this case they shall go in ``shared/``);
- upper case characters shall be avoided in script names.

As for the recipe, the diagnostic script ``./esmvaltool/diag_scripts/perfmetrics_main.ncl`` can be followed as working example.

Move preprocessing from the diagnostic script to the backend
============================================================

Many operations previously performed by the diagnostic scripts, are now included in the backend, including level extraction, regridding, masking, and multi-model statistics. If the diagnostics to be ported contains code performing any of such operations, the corresponding code has to be removed from the diagnostic script and the respective backend functionality shall be used instead.

The backend operations are fully controlled by the ``preprocessors`` section in the recipe. Here, a number of preprocessor sets can be defined, with different options for each of the operations. The sets defined in this section are applied in the ``diagnostics`` section to preprocess a given variable.

It is recommended to proceed step by step, porting and testing each operation separately before proceeding with the next one. A useful setting in the user configuration file (``config-private.yml``) called ``write_intermediary_cube`` allows writing out the variable field after each preprocessing step, thus facilitating the comparison with the old version (e.g., after CMORization, level selection, after regridding, etc.). The CMORization step of the new backend exactly corresponds to the operation performed by the old backend (and stored in the ``climo`` directory, now called ``preprec``): this shall be the very first step to be checked, by simply comparing the intermediary file produced by the new backend after CMORization with the output of the old backend in the ``climo`` directorsy (see "Testing" below for instructions).

The new backend also performs variable derivation, replacing the ``calculate`` function in the ``variable_defs`` scripts. If the recipe which is being ported makes use of derived variables, the corresponding calculation must be ported from the ``./variable_defs/<variable>.ncl`` file to ``./esmvaltool/preprocessor/_derive.py``.

Note that the Python library ``esmval_lib``, containing the ``ESMValProject`` class is no longer available in version 2. Most functionalities have been moved to the new preprocessor. If you miss a feature, please open an issue on github [https://github.com/ESMValGroup/ESMValTool/issues].

Move diagnostic- and variable-specific settings to the recipe
===============================================================

In the new version, all settings are centralized in the recipe, completely replacing the diagnostic-specific settings in ``./nml/cfg_files/`` (passed as ``diag_script_info`` to the diagnostic scripts) and the variable-specific settings in ``variable_defs/<variable>.ncl`` (passed as ``variable_info``). There is also no distinction anymore between diagnostic- and variable-specific settings: they are collectively defined in the ``scripts`` dictionary of each diagnostic in the recipe and passed as ``diag_script_info`` attributes by the new ESMValTool interface. Note that the ``variable_info`` logical still exists, but it is used to pass variable information as given in the corresponding dictionary of the recipe.

Test the recipe/diagnostic in the new version
===============================================

Once complete, the porting of the diagnostic script can be tested. Most of the diagnostic script allows writing the output in a NetCDF file before calling the plotting routine. This output can be used to check whether the results of v1.0 are correctly reproduced. As a reference for v1.0, it is recommended to use the development branch.

There are two methods for comparing NetCDF files: ``cdo`` and ``ncdiff``. The first method is applied with the command:

.. code-block:: bash

    cdo diffv old_output.nc new_output.nc

which will print a log on the stdout, reporting how many records of the file differ and the absolute/relative differences.

The second method produces a NetCDF file (e.g., ``diff.nc``) with the difference between two given files:

.. code-block:: bash

    ncdiff old_output.nc new_output.nc diff.nc

This file can be opened with ``ncview`` to visually inspect the differences.

In general, binary identical results cannot be expected, due to the use of different languages and algorithms in the two versions, especially for complex operations such as regridding. However, difference within machine precision shall be aimed at. At this stage, it is essential to test all models in the recipe and not just a subset of them.

It is also recommended to compare the graphical output (this may be necessary if the ported diagnostic does not produce a NetCDF output). For this comparison, the PostScript format shall be chosen (it can be set in the user configuration file). Two PostScript files can be compared with standard ``diff`` command in Linux:

.. code-block:: bash

   diff old_graphic.ps new_graphic.ps

but it is very unlikely to produce no differences, therefore visual inspection of the output may also be required.

Clean the code
==============

Before submitting a pull request, the code shall be cleaned to adhere to the coding standard, which are somehow stricter in v2.0. For python code, this check is performed automatically on GitHub (CircleCI and Codacy). For NCL code, this is still done manually and considers the following guidelines:

- code syntax shall be checked using ``/util/ncl-checker/pep8.py <diag>.ncl`` and all reported warnings shall be fixed;
- two-space instead of four-space indentation is now adopted for NCL as per NCL standard;
- ``load`` statements for NCL standard libraries shall be removed: these are automatically loaded since NCL v6.4.0 (see `NCL documentation <http://www.ncl.ucar.edu/current_release.shtml#PreloadedScripts6.4.0>`_);
- the description of diagnostic- and variable-specific settings shall be moved from the header of the diagnostic script to the recipe, since the settings are now defined there (see above);
- NCL ``print`` and ``printVarSummary`` statements shall be avoided and replaced by the ``log_info`` and ``log_debug`` functions;
- for error and warning statments, the ``error_msg`` function shall be used, which automatically include an exit statement.

Open a pull request
===================

Create a pull request on github to merge your branch back to ``version2_development``, provide a short description of what has been done and nominate one or more reviewers.
