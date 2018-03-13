.. _porting:

Porting namelists and diagnostics to ESMValTool v2.0
****************************************************

This guide summarizes the main steps to be taken in order to port an ESMValTool namelist and the corresponding diagnostics from v1.0 to v2.0. ESMValTool v2.0 is being developed in the public git branch REFACTORING_backend. It is strongly recommended to create a branch from REFACTORING_backend for each of the namelist to be ported and to name it REFACTORING_<namelist>. 

REFACTORING_backend contains both v1.0 and v2.0, the latter within the esmvaltool/ directory. It is therefore possible, and recommended, to run both versions of the ESMValTool within the same branch: this will facilitate testing and comparison of the two version as long as the porting process proceeds.


Convert xml to yml
==================

In ESMValTool v2.0, the main namelist (hereafter *"the namelist"*) is written in yaml format (`Yet Another Markup Language format <http://www.yaml.org/>`_). It may be useful to activate syntax highlighting for the editor in use. This improves the readability of the namelist file and facilitates the editing, especially concerning the indentations which are essential in the yml format (like in python). Instructions can be easily found online, for example for `emacs <https://www.emacswiki.org/emacs/YamlMode>`_ and `vim <http://www.vim.org/scripts/script.php?script_id=739>`_.

For each of the ESMValTool v1.0 namelists, a very first draft in yml format has already been created and is available in ./nml/. This can be used as a starting point, but keeping in mind that it has been created at a very early stage of v2.0 development and it will certainly need further changes.

The namelist file in yml shall first be moved to the esmvaltool/ directory where the v2.0 is being developed::

        mv ./nml/<namelist>.yml ./esmvaltool/namelists/<namelist>.yml


This will help to keep track of which namelists have been already ported to the new version.

The yaml namelist can now be edited and tested, starting with a few models and one diagnostics and proceed gradually. The namelist file ``namelists/namelist_perfmetrics_CMIP5.yml`` can be used as an example, as it covers most of the common cases. A corresponding version of the v1.0 namelist (xml) shall also be developed in parallel in order to compare the output of the two versions and to make sure that the porting has been successful (see "Testing" below).


Create a copy of the diag script in v2.0
========================================

As for the namelist, a copy of the diagnostic script(s) to be ported shall be created in the esmvaltool/ directory, again to allow for direct comparison of the two versions within the same branch::

    cp -i diag_scripts/<diag>.ncl esmvaltool/diag_scripts/


Note that this is not necessary for plot scripts and for the libraries in ncl/lib/, which have already been copied. Changes may however still be necessary, especially in the plot scripts which have not yet been tested with all diagnostics.


Check and apply renamings
=========================

ESMValTool v2.0 includes a completely revised interface handling the communication between the python workflow and the NCL scripts. This required renaming several variable and function names in most of the NCL scripts, which are listed in the following table and shall be applied to the diagnostic code before starting with testing.

+--------------------------------------------+-------------------------------------------+------------------+
| Name in v1.0                               | Name in v2.0                              | Affected code    |
+============================================+===========================================+==================+
| ``getenv("ESMValTool_wrk_dir")``           | ``config_user_info@work_dir``             | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``getenv(ESMValTool_att)``                 | ``diag_script_info@att`` or               | all .ncl scripts |
|                                            | ``config_user_info@att``                  |                  |
+--------------------------------------------+-------------------------------------------+------------------+
| ``xml``                                    | ``yml``                                   | all scripts      |
+--------------------------------------------+-------------------------------------------+------------------+
| ``var_attr_ref(0)``                        | ``variable_info@reference_model``         | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``var_attr_ref(1)``                        | ``variable_info@alternative_model``       | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``models``                                 | ``model_info`` or ``input_file_info``     | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``models@name``                            | ``model_info@model`` or                   | all .ncl scripts |
|                                            | ``input_file_info@model``                 |                  |
+--------------------------------------------+-------------------------------------------+------------------+
| ``verbosity``                              | ``config_user_info@log_level``            | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``isfilepresent_esmval``                   | ``fileexists``                            | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``messaging.ncl``                          | ``logging.ncl``                           | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``          | ``log_info(arg1)`` if ``arg3=1``          | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``info_output(arg1, arg2, arg3)``          | ``debug_info(arg1)`` if ``arg3>1``        | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``verbosity = config_user_info@verbosity`` | remove this statement                     | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``enter_msg(arg1, arg2, arg3)``            | ``enter_msg(arg1, arg2)``                 | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``leave_msg(arg1, arg2, arg3)``            | ``leave_msg(arg1, arg2)``                 | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``noop()``                                 | appropriate ``if-else`` statement         | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``nooperation()``                          | appropriate ``if-else`` stsatement        | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``fullpaths``                              | ``input_file_info@filename``              | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``get_output_dir(arg1, arg2)``             | ``config_user_info@plot_dir``             | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``get_work_dir``                           | ``config_user_info@work_dir``             | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``inlist(arg1, arg2)``                     | ``any(arg1.eq.arg2)``                     | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``load interface_scripts/*.ncl``           | ``load interface_scripts/interface.ncl``  | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``<varname>_info.tmp``                     | ``<varname>_info.ncl`` in ``preproc`` dir | all .ncl scripts |
+--------------------------------------------+-------------------------------------------+------------------+
| ``ncl.interface``                          | ``settings.ncl`` in ``run_dir`` and       | all .ncl scripts |
|                                            | ``interface_scripts/interface.ncl``       |                  |
+--------------------------------------------+-------------------------------------------+------------------+ 

The following changes shall also be considered: 

- ``run_dir`` (previous ``interface_data``), ``plot_dir``, ``work_dir`` are now unique to each diagnostic script, so there is no need to define specific paths within the diagnostic script to prevent file collision anymore;
- the interface functions ``interface_get_*`` and ``get_figure_filename`` are no longer available: their functionality can be easily reproduced using the ``model_info`` and ``input_file_info`` logicals and their attributes;
- there are now only 4 log levels (``debu``,``info``, ``warning``, and ``error``) instead of (infinite) numerical values.

As for the namelist, the diagnostic script ``diag_scripts/perfmetrics_main.ncl`` can be followed as working example.


Move preprocessing from the diagnostic script to the backend
============================================================

Many operations were previously performed by the diagnostic scripts, are now included in the backend, including level extraction, regridding, masking, and multi-model statistics. If the diagnostics to be ported contains code performing any of these operations, this code has to be removed from the diagnostic script and the corresponding backend functionality shall be used instead.

The backend operations are fully controlled by the ``preprocessors`` section in the namelist. Here a number of preprocessor sets can be defined, with different options for each of the operations. The sets defined in this section are referred to in the ``diagnostics`` section to process a given variable.

It is recommended to proceed step by step, porting and testing each operation separately before proceeding with the next one. A useful setting in the configuration file ``config-private.yml`` called ``write_intermediary_cube`` allows writing out the variable field after each preprocessing step, thus facilitating the comparison with the old version (e.g., after CMORization, level selection, after regridding, etc.). The CMORization step of the new backend correspond exactly to the operation performed by the old backend (and stored in the ``climo``-files): this shall be the very first step to be checked, by simply comparing the intermediary file produced by the new backend after CMORization with the output of the old backend in the ``climo`` directorsy.

The new backend also performs variable derivation, replacing the ``calculate`` function in the ``variable_defs`` scripts. If the namelist being portedmakes use of derived variables, the corresponding calculation must be ported from the ``variable_defs`` file to ``esmvaltool/preprocessor/_derive.py``.
