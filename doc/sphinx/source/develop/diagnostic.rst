Diagnostic
**********

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


Example of config dictionary
============================
To be added (use python-style code-block).
