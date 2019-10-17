.. _running:

******************
Running ESMValTool
******************

ESMValTool is mostly used as a command line tool. Whenever your
conda environment for ESMValTool is active, you can just run the command
``esmvaltool <options> <recipe>``. One of the options that must be specified
is the user configuration file, which is specified using the
option ``-c /path/to/config-user.yaml``. An
`example recipe <https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/recipes/examples/recipe_python.yml>`_
is available in the ESMValTool installation folder as
``examples/recipe_python.yml``.

This recipe finds data from CanESM2 and MPI-ESM-LR for 2000 - 2002,
extracts a single level (850 hPa), regrids it to a 1x1 degree mesh and runs
a diagnostic script that creates some plots of Air temperature and
precipitation flux. You can download the recipe from
`github <https://github.com/ESMValGroup/ESMValTool/blob/version2_development/esmvaltool/recipes/examples/recipe_python.yml>`_
and save it in your project directory as (e.g.) ``recipe_python.yml``
and then run ESMValTool with

.. code:: bash

	esmvaltool -c /path/to/config-user.yml recipe_python.yml --synda

The ``--synda`` option tells ESMValTool to use Synda to search for and download
the necessary datasets.

ESMValTool will also find recipes that are stored in its installation directory.
A copy of the example recipe is shipped with ESMValTool as:
``/path/to/installation/esmvaltool/recipes/examples/recipe_python.yml``.
Thus, the following also works:

.. code:: bash

	esmvaltool -c /path/to/config-user.yml examples/recipe_python.yml

Note that this command does not call Synda. The required data should thus be
located in the directories specified in your user configuration file.
Recall that the chapter :ref:`Configuring ESMValTool <config-user>`
provides an explanation of how to create your own config-user.yml file.

To get help on additional commands, please use

.. code:: bash

	esmvaltool --help


Available diagnostics and metrics
=================================

See Section :doc:`Recipes <../recipes/index>` for a description of all
available recipes.
