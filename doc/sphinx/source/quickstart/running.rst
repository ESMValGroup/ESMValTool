.. _running:

*******
Running
*******

ESMValTool is mostly used as a command line tool.
Whenever your Conda environment for ESMValTool is active, you can run the
command ``esmvaltool``.
See :ref:`running esmvaltool <esmvalcore:running>` in the ESMValCore
documentation for an introduction to the ``esmvaltool`` command.

Running your first recipe
=========================

There is a step-by-step tutorial available in the
`ESMValTool tutorial <https://tutorial.esmvaltool.org/>`_
on how to run your first recipe. It can be found
`here <https://tutorial.esmvaltool.org/04-recipe/index.html>`_.

An
`example recipe <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_python.yml>`_
is available in the ESMValTool installation folder as
``examples/recipe_python.yml``.

This recipe finds data from BCC-ESM1 and CanESM2 and creates two plot types:

- a global map plot that shows the monthly mean 2m surface air temperature in
  January 2000.
- a time series plot that shows the globally averaged annual mean 2m surface
  air temperature and compares it to the one in Amsterdam.

To run this recipe and automatically download the required climate data
from ESGF to the local directory ``~/climate_data``, run

.. code:: bash

	esmvaltool run examples/recipe_python.yml --search_esgf=when_missing

The ``--search_esgf=when_missing`` option tells ESMValTool to search for and
download the necessary climate data files, if they cannot be found locally.
The data only needs to be downloaded once, every following run will reuse
previously downloaded data.
If you have all required data available locally, you can run the tool with
``--search_esgf=never`` argument (the default).
Note that in that case the required data should be located in the directories
specified in the configuration (see :ref:`esmvalcore:config_option_rootpath`).
A third option ``--search_esgf=always`` is available.
With this option, the tool will first check the ESGF for the needed data,
regardless of any local data availability; if the data found on ESGF is newer
than the local data (if any) or the user specifies a version of the data that
is available only from the ESGF, then that data will be downloaded; otherwise,
local data will be used.
Recall that the chapter on :ref:`configuring ESMValTool <config>`
provides an explanation of how to set up the configuration.

See :ref:`running esmvaltool <esmvalcore:running>` in the ESMValCore
documentation for a more complete introduction to the ``esmvaltool`` command.

.. _recipes_command:

Available diagnostics and metrics
=================================

Although ESMValTool can be used to download data, analyze it using ESMValCore's
preprocessing modules, and the creation of your own analysis code, its main purpose is the
continuously growing set of diagnostics and metrics that it directly provides to
the user. These metrics and diagnostics are provided as a set of preconfigured
recipes that users can run or customize for their own analysis.
The latest list of available recipes can be found :ref:`here <esmvaltool:recipes>`.

In order to make the management of these installed recipes easier, ESMValTool
provides the ``recipes`` command group with utilities that help the users in
discovering and customizing the provided recipes.

The first command in this group allows users to get the complete list of installed
recipes printed to the console:

.. code:: bash

    esmvaltool recipes list

If the user then wants to explore any one of these recipes, they can be printed
using the following command

.. code:: bash

    esmvaltool recipes show recipe_name.yml

Note that there is no ``recipe_name.yml`` shipped with ESMValTool, replace
this with a recipes that is available, for example
`examples/recipe_python.yml <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_python.yml>`_.
Finally, to get a local copy that can then be customized and run, users can
run the following command

.. code:: bash

    esmvaltool recipes get recipe_name.yml

Note that the ``esmvaltool run recipe_name.yml`` command will first look if
``recipe_name.yml`` is the path to an existing file.
If this is the case, it will run that recipe.
If not, it will look if it is a relative path to an existing recipe with respect to the
`recipes <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/>`__
directory in your ESMValTool installation and run that.

Running multiple recipes
========================

Have a look at :ref:`running_multiple_recipes` if you are interested in running multiple
recipes in parallel.
