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
`ESMValTool tutorial <https://esmvalgroup.github.io/ESMValTool_Tutorial/>`_
on how to run your first recipe. It can be found
`here <https://esmvalgroup.github.io/ESMValTool_Tutorial/04-recipe/index.html>`_.

An
`example recipe <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/examples/recipe_python.yml>`_
is available in the ESMValTool installation folder as
``examples/recipe_python.yml``.

To run this recipe and automatically download the required climate data
from ESGF to the local directory ``~/climate_data``, run

.. code:: bash

	esmvaltool run examples/recipe_python.yml --offline=False

The ``--offline=False`` option tells ESMValTool to search for and download
the necessary datasets.
If you have data available locally, you can run the tool without the
``--offline=False`` argument (the default).
Note that in that case the required data should be located in the directories
specified in your user configuration file.
Recall that the chapter :ref:`Configuring ESMValTool <config-user>`
provides an explanation of how to create your own config-user.yml file.

Available diagnostics and metrics
=================================

If you have the ESMValTool package installed, a number of recipes
are available.
See Section :ref:`recipes` for a description of all
available recipes.

To see a list of installed recipes run

.. code:: bash

	esmvaltool recipes list

To copy an installed recipe to the current working directory, run

.. code:: bash

    esmvaltool recipes get recipe_example.yml

To view an installed recipe on the console, run

.. code:: bash

    esmvaltool recipes show recipe_example.yml

The ``esmvaltool run recipe_example.yml`` command will first look if
``recipe_example.yml`` is the path to an existing file.
If this is the case, it will run that recipe.
If not, it will look if the name matches one of the recipes
in your ESMValTool installation directory, in the subdirectory
`recipes <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/recipes/>`__
and run that.

Running multiple recipes
========================

Have a look at :ref:`running_multiple_recipes` if you are interested in running multiple
recipes in parallel.
