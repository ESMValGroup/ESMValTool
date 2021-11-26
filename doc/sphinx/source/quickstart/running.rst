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

.. _recipes_command:

Available diagnostics and metrics
=================================

Although ESMValTool can be used just to simplify the management of data
and the creation of your own analysis code, one of its main strengths is the
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

If the user then wants to explore any one of this recipes, they can be printed
using the following command

.. code:: bash

    esmvaltool recipes show recipe_name.yml

And finally, to get a local copy that can then be customized and run, users can
use the following command

.. code:: bash

    esmvaltool recipes get recipe_name.yml

Note that the ``esmvaltool run recipe_example.yml`` command will first look if
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
