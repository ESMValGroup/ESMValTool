:: _running:

******************
Running ESMValTool
******************

To run ESMValTool, use the command

.. code:: bash

	esmvaltool -c /path/to/config-user.yml examples/recipe_python.yml

This will run the example recipe_python.yml. The path to the recipe
can either be the path to a recipe file, or a path relative to the
esmvaltool/recipes directory of your installed ESMValTool. See the chapter
`User configuration file`_ for an explanation of how to create your own
config-user.yml file.

To get help on additional commands, please use

.. code:: bash

	esmvaltool --help



Available diagnostics and metrics
=================================

See :ref:`recipes` for a description of all available recipes.
