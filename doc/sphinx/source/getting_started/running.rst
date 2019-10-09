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


Some pro tips on running
========================

Re-running diagnostics
----------------------
If a diagnostic fails, you will get the message

.. code:: bash

   INFO    To re-run this diagnostic script, run:

If you run the command in the stdout you will be able to re-run the
diagnostic without having to re-run the whole preprocessor. If you add the ``-f``
argument that will force an overwrite, and it will deletes not just the failed diagnostic,
but the contents of the entire ``plots`` directory - this is risky but useful when needing to
redo the whole work. Adding ``-i`` or ``--ignore-existing`` will not delete any existing files,
and it will re-do the work related to the failed diagnostic only.


Enter interactive mode with iPython
-----------------------------------
Sometimes it is useful to enter an interactive session to have a look what's going on.
Insert a single line in the code where you want to enter IPython:
``import IPython; IPython.embed()``

This is a useful functionality because it allows the user to `fix` things on-the-fly and after
quitting the Ipython console, code execution continues as per normal.


Use multiple config-user.yml files
----------------------------------

The user profiles the configuration yaml file at run time. It's possible to 
have several configurations files. For instance, it may be practical to have one 
config file for debugging runs and another for production runs.

Create a symbolic link to the latest output directory
-----------------------------------------------------
When running multiple times the same recipe, the tool creates separate output directories
sorted by the time tag that they were created at; sometimes, when running quite a few times,
it is not straightforward to detect which one is the `latest` output directory, so a symbolic
link attached to it would make things more clear e.g.:

.. code:: bash

   recipe_example_20190905_163431
   recipe_example_20190905_163519
   recipe_example_LATEST -> recipe_example_20190905_163519


You can do that by running the tool with a recursive sift attached to it:

.. code:: bash

   esmvaltool -c ~/config-user.yml recipe_xxx.yml && \
   ln -sfT $(ls -1d ~/esmvaltool_output/recipe_xxx_* | tail -1) latest
