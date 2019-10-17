.. _protips:

************************
Some pro tips on running
************************

Re-running diagnostics
======================

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
===================================

Sometimes it is useful to enter an interactive session to have a look what's going on.
Insert a single line in the code where you want to enter IPython:
``import IPython; IPython.embed()``

This is a useful functionality because it allows the user to `fix` things on-the-fly and after
quitting the Ipython console, code execution continues as per normal.


Use multiple config-user.yml files
==================================

The user selects the configuration yaml file at run time. It's possible to
have several configurations files. For instance, it may be practical to have one
config file for debugging runs and another for production runs.

Create a symbolic link to the latest output directory
=====================================================

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


Running a dry run
=================

You can run in dry-run mode with

.. code:: bash

   esmvaltool -c ~/config-user.yml recipe_xxx.yml --dry-run


This mode activated will run through the data finding and CMOR checks and fixes
and will highlight on screen and in `run/main_log.txt` everytime certain data is
missing or there are issues with the CMOR checks; note that no data is written
to disk and no diagnostics are run; you don't have to modify your recipe in any
way to have this mode run. The information provided will help you obtain any data
that is missing and/or create fixes for the datasets and variables that failed the
CMOR checks and could not be fixed on the fly.
