.. _running:

*******
Running
*******

ESMValTool is mostly used as a command line tool. Whenever your
conda environment for ESMValTool is active, you can just run the command
``esmvaltool <options> <recipe>``. One of the options that must be specified
is the user configuration file, which is specified using the
option ``-c /path/to/config-user.yaml``. An
`example recipe <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_python.yml>`_
is available in the ESMValTool installation folder as
``examples/recipe_python.yml``.

This recipe finds data from CanESM2 and MPI-ESM-LR for 2000 - 2002,
extracts a single level (850 hPa), regrids it to a 1x1 degree mesh and runs
a diagnostic script that creates some plots of Air temperature and
precipitation flux. You can download the recipe from
`github <https://github.com/ESMValGroup/ESMValTool/blob/master/esmvaltool/recipes/examples/recipe_python.yml>`_
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


Running multiple recipes
========================

It is possible to run more tha one recipe in one go: currently this relies on the user
having access to a HPC that has ``rose`` and ``cylc`` installed since the procedure involves
installing and submitting a Rose suite. the utility that allows you to do this is
``esmvaltool/utils/rose-cylc/esmvt_rose_wrapper.py``.

Base suite:
-----------
The base suite to run esmvaltool via rose-cylc is `u-bd684`; you can find
this suite in the Met Office Rose repository at:

https://code.metoffice.gov.uk/svn/roses-u/b/d/6/8/4/trunk/

When ``rose`` will be working with python3.x, this location will become
default and the pipeline will aceess it independently of user, unless, of
course the user will specify ``-s $SUITE_LOCATION``; until then the user needs
to grab a copy of it in ``$HOME`` or specify the default location via ``-s`` option.

Environment:
------------
We will move to a unified and centrally-installed esmvaltool environment;
until then, the user will have to alter the env_setup script:

``u-bd684/app/esmvaltool/env_setup``

with the correct pointers to esmvaltool installation, if desired.

To be able to submit to cylc, you need to have the `/metomi/` suite in path
AND use a `python2.7` environment. Use the Jasmin-example below for guidance.

Jasmin-example:
---------------
This shows how to interact with rose-cylc and run esmvaltool under cylc
using this script:

.. code:: bash

   export PATH=/apps/contrib/metomi/bin:$PATH
   export PATH=/home/users/valeriu/miniconda2/bin:$PATH
   mkdir esmvaltool_rose
   cd esmvaltool_rose
   cp ESMValTool/esmvaltool/utils/rose-cylc/esmvt_rose_wrapper.py .
   svn checkout https://code.metoffice.gov.uk/svn/roses-u/b/d/6/8/4/trunk/ ~/u-bd684
   [enter Met Office password]
   [configure ~/u-bd684/rose_suite.conf]
   [configure ~/u-bd684/app/esmvaltool/env_setup]
   python esmvt_rose_wrapper.py -c config-user.yml \
   -r recipe_autoassess_stratosphere.yml recipe_OceanPhysics.yml \
   -d $HOME/esmvaltool_rose
   rose suite-run u-bd684

Note that you need to pass FULL PATHS to cylc, no `.` or `..` because all
operations are done remotely on different nodes.

A practical actual example of running the tool can be found on JASMIN:
``/home/users/valeriu/esmvaltool_rose``.
There you will find the run shell: ``run_example``, as well as an example
how to set the configuration file. If you don't have Met Office credentials,
a copy of `u-bd684` is always located in ``/home/users/valeriu/roses/u-bd684`` on Jasmin.
