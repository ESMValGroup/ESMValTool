.. _running:

******************
Running ESMValTool
******************

To run ESMValTool, use the command

.. code:: bash

	esmvaltool -c /path/to/config-user.yml examples/recipe_python.yml

This will run the example recipe_python.yml. The path to the recipe can either
be the path to a recipe file, or a path relative to the esmvaltool/recipes
directory of your installed ESMValTool. See the chapter :ref:`User
configuration file <config-user>` for an explanation of how
to create your own config-user.yml file.

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
