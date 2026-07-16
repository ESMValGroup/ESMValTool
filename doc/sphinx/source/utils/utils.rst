.. _utils:

Utilities
*********

This section provides information on tools that are useful when developing
ESMValTool.
Tools that are specific to ESMValTool live in the
`esmvaltool/utils <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/utils>`_
directory, while others can be installed using the usual package managers.

.. _pre-commit:

Pre-commit
==========

`pre-commit <https://pre-commit.com/>`__ is a handy tool that can run many
tools for checking code quality with a single command.
Usually it is used just before committing, to avoid accidentally committing
mistakes.
It knows knows which tool to run for each filetype, and therefore provides
a convenient way to check your code!


To run ``pre-commit`` on your code, go to the ESMValTool directory
(``cd ESMValTool``) and run

::

   pre-commit run

By default, pre-commit will only run on the files that have been changed,
meaning those that have been staged in git (i.e. after
``git add your_script.py``).

To make it only check some specific files, use

::

   pre-commit run --files your_script.py

or

::

   pre-commit run --files your_script.R

Alternatively, you can configure ``pre-commit`` to run on the staged files before
every commit (i.e. ``git commit``), by installing it as a `git hook <https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks>`__ using

::

   pre-commit install

Pre-commit hooks are used to inspect the code that is about to be committed. The
commit will be aborted if files are changed or if any issues are found that
cannot be fixed automatically. Some issues cannot be fixed (easily), so to
bypass the check, run

::

   git commit --no-verify

or

::

   git commit -n

or uninstall the pre-commit hook

::

   pre-commit uninstall


Note that the configuration of pre-commit lives in
`.pre-commit-config.yaml <https://github.com/ESMValGroup/ESMValTool/blob/main/.pre-commit-config.yaml>`_.

.. _nclcodestyle:

nclcodestyle
============

A tool for checking the style of NCL code, based on pycodestyle.
Install ESMValTool in development mode (``pip install -e '.[develop]'``) to make it available.
To use it, run

.. code-block:: bash

    nclcodestyle /path/to/file.ncl

.. _recipe_test_tool:

Colormap samples
================
Tool to generate colormap samples for ESMValTool's default Python and NCL colormaps.

Run

.. code-block:: bash

    esmvaltool colortables python

or

.. code-block:: bash

    esmvaltool colortables ncl

to generate the samples.

.. _running_multiple_recipes:

Running multiple recipes
========================

It is possible to run more than one recipe in one go.

This can for example be achieved by using ``rose`` and/or ``cylc``, tools
that may be available at your local HPC cluster.

In the case in which neither ``rose`` nor ``cylc`` are available at your HPC cluster,
it is possible to automatically generate job submission scripts, as well as a summary of the
job outputs using the scripts available in
`esmvaltool/utils/batch-jobs <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/batch-jobs>`__.

Using cylc
----------

A cylc suite for running all recipes is available in
`esmvaltool/utils/testing/regression <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/testing/regression>`__.
This suite is configured to work with versions of cylc older than 8.0.0 .

To prepare for using this tool:

#. Log in to a system that uses `slurm <https://slurm.schedmd.com/quickstart.html>`_
#. Make sure the required CMIP and observational datasets are available and
   their ``rootpath`` and ``drs`` is properly set up in the :ref:`configuration
   <esmvalcore:config_options>`
#. Make sure the required auxiliary data is available (see :ref:`recipe documentation <recipes>`)
#. Install ESMValTool

Next, get started with `cylc <https://cylc.github.io/cylc-doc/7.9.3/html/index.html>`_:

#. Run ``module load cylc``
#. Register the suite with cylc ``cylc register run-esmvaltool-recipes ~/ESMValTool/esmvaltool/utils/testing/regression``
#. Edit the suite if needed, this allows e.g. choosing which recipes will be run
#. Validate the suite ``cylc validate run-esmvaltool-recipes --verbose``, this will e.g. list the recipes in the suite
#. Run all recipes ``cylc run run-esmvaltool-recipes``
#. View progress ``cylc log run-esmvaltool-recipes``, use e.g. ``cylc log run-all-esmvaltool-recipes examples-recipe_python_yml.1 --stdout`` to see the log of an individual esmvaltool run. Once the suite has finished running, you will see the message "WARNING - suite stalled" in the log.
#. Stop the cylc run once everything is done ``cylc stop run-esmvaltool-recipes``.

To generate an overview page of the recipe runs, use the ``summarize.py`` :ref:`utility script <overview_page>`.

.. _utils_batch_jobs:

Using the scripts in `utils/batch-jobs`
---------------------------------------

In `utils/batch-jobs <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/batch-jobs>`_,
you can find a script to generate slurm submission scripts for all available recipes in ESMValTool,
as well as a script to parse the job outputs.

.. _utils_generate:

Using `generate.py`
...................

The script `generate.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/batch-jobs/generate.py>`_,
is a simple python script that creates slurm submission scripts, and
if configured, submits them to the HPC cluster. It has been tested in `DKRZ's Levante cluster <https://docs.dkrz.de/doc/levante/index.html>`_.

The following parameters have to be set in the script in order to make it run:

* ``env``, *str*: Name of the conda environment in which `esmvaltool` is installed.
* ``mail``, *bool*: Whether or not to receive mail notifications when a submitted job fails or finishes successfully. Default is ``False``.
* ``submit``, *bool*: Whether or not to automatically submit the job after creating the launch script. Default value is ``False``.
* ``account``, *str*: Name of the DKRZ account in which the job will be billed.
* ``outputs``, *str*: Name of the directory in which the job outputs (.out and .err files) are going to be saved. The outputs will be saved in `/home/user/<outputs>`.
* ``conda_path``, *str*: Full path to the `miniforge3/etc/profile.d/conda.sh` executable.

Optionally, the following parameters can be edited:

* ``config_dir``, *str*: Path to :ref:`configuration directory <esmvalcore:config_yaml_files>`, by default ``~/.config/esmvaltool/``.
* ``partition``, *str*: Name of the DKRZ partition used to run jobs. Default is ``interactive`` to minimize computing cost compared to ``compute`` for which nodes cannot be shared.
* ``memory``, *str*: Amount of memory requested for each run. Default is ``64G`` to allow to run 4 recipes on the same node in parallel.
* ``time``, *str*: Time limit. Default is ``04:00:00`` to increase the job priority. Jobs can run for up to 8 hours and 12 hours on the compute and interactive partitions, respectively.
* ``default_max_parallel_tasks``, *int*: Default is ``8`` which works for most recipes. For other cases, an entry needs to be made to the ``MAX_PARALLEL_TASKS`` dictionary (see below).

The script will generate a submission script for all recipes using by default the ``interactive`` queue and with a time limit of 4h. In case a recipe
may require of additional resources, they can be defined in the ``SPECIAL_RECIPES`` dictionary. The recipe name has to be given as a ``key`` in which the
values are another dictionary.
The latter are used to specify the ``partition`` in which to submit the recipe, the new ``time`` limit and other ``memory`` requirements
given by the slurm flags ``--mem``, ``--constraint`` or ``--ntasks``. In general, an entry in ``SPECIAL_RECIPES`` should be set as:

.. code-block:: python

   SPECIAL_RECIPES = {
    'recipe_name': {
        'partition': '#SBATCH --partition=<name_of_the_partition>',
        'time': '#SBATCH --time=<custom_time_limit>',
        'memory': '#SBATCH --mem=<custom_memory_requirement>' # --constraint or --nstasks can be used instead.
        },
   }

Some recipes can only be run with a number of tasks less than ``default_max_parallel_tasks`` for various reasons (memory issues, diagnostic issues, CMIP3 data used).
These recipes need to be added to the ``MAX_PARALLEL_TASKS`` dictionary with a specific ``max_parallel_tasks`` value.

Note that the script has been optimized to use standard SLURM settings to run most recipes while minimizing the computational cost of the jobs and tailored runtime settings for resource-intensive recipes.
It is only necessary to edit this script for recipes that have been added since the last release and cannot be run with the default settings.

In the case in which ``submit`` is set to ``True``, but you want to exclude certain recipes from being submitted, their name can be added in the ``exclude`` list:

.. code-block:: python

   exclude = ['recipe_to_be_excluded_1', 'recipe_to_be_excluded_2']

.. _utils_parse:

Using `parse_recipes_outputs`
.............................

You can run this script (simply as a standalone Python script) after all recipes have been run, to gather a bird's eye view
of the run status for each recipe; running the script provides you with a Markdown-formatted list of recipes that succeeded,
recipes that failed due to a diagnostic error, and recipes that failed due to missing data (the two most common causes for
recipe run failure). You should provide the location of the output log files from SLURM (``*.out`` and ``*.err``) to the
script as well as a list of all available recipes. To generate the list, run the command:

.. code-block:: bash

   for recipe in $(esmvaltool recipes list | grep '\.yml$'); do echo $(basename "$recipe"); done > all_recipes.txt

To keep the script execution fast, it is recommended to use ``log_level: info`` in the configuration so that SLURM
output files are rather small.

.. _overview_page:

Overview of recipe runs
=======================

To create overview webpages of a set of recipe runs, run:

.. code-block:: python

   python esmvaltool/utils/testing/regression/summarize.py ~/esmvaltool_output/

This will generate 2 html files:

-  ``index.html`` that displays a summary of each recipe run, with a title and
   a representative plot, a short description of the aim of the recipe, and
   links to each individual run.
-  ``debug.html`` that provides an overview table of successful and failed runs
   with links to each individual run, and computing resources used for each run.

.. _compare_recipe_runs:

Comparing recipe runs
=====================

A command-line tool is available for comparing one or more recipe runs to
known good previous run(s).
This tool uses `xarray <https://docs.xarray.dev/en/stable/>`_ to compare NetCDF
files and difference hashing provided by
`imagehash <https://pypi.org/project/ImageHash/>`_ to compare PNG images.
All other file types are compared byte for byte.

The package imagehash_ is included in the standard ESMValTool environment.
However, if you are using an environment that does not include imagehash,
install it as follows:

.. code-block:: bash

   pip install imagehash

Next, run the ``esmvaltool develop compare`` command:

.. code-block:: bash

    esmvaltool develop compare ~/reference_output/recipe_python_20220109_110112/ ~/output/recipe_python_20220310_180417/

where the first argument is a directory containing the output
from a reference run and the second argument is a directory
containing output from a test run to compare against the reference.
This command will exit with a return code of 0 if the run outputs :w
are the same,
or 1 if they are different.

To compare results from a set of multiple runs to a reference set,
e.g. between the current version of the code to the previous version,
there is a wrapper script that loops over multiple comparisons.
To use this go to the location where ESMValTool is installed and run e.g.:

.. code-block:: bash

    python esmvaltool/utils/testing/regression/compare.py /shared/esmvaltool/v2.4.0 /shared/esmvaltool/v2.5.0

To get more information on how a result is different, use the
``--verbose`` flag with either of the above commands.

.. _draft_release_notes.py:

draft_release_notes.py
======================

`draft_release_notes.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/utils/draft_release_notes.py>`__
is a script for drafting release notes based on the titles and labels of
the GitHub pull requests that have been merged since the previous release.

To use it, install the package pygithub_:

.. code-block:: bash

   pip install pygithub

Create a `GitHub access token`_ (leave all boxes for additional
permissions unchecked) and store it in the file ``~/.github_api_key``.

Edit the script and update the date and time of the previous release and run
the script:

.. code-block:: bash

   python esmvaltool/utils/draft_release_notes.py ${REPOSITORY}

``REPOSITORY`` can be either ``esmvalcore`` or ``esmvaltool`` depending on the
release notes you want to create.

Review the resulting output (in ``.rst`` format) and if anything needs changing,
change it on GitHub and re-run the script until the changelog looks acceptable.
In particular, make sure that pull requests have the correct label, so they are
listed in the correct category.
Finally, copy and paste the generated content at the top of the changelog.

Converting Version 1 Namelists to Version 2 Recipes
===================================================

The
`xml2yml <https://github.com/ESMValGroup/ESMValTool/tree/main/esmvaltool/utils/xml2yml>`_
converter can turn the old xml namelists into new-style yml
recipes. It is implemented as a xslt stylesheet that needs a processor
that is xslt 2.0 capable. With this, you simply process your old
namelist with the stylesheet xml2yml.xsl to produce a new yml recipe.

After the conversion you need to manually check the mip information in
the variables! Also, check the caveats below!

Howto
-----

One freely available processor is the Java based
`saxon <http://saxon.sourceforge.net/>`__. You can download the free he
edition
`here <https://sourceforge.net/projects/saxon/files/latest/download?source=files>`__.
Unpack the zip file into a new directory. Then, provided you have Java
installed, you can convert your namelist simply with:

::

   java -jar $SAXONDIR/saxon9he.jar -xsl:xml2yml.xsl -s:namelist.xml -o:recipe.yml

Caveats/Known Limitations
-------------------------

-  At the moment, not all model schemes (OBS, CMIP5, CMIP5_ETHZ…) are
   supported. They are, however, relatively easy to add, so if you need
   help adding a new one, please let me know!
-  The documentation section (namelist_summary in the old file) is not
   automatically converted.
-  In version 1, one could name an exclude, similar to the reference
   model. This is no longer possible and the way to do it is to include
   the models with another ``additional_models`` tag in the variable
   section. That conversion is not performed by this tool.

Authored by **Klaus Zimmermann**, direct questions and comments to
klaus.zimmermann@smhi.se

.. _GitHub access token: https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
.. _pygithub: https://pygithub.readthedocs.io/en/latest/introduction.html


Extracting a list of input files from the provenance
====================================================

There is a small tool available to extract just the list of input files used to generate
a figure from the ``*_provenance.xml`` files (see :ref:`recording-provenance` for more
information).

To use it, install ESMValTool from source and run

.. code-block:: bash

    python esmvaltool/utils/prov2files.py /path/to/result_provenance.xml

The tool is based on the `prov <https://prov.readthedocs.io/en/latest/readme.html>`_
library, a useful library for working with provenance files.
With minor adaptations, this script could also print out global attributes
of the input NetCDF files, e.g. the tracking_id.

Recipe Test Workflow (RTW)
==========================

.. include:: RTW/common.txt

The Recipe Test Workflow (|RTW|) is a workflow that is used to regularly run
recipes so issues can be discovered during the development process sooner
rather than later.

.. toctree::
   :maxdepth: 1

   RTW/index
