.. _detailed_release_procedure:

Release: recipes runs and comparison
====================================

The release procedure for ESMValTool is a fairly involved process (at the moment), so it
is important to be very well organized and to have documented each procedural steps, so that
the next release manager can follow said steps, and finalize the release without any delays.

The workflow below assumes an ESMValCore release candidate, or a completed stable release, have been released
and deployed on conda-forge and PyPI; it also assumes the release manager has access to accounts on `DKRZ/Levante
<https://docs.dkrz.de/>`_.

Below is a list of steps that the release manager, together with the previous release manager, should go through before the actual release;
these include testing the new code by running all available recipes in the ``main`` branch, and comparing the output against
the previous release.

Open an issue on GitHub
-----------------------

First, open an issue on GitHub where the testing workflow before the release is documented (see example https://github.com/ESMValGroup/ESMValTool/issues/2881).
Name it something relevant like "Recipe testing and comparison for release 2.x.x", and populate the issue description with information
about where the testing is taking place, what tools are used, and what versions, here are some suggestions:


- path to the output directories on DKRZ/Levante

We should document various utilities' versions so that the work can be reproduced in case there
is an issue, or release work needs to be picked up mid-release by another release manager:

- documenting `conda`/`mamba` versions:

.. code-block:: bash

  mamba --version

- documenting `git` branch and its state:

.. code-block:: bash

  git status

Furthermore, the runtime environment needs to be documented: make a copy of the environment file,
and attach it in the release testing issue; to record the environment in a yaml file use e.g.

.. code-block:: bash

  conda env export > ToolEnv_v2xx_Test.txt

Modifications to configuration files need to be documented as well.
To test recipes, it is recommended to only use the default options and DKRZ data directories, simply by uncommenting
the DKRZ-Levante block of a :ref:`newly generated configuration file <esmvalcore:config_yaml_files>`.

Submit run scripts - test recipe runs
-------------------------------------

We are now ready to start running all the available recipes, to compare output against previous release. Running is currently done
via batch scripts submitted to a schedulers (SLURM). Generate the submission scripts using the ``generate.py`` :ref:`utility Python script <utils_generate>`.

You will have to set the name of your environment, your email address (if you want to get email notifications for successful/failed jobs) and the name of the directory you want to store the log files of the jobs. A compute project from which resources are billed needs to be set, and the default partition is set to `interactive`.
More information on running jobs with SLURM on DKRZ/Levante can be found in the DKRZ `documentation
<https://docs.dkrz.de/doc/levante/running-jobs/index.html>`_.

You can also specify the path to your configuration directory where ``max_parallel_tasks`` can be set in a YAML file. The script was found to work well with ``max_parallel_tasks=8``. Some recipes need to be run with ``max_parallel_tasks=1`` (large memory requirements, CMIP3 data, diagnostic issues, ...). These recipes are listed in `ONE_TASK_RECIPES`.

Some recipes need other job requirements, you can add their headers in the `SPECIAL_RECIPES` dictionary. Otherwise the header will be written following the template that is written in the lines below. If you want to exclude recipes, you can do so by uncommenting the `exclude` lines.

Before submitting all jobs, it is recommended to try the batch script generation with ``submit = False`` and check the generated files. If recipes with special runtime requirements have been added to ESMValTool since the previous release, these may need to be added to `SPECIAL_RECIPES` and/or to `ONE_TASK_RECIPES`.
Other recipes should run successfully with the default SLURM settings set in this script.

The launch scripts will be saved in the same directory you execute the script from. These are named like ``launch_recipe_<name>.sh``.
To submit these scripts to the SLURM scheduler, use the ``sbatch launch_recipe_<name>.sh`` command. You can check the status of your BATCH queue by invoking:

.. code-block:: bash

  squeue -u $USER

Also, for computationally-heavy recipes, you can require more memory and/or time, see e.g. edited batch header below
(note the `compute` partition which is used for such heavy runs):

.. code-block:: bash

  #SBATCH --partition=compute
  #SBATCH --time=08:00:00
  #SBATCH --mem=0
  #SBATCH --constraint=512G

.. note::

  On DKRZ/Levante, a user can't have more than 20 SLURM jobs running at a time.
  As soon as a job is finished, the next one should start. More information on the job handling at DKRZ `here
  <https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html#levante-partitions-and-limits>`_.
  Also note that the ``--mem=0`` argument needs be specified if any of the ``--constraint`` arguments are
  used for memory requests, so that the node's full memory is allocated.

Analyse the results
-------------------

Once all jobs are completed, assemble some statistics so that issues with certain recipes
can be followed-up, and document this information in the release issue, such as:

- number of successfully run recipes
- number of failed recipes with preprocessor errors (can they be fixed? Can the fixes be included in the release?)
- number of failed recipes with diagnostic errors (can they be fixed? Can the fixes be included in the release?)
- number of recipes that are missing data
- number of recipes that have various other issues (and document them)

To parse the output of all these runs, use the ``parse_recipes_output.py`` :ref:`utility Python script <utils_parse>`.
It is recommended to run the recipes with `log_level: info` in your config file to enable the parsing script to run fast.

Share the results with the community
------------------------------------

Create the debug.html and index.html overview webpages by running the :ref:`utility script <overview_page>`
in the directory containing the recipe runs.
These two files, together with the recipe output, need to be copied to the disk of a virtual machine (VM)
used to display recipe output in `webpages
<https://esmvaltool.dkrz.de/shared/esmvaltool/>`_.
Do not store final release results on the VM including `/preproc/` dirs, the total
size for all the recipes output, including `/preproc/` dirs is in the 4.5TB ballpark,
much too high for the VM storage capacity! Therefore, we would recommend using the option
to remove preprocessing directories upon recipe running successfully ``--remove-preproc-dir=True``
at runtime, or set ``remove_preproc_dir: true`` in the configuration file.

Login and access to the DKRZ esmvaltool VM - results from recipe runs
are stored on the VM; log in to the Levante head node and then continue to the VM with:

.. code-block:: bash

  ssh user@esmvaltool.dkrz.de

where `user` is your DKRZ/Levante user name.
Then create a new subdirectory in ``/shared/esmvaltool/`` that will contain recipe output.
This should be named like the ESMValCore version used for the testing, e.g. ``v2.x.xrcx`` or  ``v2.x.x``.
Recipe output can be copied by doing from the VM:

.. code-block:: bash

  nohup rsync --exclude preproc/ -rlt /path_to_testing/esmvaltool_output/* /shared/esmvaltool/v2.x.x/

By copying the debug.html and index.html files into /shared/esmvaltool/v2.x.x/, the output
becomes available online, see for `example
<https://esmvaltool.dkrz.de/shared/esmvaltool/v2.7.0>`_.
Before copying the recipe output to the VM, you may want to clean up your directory containing
the results by removing any large ``preproc`` directories of failed runs and only keeping the last run for each recipe.
This will help generating a clearer overview webpage.
Note that the ``summarize.py`` script needs to be rerun if recipe runs were added or deleted
from your testing directory.

Link the overview webpage to the release issue.
This makes it much easier to ask for feedback from recipe developers and analyse failures.

Results produced with the final ESMValCore release candidate should be put in a VM directory
named after the version number, e.g. ``v2.x.x``.
Once the release process is over, test results produced with previous release candidates can be deleted to save space on the VM.

.. note::

  If you wrote recipe runs output to Levante's `/scratch` partition, be aware that
  the data will be removed after two weeks, so you will have to quickly move the
  output data to the VM, using the ``nohup`` command above.

Running the comparison
----------------------

To compare the newly produced output from running all recipes, follow these steps below.

Access to the DKRZ esmvaltool VM, then install miniconda on the VM, and
if you have a Miniconda installer already downloaded in your Levante $HOME

.. code-block:: bash

  scp Miniconda3-py39_4.12.0-Linux-x86_64.sh user@esmvaltool.dkrz.de:/mnt/esmvaltool_disk2/work/<username>

.. warning::

  conda environments should not be created in the home directory because it is on a very small disk,
  but rather in a directory with your username under `/mnt/esmvaltool_disk2/work/<username>`

Next, we need to set up the input files

The ``/work`` partition is visible from the VM so you can run the compare tool straight on the VM.

The steps to running the compare tool on the VM are the following:

- run date: log the run date here
- conda env: log the name of the conda environment you are using
- ESMValTool branch: log the name of the code branch you are using (e.g. `v2.8.x`)
- prerequisite - install `imagehash`: `pip install imagehash`
- reference run (v2.7.0; previous stable release): `export reference_dir=/work/bd0854/b382109/v270` (contains `preproc/` dirs too, 122 recipes)
- current run (v2.8.0): `export current_dir=path_to_current_run`
- run the :ref:`comparison script<compare_recipe_runs>` with:

.. code-block:: bash

  nohup python ESMValTool/esmvaltool/utils/testing/regression/compare.py --reference $reference_dir --current $current_dir > compare_v280_output.txt

Copy the comparison txt file to the release issue.
Some of the recipes will appear as having identical output to the one from previous release.
However, others will need human inspection.
Ask the recipe maintainers (`@ESMValGroup/esmvaltool-recipe-maintainers`_) and ESMValTool Development Team (`@ESMValGroup/esmvaltool-developmentteam`_) to provide assistance in checking the results.
Here are some guidelines on how to perform the human inspection:

- look at plots from current run vs previous release run: most of them will be identical, but if Matplotlib
  has changed some plotting feature, images may look slightly different so the comparison script may report them
  if the difference is larger than the threshold - but Mark I eyeball inspection will show they are identical
- other plots will differ due to changes in plot settings (different colours, axes etc) due to updated settings from the
  diagnostic developers: if they look similar enough, then it's fine
- report (and subsequently open issues) if you notice major differences in plots; most times a simple comment on the
  release issue, whereby you tag the diagnostic developers leads to them having a look at the plots and OK-ing them; if that's
  not the case, then open a separate issue. You can example of release issues containing overview lists and tables
  of failures and problems in `2881 <https://github.com/ESMValGroup/ESMValTool/issues/2881>`_
  and `3076 <https://github.com/ESMValGroup/ESMValTool/issues/3076>`_.

Appendix
--------

Here you can find a list of utility scripts used to run recipes and analyse the results:

- :ref:`Python scripts<utils_batch_jobs>` that create slurm submission scripts and parse slurm log files.
- :ref:`Python script<compare_recipe_runs>` that compares one or more recipe runs to known good previous run(s).
- :ref:`Python script<overview_page>`  that creates the ``index.html`` and ``debug.html`` overview pages.

.. _`@ESMValGroup/esmvaltool-recipe-maintainers`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-recipe-maintainers
.. _`@ESMValGroup/esmvaltool-developmentteam`: https://github.com/orgs/ESMValGroup/teams/esmvaltool-developmentteam
