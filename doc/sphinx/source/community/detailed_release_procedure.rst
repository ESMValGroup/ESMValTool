.. _detailed-release-procedure:

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


Where the work is done e.g. work done on DKRZ/Levante:

- submit files in `/home/b/b382109/submit` (for v2.7.0)
- output in `/scratch/b/b382109/esmvaltool_output` (for v2.7.0)

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
`conda env export > ToolEnv270Test.txt`.

Modifications to configuration files need to be documented as well.
To test recipes, it is recommended to only use the default DKRZ data directories, simply by uncommenting
the DKRZ-Levante block of the default ``config-user.yml`` file.

Submit run scripts - test recipe runs
-------------------------------------

We are now ready to start running all the available recipes, to compare output against previous release. Running is currently done
via batch scripts submitted to a schedulers (SLURM). Generate the submission scripts using the ``generate.py`` Python script;
you can find a copy of the script either in `/home/b/b382109/Tool_Release_270_Scripts` or
in the draft Pull Request https://github.com/ESMValGroup/ESMValTool/pull/2883.

You will have to set the name of your environment, your email address (if you want to get email notifications for successful/failed jobs) and the name of the directory you want to store the job outputs. The name of the account is the same (`bk1088`), and the default partition is set to `compute`.
More information on running jobs with SLURM on DKRZ/Levante can be found `here
<https://docs.dkrz.de/doc/levante/running-jobs/index.html>`_.

.. note::

  It remains to be checked see if the project name `bk1088` will still
  be active after the termination of IS-ENES3

Some recipes need other job requirements, you can add their headers in the `SPECIAL_RECIPES` dictionary. Otherwise the header will be written following the template that is written in the lines below. If you want to exclude recipes, you can do so by uncommenting the exclude lines.

The launch scripts will be saved in the same directory you execute the script from. You can find the ones used for the v2.7.0 release in `/home/b/b382109/submit`

To submit the scripts use the `sbatch` submit scripts (that make use of the SLURM scheduler) produced by `generate.py`,
and execute them like any other shell script. You can check the status of your BATCH queue by invoking:

.. code-block:: bash

  squeue -u $USER

Also, for computationally-heavy recipes, you can require more memory and/or time, see e.g. edited batch header below
(note the `compute` partition which is used for such heavy runs):

.. code-block:: bash

  #SBATCH --partition=compute
  #SBATCH --time=08:00:00
  #SBATCH --constraint=512G

.. note::

  On DKRZ/Levante, a user can't have more than 20 SLURM jobs running at a time.
  As soon as a job is finished, the next one should start. More information on the job handling at DKRZ `here
  <https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html#levante-partitions-and-limits>`_.

Once all jobs are completed, assemble some statistics so that issues with certain recipes
can be followed-up, and document this information in the release issue, such as:

- number of successfully run recipes
- number of failed recipes with preprocessor errors (can they be fixed? Can the fixes be included in the release?)
- number of failed recipes with diagnostic errors (can they be fixed? Can the fixes be included in the release?)
- number of recipes that are missing data
- number of recipes that have various other issues (and document them)

To parse the output of all these runs use the `parse_recipes_output.py` Python script, included at the
same locations where the generation script is.

Running the comparison
----------------------

To compare the newly produced output from running all recipes, follow these steps below.

Login and access to the DKRZ esmvaltool virtual machine (VM) - results from recipe runs
are stored on the VM; login with:

.. code-block:: bash

  ssh user@esmvaltool.dkrz.de

where `user` is your DKRZ/Levante user name; then get and install miniconda on the VM, and
if you have a Miniconda installer already downloaded in your Levante $HOME

.. code-block:: bash

  scp Miniconda3-py39_4.12.0-Linux-x86_64.sh user@esmvaltool.dkrz.de:~

.. warning::

  conda environments should not be created in the home directory because it is on a very small disk,
  but rather in a directory with your username under `/mnt/esmvaltool_disk2/work/`

Next, we need to set up the input files

.. note::

  If you wrote recipe runs output to Levante's `/scratch` partition, be aware that
  the data will be removed after two weeks, so you will have to move the output data
  to the VM, using a process that's not killed by a logoff e.g. a `nohup` job. Also add
  the link to the issue/discussion so anyone can see the results.
  This makes it much easier to ask for feedback from recipe developers.

  .. code-block:: bash

    nohup cp -r /scratch/b/$USER/esmvaltool_output/* /shared/esmvaltool/v2.x.x/

  where `bd0854/b382109` is the project location in `work`


The `/work` partition is visible by the VM so you can run the compare tool straight on the VM.

Do not store final release results on the VM including `/preproc/` dirs, the total
size for all the recipes output, including `/preproc/` dirs is in the 4.5TB ballpark,
much too high for the VM storage capacity! Therefore we would recommend using the option
to remove preprocessing directories upon recipe running successfully `--remove-preproc-dir=True`
at runtime, or set `remove_preproc_dir: true` in the configuration file.

The steps to running the compare tool on the VM are the following:

- run date: log the run date here
- conda env: log the name of the conda environment you are using
- ESMValTool branch: log the name of the code branch you are using (e.g. `v2.8.x`)
- prerequisite - install `imagehash`: `pip install imagehash`
- reference run (v2.7.0): `export reference_dir=/work/bd0854/b382109/v270` (contains `preproc/` dirs too, 122 recipes)
- current run (v2.8.0): `export current_dir=path_to_current_run`
- command to run: `nohup python ESMValTool/esmvaltool/utils/testing/regression/compare.py --reference $reference_dir --current $current_dir > compare_v280_output.txt`

Some of the recipes will appear as having identical output to the one from previous release. However, others
will need human inspection. Ask the recipe maintainers (`@ESMValGroup/esmvaltool-recipe-maintainers`_) and ESMValTool Development Team (`@ESMValGroup/esmvaltool-developmentteam`_) to provide assistance in checking the results.
Here are some guidelines on how to perform the human inspection:

- look at plots from current run vs previous release run: most of them will be identical, but if Matplotlib
  has changed some plotting feature, images will have slightly different metadata so the comparison script will report them
  as different - but Mark I eyeball inspection will show they are identical
- other plots will differ due to changes in plot settings (different colours, axes etc) due to updated settings from the
  diagnostic developers: if they look similar enough, then it's fine
- report (and subsequently open issues) if you notice major differences in plots; most times a simple comment on the
  release issue, whereby you tag the diagnostic developers leads to them having a look at the plots and OK-ing them; if that's
  not the case, then open a separate issue

Appendix
--------

Here you can find a list of useful files and directories:

- Formatted list of current recipes (as of v2.7.0) to be used with Markdown entries (on DKRZ/Lvante) at `/home/b/b382109/Tool_Release_270_Scripts/all_recipes.md` or in the draft Pull Request https://github.com/ESMValGroup/ESMValTool/pull/2883
- last release (v2.7.0) submit scripts on DKRZ/Levante `/home/b/b382109/submit`
- Miniconda3 installer file on DKRZ/Levante `/home/b/b382109/Miniconda3-py39_4.12.0-Linux-x86_64.sh` (remember to immediately update conda after using it, it is fairly old, from May 2022)
- list of Autoassess reference files and masks on DKRZ/Levante `/home/b/b382109/autoassess_files`
