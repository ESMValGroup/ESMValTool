.. _detailed-release-procedure:

Overview
========

The release procedure for ESMValTool is a fairly involved process (at the moment), so it
is important to be very well organized and to have documented each procedural steps, so that
the next release manager can follow said steps, and finalize the release without any delays.

The workflow below assumes an ESMValCore release candidate, or a completed stable release, have been released
and deployed on conda-forge and PyPI; it also assumes the release manager has access to accounts on LEVANTE [TODO insert info]

Procedural workflow
===================

Open an issue on GitHub
-----------------------

First, open an issue on GitHub where the release workflow is documented (see example https://github.com/ESMValGroup/ESMValTool/issues/2881).
Name it something relevant like "Recipe testing and comparison for release 2.x.x", and populate the isue description with information
about where the testing is taking place, what tools are used, and what versions, here is a template in Markdown:

Documenting system parameters
-----------------------------

Work done on DKRZ/Levante:

- submit files in `/home/b/b382109/submit`
- output in `/scratch/b/b382109/esmvaltool_output`

# System and Settings

We should document various versions so that the work can be reproduced in case there
is an issue, or release work needs to be picked up mid-release by another release manager:

- documentic `conda`/`mamba` versions:

.. code-block:: bash
  mamba --version

- documentic the `git` branch and its state:

.. code-block:: bash
  git status

Creating and documenting the runtime environment
------------------------------------------------

On Levante:

.. code-block:: bash
  mamba env create -n tool270Test -f environment.yml
  conda activate tool270Test

Make a copy of the environment file, and attach it in the release testing issue; to
record the environment in a yaml file use e.g. `conda env export > ToolEnv270Test.txt`.

Extraneous file movements
-------------------------

Mention any special data movement one needs to perfom to run the tests e.g.:
for v2.7.0 the release manager moved the autoassess-specific files to
`/home/b/$USER/autoassess_files` on DKRZ/Levante.

List modifications to configuration files
-----------------------------------------

Example - for v2.7.0 the release manager had to modify `config-user.yml` file by
adding a number of extra paths for DKRZ-specific data pools:

.. code-block:: yaml
  CMIP6:
    - /work/bd0854/DATA/ESMValTool2/CMIP6_DKRZ
    - /work/bd0854/DATA/ESMValTool2/download/CMIP6
  CMIP5:
    - /work/bd0854/DATA/ESMValTool2/CMIP5_DKRZ
    - /work/bd0854/b309141/additional_CMIP5
    - /work/bd0854/DATA/ESMValTool2/download/cmip5/output1
    - /work/bd0854/DATA/ESMValTool2/download/cmip5

Submit run scripts - test recipe runs
-------------------------------------

Submit the batch scripts that will run all recipes. Assemble some statistics so that issues with certain recipes
can be followed-up:

- number of successfully run recipes
- number of failed recipes with Diagnostic error (can they be fixed? Can the fixes be included in the release?)
- number of recipes that are missing data
- number of recipes that have various other issues (and document them)

To submit the scripts use the `sbatch` submit scripts (that make use of the SLURM scheduler)
and execute them like any other shell script. You can check the status of your BATCH queue by invoking:

.. code-block:: bash
  squeue -u b382109

Also, for computationally-heavy recipes, you can require more memory and/or time, see e.g. edited batch header below
(note the `compute` partition which is used for such heavy runs):

.. code-block:: bash
  #SBATCH --partition=compute
  #SBATCH --time=08:00:00
  #SBATCH --constraint=512G

.. note::
  On DKRZ/Levante, a user can't have more than 20 SLURM jobs running at a time.
  As soon as a job is finished, the next one should start

Running the comparison
======================

Login and access to the DKRZ esmvaltool VM - results from recipe runs are stored on the VM; login with:

.. code-block:: bash
  ssh user@esmvaltool.dkrz.de

where `user` is your DKRZ/Levante user name; then get and install miniconda on VM, and
if you already have a Miniconda installer already downloaded in your Levante $HOME

.. code-block:: bash
  scp Miniconda3-py39_4.12.0-Linux-x86_64.sh user@esmvaltool.dkrz.de:~

Next, we need to set up the input files

.. note::
  If you wrote recipe runs output to Levante's `/scratch` partition, be aware that
  the data will be removed after two weeks, so you will have to move the output data
  to the `/work` partition, via e.g. a `nohup` job:

  .. code-block:: bash
    nohup cp -r /scratch/b/$USER/esmvaltool_output/* /work/bd0854/b382109/v2xx

  where `bd0854/b382109` is the project location in `work`


The `/work` partition is visible by the VM so you can run the compare tool straight on the VM.

Do not store final release results on the VM including `/preproc/` dirs, the total
size for all the recipes output, including `/preproc/` dirs is in the 4.5TB ballpark,
much too high for the VM storage capacity!

The steps to running the compare tool at VM are the following:

- run date: log the run date here
- conda env: log the name of the conda environment you are using
- ESMValTool branch: log the name of the code branch you are using (e.g. `v2.8.x`)
- prerquisite - install `imagehash`: `pip install imagehash`
- reference run (v2.7.0): `export reference_dir=/work/bd0854/b382109/v270` (contains `preproc/` dirs too, 122 recipes)
- current run (v2.8.0): `export current_dir=path_to_current_run`
- command to run: `nohup python ESMValTool/esmvaltool/utils/testing/regression/compare.py reference_dir current_dir > compare_v280_output.txt`


Appendix
========

Formatted list of current recipes (as of v2.7.0) to be used with Markdown entries:


