.. _episode_5_recipe:

*************************************
Running your first recipe
*************************************

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: 15 min
         * Exercises: 15 min

      .. grid-item-card:: Questions

         * How to run a recipe?
         * What happens when I run a recipe?

      .. grid-item-card:: Objectives

         * Run an existing ESMValTool recipe
         * Examine the log information
         * Navigate the output created by ESMValTool
         * Make small adjustments to an existing recipe

      .. grid-item-card:: Compatibility
         :columns: 12

         ESMValTool v2.15.0


Running an existing recipe
==========================

This episode describes how ESMValTool recipes work, how to run a recipe and how
to explore the recipe output. By the end of this episode, you should be able to
run your first recipe, look at the recipe output, and make small modifications.

The recipe format has briefly been introduced in the
:doc:`episode_1_introduction` episode. To see all the recipes that are shipped
with ESMValTool, type:

.. code-block:: bash

   esmvaltool recipes list

We will start by running
`examples/recipe_python.yml <https://docs.esmvaltool.org/en/latest/recipes/recipe_examples.html>`_.

.. code-block:: bash

   esmvaltool run examples/recipe_python.yml

or, if you have the user configuration file in your current directory, then:

.. code-block:: bash

   esmvaltool run --config_dir . examples/recipe_python.yml

If everything is okay, you should see a lot of output from ESMValTool on the
command line. The final message should be ``Run was successful``. The exact output
varies depending on your machine, but it should look similar to the example log
output below.

.. dropdown:: Output
   :color: secondary
   :icon: eye

   .. code-block:: text

      2026-10-15 12:04:36,964 UTC [194240] INFO
      ______________________________________________________________________
                _____ ____  __  ____     __    _ _____           _
              | ____/ ___||  \/  \ \   / /_ _| |_   _|__   ___ | |
              |  _| \___ \| |\/| |\ \ / / _` | | | |/ _ \ / _ \| |
              | |___ ___) | |  | | \ V / (_| | | | | (_) | (_) | |
              |_____|____/|_|  |_|  \_/ \__,_|_| |_|\___/ \___/|_|
      ______________________________________________________________________

      Earth System Model Evaluation Tool

      A community tool for the evaluation of Earth system models.

      https://esmvaltool.org

      The Earth System Model Evaluation Tool (ESMValTool) is a community
      diagnostics and performance metrics tool for the evaluation of Earth
      System Models (ESMs) that allows for routine comparison of single or
      multiple models, either against predecessor versions or against
      observations.

      Tutorial: https://tutorial.esmvaltool.org
      Documentation: https://docs.esmvaltool.org
      Contact: esmvaltool-dev@listserv.dfn.de

      If you find this software useful for your research, please cite it using
      https://doi.org/10.5281/zenodo.3387139 for ESMValCore or
      https://doi.org/10.5281/zenodo.3401363 for ESMValTool or
      any of the reference papers listed at https://esmvaltool.org/references/.

      Have fun!

      2026-10-15 12:04:36,966 UTC [194240] INFO    Package versions
      2026-10-15 12:04:36,966 UTC [194240] INFO    ----------------
      2026-10-15 12:04:36,966 UTC [194240] INFO    ESMValCore: 2.15.0
      2026-10-15 12:04:36,967 UTC [194240] INFO    ESMValTool: 2.16.0.dev75+g6afe583bf.d20261015
      2026-10-15 12:04:36,967 UTC [194240] INFO    ----------------
      2026-10-15 12:04:36,967 UTC [194240] INFO    Reading configuration files from:
      /home/b/USERNAME/ESMValTool/ESMValTool/.pixi/envs/default/lib/python3.14/site-packages/esmvalcore/config/configurations/defaults (defaults)
      /home/b/USERNAME/.config/esmvaltool (default user configuration directory)
      2026-10-15 12:04:36,967 UTC [194240] INFO    Writing program log files to:
      /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/main_log.txt
      /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/main_log_debug.txt
      /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/cmor_log.txt
      2026-10-15 12:04:37,593 UTC [194240] INFO    Starting the Earth System Model Evaluation Tool at time: 2026-10-15 12:04:37 UTC
      2026-10-15 12:04:37,593 UTC [194240] INFO    ----------------------------------------------------------------------
      2026-10-15 12:04:37,593 UTC [194240] INFO    RECIPE   = /home/b/USERNAME/ESMValTool/ESMValTool/esmvaltool/recipes/examples/recipe_python.yml
      2026-10-15 12:04:37,593 UTC [194240] INFO    RUNDIR     = /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run
      2026-10-15 12:04:37,593 UTC [194240] INFO    WORKDIR    = /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/work
      2026-10-15 12:04:37,593 UTC [194240] INFO    PREPROCDIR = /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/preproc
      2026-10-15 12:04:37,594 UTC [194240] INFO    PLOTDIR    = /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/plots
      2026-10-15 12:04:37,594 UTC [194240] INFO    ----------------------------------------------------------------------
      2026-10-15 12:04:37,665 UTC [194240] INFO    Running tasks using at most 256 processes
      2026-10-15 12:04:37,665 UTC [194240] INFO    If your system hangs during execution, it may not have enough memory for keeping this number of tasks in memory.
      2026-10-15 12:04:37,665 UTC [194240] INFO    If you experience memory problems, try reducing 'max_parallel_tasks' in your configuration.
      2026-10-15 12:04:51,651 UTC [194240] INFO    Creating tasks from recipe
      2026-10-15 12:04:51,651 UTC [194240] INFO    Creating tasks for diagnostic map
      2026-10-15 12:04:51,651 UTC [194240] INFO    Creating diagnostic task map/script1
      2026-10-15 12:04:51,652 UTC [194240] INFO    Creating preprocessor task map/tas
      2026-10-15 12:04:51,652 UTC [194240] INFO    Creating preprocessor 'to_degrees_c' task for variable 'tas'
      2026-10-15 12:04:51,788 UTC [194240] INFO    Found input files for Dataset: tas, Amon, CMIP6, BCC-ESM1, historical, r1i1p1f1, mon, gn, v20181214
      2026-10-15 12:04:51,917 UTC [194240] INFO    Found input files for Dataset: tas, atmos, CMIP7, CanESM6-0-MR, historical, r1i1p1f1, tavg-h2m-hxy-u, mon, glb, g150, v20260629
      2026-10-15 12:04:51,918 UTC [194240] INFO    PreprocessingTask map/tas created.
      2026-10-15 12:04:51,918 UTC [194240] INFO    Creating tasks for diagnostic timeseries
      2026-10-15 12:04:51,918 UTC [194240] INFO    Creating diagnostic task timeseries/script1
      2026-10-15 12:04:51,919 UTC [194240] INFO    Creating preprocessor task timeseries/tas_amsterdam
      2026-10-15 12:04:51,919 UTC [194240] INFO    Creating preprocessor 'annual_mean_amsterdam' task for variable 'tas_amsterdam'
      2026-10-15 12:04:52,052 UTC [194240] INFO    Found input files for Dataset: tas, Amon, CMIP6, BCC-ESM1, historical, r1i1p1f1, mon, gn, v20181214
      2026-10-15 12:04:52,191 UTC [194240] INFO    Found input files for Dataset: tas, atmos, CMIP7, CanESM6-0-MR, historical, r1i1p1f1, tavg-h2m-hxy-u, mon, glb, g150, v20260629
      2026-10-15 12:04:52,193 UTC [194240] INFO    PreprocessingTask timeseries/tas_amsterdam created.
      2026-10-15 12:04:52,193 UTC [194240] INFO    Creating preprocessor task timeseries/tas_global
      2026-10-15 12:04:52,193 UTC [194240] INFO    Creating preprocessor 'annual_mean_global' task for variable 'tas_global'
      2026-10-15 12:04:52,435 UTC [194240] INFO    Found input files for Dataset: tas, Amon, CMIP6, BCC-ESM1, historical, r1i1p1f1, mon, gn, v20181214, supplementaries: areacella, fx, 1pctCO2, fx, v20190613
      2026-10-15 12:04:52,723 UTC [194240] INFO    Found input files for Dataset: tas, atmos, CMIP7, CanESM6-0-MR, historical, r1i1p1f1, tavg-h2m-hxy-u, mon, glb, g150, v20260629, supplementaries: areacella, ti-u-hxy-u, fx
      2026-10-15 12:04:52,723 UTC [194240] INFO    PreprocessingTask timeseries/tas_global created.
      2026-10-15 12:04:52,724 UTC [194240] INFO    These tasks will be executed: map/tas, timeseries/script1, timeseries/tas_global, map/script1, timeseries/tas_amsterdam
      2026-10-15 12:04:52,739 UTC [194240] INFO    Wrote recipe with version numbers and wildcards to:
      file:///scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/recipe_python_filled.yml
      2026-10-15 12:04:52,740 UTC [194240] INFO    Downloading missing data (this may take a while...). Details can be found in the debug log at /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/main_log_debug.txt
      2026-10-15 12:04:53,002 UTC [194240] INFO    Successfully downloaded missing data
      2026-10-15 12:04:54,490 UTC [194240] INFO    Using Dask distributed scheduler (address: tcp://127.0.0.1:33777, dashboard link: http://127.0.0.1:8787/status)
      2026-10-15 12:04:54,490 UTC [194240] INFO    Running 5 tasks using 5 processes
      2026-10-15 12:04:54,764 UTC [194240] INFO    Progress: 3 tasks running, 2 tasks waiting for ancestors, 0/5 done
      2026-10-15 12:05:01,704 - distributed.nanny - WARNING - Restarting worker (status=Status.running)
      WARNING:esmvalcore.cmor._fixes.fix.genericfix:Long name changed from 'Grid-Cell Area for Atmospheric Variables' to 'Grid-Cell Area for Atmospheric Grid Variables'
      (for file /work/PROJECT/DATA/ESMValTool2/download/CMIP6/CMIP/BCC/BCC-ESM1/1pctCO2/r1i1p1f1/fx/areacella/gn/v20190613/areacella_fx_BCC-ESM1_1pctCO2_r1i1p1f1_gn.nc)
      2026-10-15 12:05:05,461 UTC [200348] INFO    map/tas                       [########################################] | 100% Completed | 2.01 s
      2026-10-15 12:05:05,468 UTC [200348] INFO    Successfully completed task map/tas (priority 1) in 0:00:05.678506
      2026-10-15 12:05:05,874 UTC [194240] INFO    Progress: 2 tasks running, 2 tasks waiting for ancestors, 1/5 done
      2026-10-15 12:05:05,975 UTC [194240] INFO    Progress: 3 tasks running, 1 tasks waiting for ancestors, 1/5 done
      2026-10-15 12:05:08,657 UTC [200349] INFO    timeseries/tas_amsterdam      [########################################] | 100% Completed | 3.17 s
      2026-10-15 12:05:08,680 UTC [200349] INFO    Successfully completed task timeseries/tas_amsterdam (priority 3) in 0:00:08.859076
      2026-10-15 12:05:09,095 UTC [194240] INFO    Progress: 2 tasks running, 1 tasks waiting for ancestors, 2/5 done
      2026-10-15 12:05:10,831 UTC [200350] INFO    timeseries/tas_global         [########################################] | 100% Completed | 2.15 s
      2026-10-15 12:05:10,844 UTC [200350] INFO    Successfully completed task timeseries/tas_global (priority 4) in 0:00:10.712513
      2026-10-15 12:05:11,304 UTC [194240] INFO    Progress: 1 tasks running, 1 tasks waiting for ancestors, 3/5 done
      2026-10-15 12:05:11,404 UTC [194240] INFO    Progress: 2 tasks running, 0 tasks waiting for ancestors, 3/5 done
      2026-10-15 12:05:16,508 UTC [194240] INFO    Progress: 1 tasks running, 0 tasks waiting for ancestors, 4/5 done
      2026-10-15 12:05:22,715 UTC [194240] INFO    Progress: 0 tasks running, 0 tasks waiting for ancestors, 5/5 done
      2026-10-15 12:05:22,716 UTC [194240] INFO    Successfully completed all tasks.
      2026-10-15 12:05:23,465 UTC [194240] INFO    Wrote recipe with version numbers and wildcards to:
      file:///scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/run/recipe_python_filled.yml
      2026-10-15 12:05:23,737 UTC [194240] INFO    Wrote recipe output to:
      file:///scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436/index.html
      2026-10-15 12:05:23,738 UTC [194240] INFO    It looks like you are connected to a remote machine via SSH. To show the output html file, you can try the following command on your local machine:
      server=136.172.124.7 && port=31415 && ssh -t -L ${port}:localhost:${port} USERNAME@${server} /home/b/USERNAME/ESMValTool/ESMValTool/.pixi/envs/default/bin/python -m http.server ${port} -d /scratch/b/USERNAME/esmvaltool_output/recipe_python_20261015_120436
      Then visit http://localhost:31415 in your browser
      2026-10-15 12:05:23,738 UTC [194240] INFO    If the port 31415 is already in use, you can replace it with any other free one (e.g., 12789). If you are connected through a jump host, replace the server IP address 136.172.124.7 with your SSH server name
      2026-10-15 12:05:23,738 UTC [194240] INFO    Ending the Earth System Model Evaluation Tool at time: 2026-10-15 12:05:23 UTC
      2026-10-15 12:05:23,738 UTC [194240] INFO    Time for running the recipe was: 0:00:46.145520
      2026-10-15 12:05:23,984 UTC [194240] INFO    Maximum memory used (estimate): 11.2 GB
      2026-10-15 12:05:23,985 UTC [194240] INFO    Sampled every second. It may be inaccurate if short but high spikes in memory consumption occur.
      2026-10-15 12:05:23,988 UTC [194240] INFO    Run was successful

.. admonition:: Pro tip: ESMValTool search paths
   :class: note

   You might wonder how ESMValTool was able to find the recipe file, even though
   it is not in your working directory. All the recipe paths printed from
   ``esmvaltool recipes list`` are relative to ESMValTool's installation location.
   This is where ESMValTool will look if it cannot find the file by following the
   path from your current working directory.


Investigating the log messages
==============================

Let's dissect what's happening here.

.. admonition:: Output files and directories
   :class: admonition-todo

   After the banner and general information, the output starts with some important
   locations.

   1. Did ESMValTool use the right config files?
   2. What is the path to the example recipe?
   3. What is the main output folder generated by ESMValTool?
   4. Can you guess what the different output directories are for?
   5. ESMValTool creates two log files. What is the difference?

   .. dropdown:: Answers
      :color: secondary
      :icon: eye

      1. The list of directories involves the directory with all default config files
         (for example ``../ESMValCore/esmvalcore/config/configurations/defaults``)
         and the default user configuration directory with our edited files from the
         previous episode (default: ``~/.config/esmvaltool``).
      2. ESMValTool found the recipe in its installation directory, for example
         ``/home/users/username/mambaforge/envs/esmvaltool/esmvaltool/recipes/examples/``
         or on a server, e.g.
         ``/apps/jasmin/community/esmvaltool/ESMValTool_<version>/esmvaltool/recipes/examples/recipe_python.yml``.
      3. ESMValTool creates a time-stamped output directory for every run. In this
         case, it should be something like ``recipe_python_YYYYMMDD_HHMMSS``. This
         folder is created inside the output directory specified in the previous
         episode: ``~/esmvaltool_tutorial/esmvaltool_output``.
      4. There should be four output folders:

         * ``run/``: stores general information about the run, such as log messages
           and a copy of the recipe file.
         * ``work/``: stores output files that are not figures.
         * ``preproc/``: stores pre-processed data.
         * ``plots/``: stores output figures.

      5. The log files are:

         * ``main_log.txt``: a copy of the command-line output.
         * ``main_log_debug.txt``: more detailed information for debugging.

.. admonition:: Debugging: No ``preproc`` directory?
   :class: note

   If you are missing the ``preproc`` directory, then your ``config-user.yml`` file
   has the value ``remove_preproc_dir`` set to ``true`` (this is used to save disk
   space). Please set this value to ``false`` and run the recipe again.

After the output locations, there are two main sections that can be distinguished
in the log messages:

* Creating tasks
* Executing tasks

.. admonition:: Analyse the tasks
   :class: admonition-todo

   List all the tasks that ESMValTool is executing for this recipe. Can you guess
   what this recipe does?

   .. dropdown:: Answer
      :color: secondary
      :icon: eye

      Just after all the "creating tasks" lines and before "executing tasks", we find
      the following output:

      .. code-block:: text

         [3098040] INFO    These tasks will be executed: timeseries/tas_global,
         map/tas, map/script1, timeseries/tas_amsterdam, timeseries/script1

      So there are three tasks related to time series: global temperature,
      Amsterdam temperature, and a script (``tas``: near-surface air temperature).
      There are also two tasks related to a map: something with temperature and
      another script.


Examining the recipe file
=========================

To get more insight into what is happening, we will have a look at the recipe
file itself. Use the following command to copy the recipe to your working
directory:

.. code-block:: bash

   esmvaltool recipes get examples/recipe_python.yml

Now you should see the recipe file in your working directory (type ``ls`` to
verify). Use ``nano`` to open it:

.. code-block:: bash

   nano recipe_python.yml

For reference, you can also view the recipe in the box below.

.. dropdown:: ``recipe_python.yml``
   :color: secondary
   :icon: eye

   .. code-block:: yaml

      # ESMValTool
      # recipe_python.yml
      #
      # See https://docs.esmvaltool.org/en/latest/recipes/recipe_examples.html
      # for a description of this recipe.
      #
      # See https://docs.esmvaltool.org/projects/esmvalcore/en/latest/recipe/overview.html
      # for a description of the recipe format.
      ---
      documentation:
        description: |
          Example recipe that plots a map and timeseries of temperature.

        title: Recipe that runs an example diagnostic written in Python.

        authors:
          - andela_bouwe
          - righi_mattia

        maintainer:
          - schlund_manuel

        references:
          - acknow_project

        projects:
          - esmval
          - c3s-magic

      datasets:
        - {dataset: BCC-ESM1, project: CMIP6, mip: Amon, exp: historical, ensemble: r1i1p1f1, grid: gn}
        - {dataset: CanESM6-0-MR, institute: CCCma, project: CMIP7, mip: atmos, branding_suffix: tavg-h2m-hxy-u, region: glb, ensemble: r1i1p1f1, grid: g150}

      preprocessors:
        # See https://docs.esmvaltool.org/projects/esmvalcore/en/latest/recipe/preprocessor.html
        # for a description of the preprocessor functions.

        to_degrees_c:
          convert_units:
            units: degrees_C

        annual_mean_amsterdam:
          extract_location:
            location: Amsterdam
            scheme: linear
          annual_statistics:
            operator: mean
          multi_model_statistics:
            statistics:
              - mean
            span: overlap
          convert_units:
            units: degrees_C

        annual_mean_global:
          area_statistics:
            operator: mean
          annual_statistics:
            operator: mean
          convert_units:
            units: degrees_C

      diagnostics:

        map:
          description: Global map of temperature in January 2000.
          themes:
            - phys
          realms:
            - atmos
          variables:
            tas:
              preprocessor: to_degrees_c
              exp: historical
              frequency: mon
              timerange: 2000/P1M
              caption: |
                Global map of {long_name} in January 2000 according to {dataset}.
          scripts:
            script1:
              script: examples/diagnostic.py
              quickplot:
                plot_type: pcolormesh
                cmap: RdYlBu_r

        timeseries:
          description: Annual mean temperature in Amsterdam and global mean since 1850.
          themes:
            - phys
          realms:
            - atmos
          variables:
            tas_amsterdam:
              short_name: tas
              exp: historical
              frequency: mon
              preprocessor: annual_mean_amsterdam
              timerange: 1850/2000
              caption: Annual mean {long_name} in Amsterdam according to {dataset}.
            tas_global:
              short_name: tas
              exp: historical
              frequency: mon
              preprocessor: annual_mean_global
              timerange: 1850/2000
              caption: Annual global mean {long_name} according to {dataset}.
          scripts:
            script1:
              script: examples/diagnostic.py
              quickplot:
                plot_type: plot

Do you recognize the basic recipe structure that was introduced in the
:doc:`episode_1_introduction` episode?

* ``documentation`` with relevant citation information
* ``datasets`` that should be analysed
* ``preprocessors`` with groups of common preprocessing steps
* ``diagnostics`` scripts performing more specific evaluation steps

.. admonition:: Analyse the recipe
   :class: admonition-todo

   Try to answer the following questions:

   1. Who wrote this recipe?
   2. Who should be approached if there is a problem with this recipe?
   3. How many datasets are analysed?
   4. What does the preprocessor called ``annual_mean_global`` do?
   5. Which script is applied for the diagnostic called ``map``?
   6. Can you link specific lines in the recipe to the tasks that we saw before?
   7. How is the location of the city specified?
   8. How is the temporal range of the data specified?

   .. dropdown:: Answers
      :color: secondary
      :icon: eye

      1. The example recipe is written by Bouwe Andela and Mattia Righi.
      2. Manuel Schlund is listed as the maintainer of this recipe.
      3. Two datasets are analysed:

         * CMIP6 data from the model ``BCC-ESM1``
         * CMIP5 data from the model ``bcc-csm1-1``

      4. The preprocessor ``annual_mean_global`` computes an area mean as well as
         annual means.
      5. The diagnostic called ``map`` executes a script referred to as ``script1``.
         This is a Python script named ``examples/diagnostic.py``.
      6. There are two diagnostics: ``map`` and ``timeseries``. Under the diagnostic
         ``map`` we find two tasks:

         * a preprocessor task called ``tas`` applying the preprocessor ``to_degrees_c``
           to the variable ``tas``;
         * a diagnostic task called ``script1`` applying the script
           ``examples/diagnostic.py`` to the preprocessed data (``map/tas``).

         Under the diagnostic ``timeseries`` we find three tasks:

         * a preprocessor task called ``tas_amsterdam`` applying the preprocessor
           ``annual_mean_amsterdam`` to the variable ``tas``;
         * a preprocessor task called ``tas_global`` applying the preprocessor
           ``annual_mean_global`` to the variable ``tas``;
         * a diagnostic task called ``script1`` applying the script
           ``examples/diagnostic.py`` to the preprocessed data (``timeseries/tas_global``
           and ``timeseries/tas_amsterdam``).

      7. The ``extract_location`` preprocessor is used to get data for a specific
         location here. ESMValTool interpolates to the location based on the chosen
         scheme.
      8. The ``timerange`` tag is used to extract data from a specific time period.
         The start time is ``01/01/2000`` and the span of time to calculate means is
         ``1 Month`` given by ``P1M``.

.. admonition:: Pro tip: short names and variable groups
   :class: note

   The preprocessor tasks in ESMValTool are called ``variable groups``. For the
   diagnostic ``timeseries``, we have two variable groups: ``tas_amsterdam`` and
   ``tas_global``. Both operate on the variable ``tas`` (as indicated by the
   ``short_name``), but they apply different preprocessors. For the diagnostic
   ``map``, the variable group itself is named ``tas``, and you will notice that
   we do not explicitly provide the ``short_name``. This is a shorthand built into
   ESMValTool.

.. admonition:: Output files
   :class: admonition-todo

   Have another look at the output directory created by the ESMValTool run.
   Which files and folders are created by each task?

   .. dropdown:: Answer
      :color: secondary
      :icon: eye

      * ``map/tas`` creates ``/preproc/map/tas``, which contains preprocessed data
        for each input dataset, a ``metadata.yml`` file describing the contents of
        these datasets and provenance information in the form of ``.xml`` files.
      * ``timeseries/tas_global`` creates ``/preproc/timeseries/tas_global``, which
        contains preprocessed data for each input dataset, a ``metadata.yml`` file,
        and provenance information in the form of ``.xml`` files.
      * ``timeseries/tas_amsterdam`` creates ``/preproc/timeseries/tas_amsterdam``,
        which contains preprocessed data for each input dataset, plus a combined
        ``MultiModelMean``, a ``metadata.yml`` file and provenance files.
      * ``map/script1`` creates ``/run/map/script1`` with general information and a
        log of the diagnostic script run. It also creates ``/plots/map/script1`` and
        ``/work/map/script1``, which contain output figures and output datasets,
        respectively. For each output file, there is corresponding provenance
        information in the form of ``.xml``, ``.bibtex`` and ``.txt`` files.
      * ``timeseries/script1`` creates ``/run/timeseries/script1`` with general
        information and a log of the diagnostic script run. It also creates
        ``/plots/timeseries/script1`` and ``/work/timeseries/script1``, which contain
        output figures and output datasets, respectively. For each output file,
        there is corresponding provenance information in the form of ``.xml``,
        ``.bibtex`` and ``.txt`` files.

.. admonition:: Pro tip: diagnostic logs
   :class: note

   When you run ESMValTool, any log messages from the diagnostic script are not
   printed on the terminal. But they are written to the ``log.txt`` files in the
   folder ``/run/[DIAG_NAME]/log.txt``.

   ESMValTool does print a command that can be used to re-run a diagnostic script.
   When you use this, the output will be printed to the command line.


Modifying the example recipe
============================

Let's make a small modification to the example recipe. Notice that, now that you
have copied and edited the recipe, you can use

.. code-block:: bash

   esmvaltool run recipe_python.yml

to refer to your local file rather than the default version shipped with
ESMValTool.

.. admonition:: Change your location
   :class: admonition-todo

   Modify and run the recipe to analyse the temperature for your own location.

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      In principle, you only have to modify the location in the preprocessor called
      ``annual_mean_amsterdam``. However, it is good practice to also replace all
      instances of ``amsterdam`` with the correct name of your location. Otherwise
      the log messages and output will be confusing. You are free to modify the
      names of preprocessors or diagnostics.

      In the diff below you will see the changes we have made to the file. The top
      two lines are the filenames and the lines like ``@@ -39,9 +39,9 @@`` represent
      the line numbers in the original and modified file, respectively.

      .. code-block:: diff

         --- recipe_python.yml
         +++ recipe_python_london.yml
         @@ -39,9 +39,9 @@
              convert_units:
                units: degrees_C

         -  annual_mean_amsterdam:
         +  annual_mean_london:
              extract_location:
         -      location: Amsterdam
         +      location: London
                scheme: linear
              annual_statistics:
                operator: mean
         @@ -83,7 +83,7 @@
                cmap: Reds

           timeseries:
         -    description: Annual mean temperature in Amsterdam and global mean since 1850.
         +    description: Annual mean temperature in London and global mean since 1850.
             themes:
               - phys
             realms:
         @@ -92,9 +92,9 @@
               tas_amsterdam:
                 short_name: tas
                 mip: Amon
         -        preprocessor: annual_mean_amsterdam
         +        preprocessor: annual_mean_london
                 timerange: 1850/2000
         -        caption: Annual mean {long_name} in Amsterdam according to {dataset}.
         +        caption: Annual mean {long_name} in London according to {dataset}.
               tas_global:
                 short_name: tas
                 mip: Amon


Congratulations!
----------------

You have run your first ESMValTool recipe and learned how to inspect the relevant
log messages and output directories. You can now start making small changes to an
existing recipe to adapt it to your own use case.

.. admonition:: Key points
   :class: important

   * ESMValTool recipes work out of the box if the input data are available.
   * There are strong links between the recipe, the log file, and the output folders.
   * Recipes can easily be modified to re-use existing code for your own use case.
