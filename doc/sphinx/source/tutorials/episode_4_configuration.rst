.. _episode_4_configuration:

***************
Configuration
***************

.. admonition:: Overview
   :class: note

   .. grid:: 2
      :gutter: 1
      :margin: 3 3 0 5

      .. grid-item-card:: Timings
         :columns: 12

         * Teaching: 10 min
         * Exercises: 15 min

      .. grid-item-card:: Questions

         *  What is the user configuration file and how should I use it?

      .. grid-item-card:: Objectives

         * Understand how ESMValTool is configured
         * Prepare a personalized ESMValTool configuration
         * Configure ESMValTool to use stored climate data and to download climate data

      .. grid-item-card:: Compatibility
         :columns: 12

         ESMValTool v2.15.0


Configuring ESMValTool via YAML files
======================================

ESMValTool provides a set of predefined configuration files.
These include the files specifying the default configuration values,
but also machine-specific files that include data sources for various
HPC systems.

To show all available files, run

.. code-block:: bash

   esmvaltool config list

All configuration files are `YAML files <https://yaml.org/spec/1.2/spec.html>`_.

To customize your configuration via YAML files, you can copy one of the existing files.
For example, to copy the file containing the default values for many options, run

.. code-block:: bash

   esmvaltool config copy defaults/config-user.yml

The default configuration file will be downloaded to the default location:
``~/.config/esmvaltool/config-user.yml``, where ``~`` is the
path to your home directory. Note that files and directories starting with a
period are "hidden", to see the ``.config`` directory in the terminal use
``ls -la ~``. Note, if a configuration file by that name already exists in the
default location, the ``config copy`` command will not update the file as
ESMValTool will not overwrite the file. You will have to move the file first if
you want an updated copy of the default user configuration file.

We run a text editor called ``nano`` to have a look inside the configuration file
and then modify it if needed:

.. code-block:: bash

   nano ~/.config/esmvaltool/config-user.yml

If ``nano`` does not work on your system, or if you prefer a different editor,
any other editor can be used, e.g. ``vim``.

This file contains the information for:

- Output settings
- Destination directory
- Auxiliary data directory
- Number of tasks that can be run in parallel
- ...

.. admonition:: Text editor side note
   :class: discussion

   No matter what editor you use, you will need to know where it searches for
   and saves files. If you start it from the shell, it will (probably) use your
   current working directory as its default location. We use ``nano`` in examples
   here because it is one of the least complex text editors.
   Press :kbd:`ctrl` + :kbd:`O` to save the file,
   and then :kbd:`ctrl` + :kbd:`X` to exit ``nano``.


Destination directory
=====================

The example configuration file contains the option ``output_dir``, which is the
rootpath where ESMValTool will store its output folders containing e.g. figures,
data, logs, etc. With every run, ESMValTool automatically generates a new output
folder determined by recipe name, and date and time using the format:
YYYYMMDD_HHMMSS.

.. admonition:: Set the destination directory
   :class: admonition-todo

   Let's name our destination directory ``esmvaltool_output`` in the current
   directory. ESMValTool should write the output to this path, so make sure you
   have the disk space to write output to this directory.
   How do we set this in the `config-user.yml`?

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      We use ``output_dir`` entry in the ``config-user.yml`` file as:

      .. code-block:: yaml

         output_dir: ./esmvaltool_output

      If the ``esmvaltool_output`` does not exist, ESMValTool will generate it for you.


Output settings
===============

Additionally you can configure the output settings that
inform ESMValTool about your preference for output.
Most of these settings are fairly self-explanatory.

.. admonition:: Saving preprocessed data
   :class: admonition-todo

   Later in this tutorial, we will want to look at the contents of the
   ``preproc`` folder.
   This folder contains preprocessed data and is removed by default when
   ESMValTool is run.
   In the configuration, which settings can be modified to prevent
   this from happening?

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      If the option ``remove_preproc_dir`` is set to ``false``, then the
      ``preproc/`` directory contains all the pre-processed data and the metadata
      interface files.
      If the option ``save_intermediary_cubes`` is set to ``true`` then data will
      also be saved after each preprocessor step in the folder ``preproc``. Note
      that saving all intermediate results to file will result in a considerable
      slowdown, and can quickly fill your disk.


Other settings
==============

.. admonition:: Auxiliary data directory
   :class: note

   The ``auxiliary_data_dir`` setting is the path where any required additional
   auxiliary data files are stored. This location allows us to tell the diagnostic
   script where to find the files if they can not be downloaded at runtime. This
   option should not be used for model or observational datasets, but for data
   files (e.g. shape files) used in plotting such as coastline descriptions and
   if you want to feed some additional data (e.g. shape files) to your recipe.

   .. code-block:: yaml

      auxiliary_data_dir: ~/auxiliary_data

   See more information `here <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html?highlight=auxiliary_data#top-level-configuration-options>`_
   in the ESMValTool documentation.


.. admonition:: Number of parallel tasks
   :class: note

   This option enables you to perform parallel processing. You can choose the
   number of tasks in parallel as 1/2/3/4/... or you can set it to ``null``. That
   tells ESMValTool to use the maximum number of available CPUs. For the purpose
   of the tutorial, please set ESMValTool use only 1 cpu:

   .. code-block:: yaml

      max_parallel_tasks: 1

   In general, if you run out of memory, try setting ``max_parallel_tasks`` to 1.
   Then, check the amount of memory you need for that by inspecting the file
   ``run/resource_usage.txt`` in the output directory. Using the number there you
   can increase the number of parallel tasks again to a reasonable number for the
   amount of memory available in your system.


Customizing your configuration
===============================

By default, configuration files are read from the directory ``~/.config/esmvaltool``.
This can be changed via the ``ESMVALTOOL_CONFIG_DIR`` environment variable.
In addition another custom configuration directory can be specified via the
``--config_dir`` command line argument.
We will learn how to do this in the next lesson on recipes.

It is possible to have several configuration files with different purposes, for
example: ``dask_options.yml``, ``data_sources.yml``.
In this case, ESMValTool searches for all YAML files within each of the
configuration directories and merges them together. How this is done is explained
in the `Configuration section <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#yaml-files>`_
in the ESMValTool documentation.

To show the final configuration that is actually used when running ESMValTool, you can use

.. code-block:: bash

   esmvaltool config show


Rootpath to input data
======================

ESMValTool uses several categories (in ESMValTool, these are referred to as projects)
for input data based on their source (e.g.
CMIP6, CMIP5, obs4mips, OBS6, OBS). For example, CMIP is used for a dataset from
the Climate Model Intercomparison Project whereas OBS may be
used for an observational dataset.
More information about the projects used in ESMValTool is available in the
`documentation <https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/find_data.html>`_.
The ``data`` section for each project in the configuration
files defines sources of input data. The easiest way to get started with these is to
copy one of the example configuration files and tailor it to your needs.

When using ESMValTool on your own machine, the recommended setup can
be obtained by running the command

.. code-block:: bash

   esmvaltool config copy data-local-esmvaltool.yml

After the file ``data-local-esmvaltool.yml`` has been copied to your configuration
directory ``~/.config/esmvaltool/``, you can update the ``rootpath`` and the
``dirname_template`` to match your file locations. The ``rootpath`` specifies the
directories where ESMValTool will look for input data of the specific project. The
``dirname_template`` setting describes the file structure for each project.

If you are working on a HPC system, there are also several configurations
for popular HPC systems available that you can use instead, e.g. JASMIN, DKRZ,
ETH, and IPSL. To list the available example files, run the command:

.. code-block:: bash

   esmvaltool config list data-hpc

To load the configuration suitable for the HPC system at DKRZ, run:

.. code-block:: bash

   esmvaltool config copy data-hpc-dkrz.yml

It is also possible to ask ESMValTool to download climate model data as needed. When
running ESMValTool you can automatically download the files required to run a recipe
from ESGF for the projects CMIP3, CMIP5, CMIP6, CORDEX, and obs4MIPs. For this,
copy the appropriate configuration file by running

.. code-block:: bash

   esmvaltool config copy data-intake-esgf.yml

Additionally, it is necessary to configure
`intake-esgf <https://intake-esgf.readthedocs.io/en/latest/configure>`_.
For this you need to copy the file ``conf.yaml`` (see below) into the directory ``~/.config/intake-esgf`` and
update the ``local_cache`` and ``esg_dataroot`` with your desired download directory in this intake-esgf
configuration file. The updated file should look like this:

.. dropdown:: conf.yaml
      :color: secondary
      :icon: eye

      .. code-block:: yaml

         additional_df_cols: []
         break_on_error: true
         confirm_download: false
         default_open_kwargs:
         chunks: auto
         compat: override
         coords: minimal
         data_vars: minimal
         default_project: CMIP7
         download_db: ~/.config/intake-esgf/download.db
         esg_dataroot:
         - /p/css03/esgf_publish
         - /eagle/projects/ESGF2/esg_dataroot
         - /global/cfs/projectdirs/m3522/cmip6/
         - /glade/campaign/collections/cmip.mirror
         globus_indices:
         ESGF2-US-1.5-Catalog: true
         local_cache:
         - ~/.esgf/
         logfile: ~/.config/intake-esgf/esgf.log
         num_threads: 6
         print_log_on_error: false
         requests_cache:
         cache_name: intake-esgf/requests-cache.sqlite
         expire_after: 3600
         use_cache_dir: true
         slow_download_threshold: 0.5
         solr_indices:
         esg-dn1.nsc.liu.se: false
         esgf-data.dkrz.de: false
         esgf-node.ipsl.upmc.fr: false
         esgf.ceda.ac.uk: false
         esgf.nci.org.au: false
         stac_indices:
         discovery.east.esgf.io: true
         discovery.west.esgf.io: true


.. admonition:: Set the correct rootpaths
   :class: admonition-todo

   In this tutorial, we will work with data from
   `CMIP5 <https://esgf-node.llnl.gov/projects/cmip5/>`_
   and `CMIP6 <https://esgf-node.llnl.gov/projects/cmip6>`_.
   How can we modify the ``rootpath`` to make sure the data path is set correctly
   for both CMIP5 and CMIP6?

   Note: to get the data, check the instructions in :doc:`tutorial_setup`.

   .. dropdown:: Solution
      :color: secondary
      :icon: eye

      - Are you working on your own local machine?

        You need to copy ``data-local-esmvaltool.yml`` into your configuration directory
        and specify the root path of the folder where the data is available (e.g.,
        [YOUR_CLIMATE_DATA_DIR]) as:

        .. code-block:: yaml

           projects:
           ...
             CMIP6:
               data:
                 local:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project}/{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{version}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}_{grid}*.nc"
             CMIP5:
               data:
                 local:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project.lower}/{product}/{institute}/{dataset}/{exp}/{frequency}/{modeling_realm}/{mip}/{ensemble}/{version}/{short_name}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}*.nc"

      - Are you working on your local machine and you want to download missing data using ESMValTool?

        You need to configure ``intake-esgf`` (see above) and add the root path of the folder where the data
        has been downloaded to in ``data-local-esmvaltool.yml`` as specified in the ``esgf-cache``.

        .. code-block:: yaml

           projects:
           ...
             CMIP6:
               data:
                 local:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project}/{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{version}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}_{grid}*.nc"
                 esgf-cache:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project}/{activity}/{institute}/{dataset}/{exp}/{ensemble}/{mip}/{short_name}/{grid}/{version}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}_{grid}*.nc"
             CMIP5:
               data:
                 local:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project.lower}/{product}/{institute}/{dataset}/{exp}/{frequency}/{modeling_realm}/{mip}/{ensemble}/{version}/{short_name}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}*.nc"
                 esgf-cache:
                   type: esmvalcore.io.local.LocalDataSource
                   rootpath: [YOUR_CLIMATE_DATA_DIR]
                   dirname_template: "{project.lower}/{product}/{institute}/{dataset}/{exp}/{frequency}/{modeling_realm}/{mip}/{ensemble}/{version}"
                   filename_template: "{short_name}_{mip}_{dataset}_{exp}_{ensemble}*.nc"

      - Are you working on a computer cluster like Jasmin or DKRZ?

        Site-specific path to the data for JASMIN/DKRZ/ETH/IPSL
        are already available in specific configuration files. You
        need to copy this file in your configuration directory.
        For example, on DKRZ, run:

        .. code-block:: bash

           esmvaltool config copy data-hpc-dkrz.yml

      - For more information about configuring the data sources, see also the
        `relevant section <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#project-specific-configuration>`_
        in the ESMValTool documentation.


.. admonition:: Configuration via command line
   :class: note

   In addition, all configuration options can also be specified via the command
   line and those settings will overwrite any setting given by the YAML files. You
   can find more information in the
   `ESMValTool documentation <https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#command-line-arguments>`_.


.. admonition:: Key points
   :class: important

   * ESMValTool can be configured through YAML files located in ``~/.config/esmvaltool`` or command line arguments
   * The final configuration is created by merging the contents of all YAML files and command line arguments
   * Users can choose to use one big configuration file, or spread its contents among many small configuration files
   * ESMValTool can be configured to automatically download climate data from ESGF
