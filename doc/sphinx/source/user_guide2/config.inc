.. _config:

*******************
Configuration files
*******************

There are several configuration files in ESMValTool:

  - config-user.yml
  - config-developer.yml
  - config-references.yml
  - environment.yml

This section describes those files.

User config file
================

The ``config-user.yml`` configuration file contains all the global level
information needed by ESMValTool. The following shows the default settings from
the ``config-user.yml`` file.

.. code-block:: yml

    # Diagnostics create plots? [true]/false
    write_plots: true
    # Diagnositcs write NetCDF files? [true]/false
    write_netcdf: true
    # Set the console log level debug, [info], warning, error
    log_level: info
    # verbosity is deprecated and will be removed in the future
    # verbosity: 1
    # Exit on warning? true/[false]
    exit_on_warning: false
    # Plot file format? [ps]/pdf/png/eps/epsi
    output_file_type: pdf
    # Destination directory
    output_dir: ./esmvaltool_output
    # Auxiliary data directory (used for some additional datasets)
    auxiliary_data_dir: ./auxiliary_data
    # Use netCDF compression true/[false]
    compress_netcdf: false
    # Save intermediary cubes in the preprocessor true/[false]
    save_intermediary_cubes: false
    # Remove the preproc dir if all fine
    remove_preproc_dir: true
    # Run at most this many tasks in parallel null/[1]/2/3/4/..
    # Set to null to use the number of available CPUs.
    # Make sure your system has enough memory for the specified number of tasks.
    max_parallel_tasks: 1
    # Path to custom config-developer file, to customise project configurations.
    # See config-developer.yml for an example. Set to None to use the default
    config_developer_file: null
    # Get profiling information for diagnostics
    # Only available for Python diagnostics
    profile_diagnostic: false

    # Rootpaths to the data from different projects (lists are also possible)
    rootpath:
      CMIP5: [~/cmip5_inputpath1, ~/cmip5_inputpath2]
      OBS: ~/obs_inputpath
      default: ~/default_inputpath

    # Directory structure for input data: [default]/BADC/DKRZ/ETHZ/etc
    # See config-developer.yml for definitions.
    drs:
      CMIP5: default

Most of these settings are fairly self-explanatory, ie:

.. code-block:: yml

    # Diagnostics create plots? [true]/false
    write_plots: true
    # Diagnositcs write NetCDF files? [true]/false
    write_netcdf: true

The ``write_plots`` setting is used to inform ESMValTool about your preference
for saving figures. Similarly, the ``write_netcdf`` setting is a boolean which
turns on or off the writing of netCDF files.

.. code-block:: yml

    # Auxiliary data directory (used for some additional datasets)
    auxiliary_data_dir: ./auxiliary_data

The ``auxiliary_data_dir`` setting is the path to place any required
additional auxiliary data files. This method was necessary because certain
Python toolkits such as cartopy will attempt to download data files at run
time, typically geographic data files such as coastlines or land surface maps.
This can fail if the machine does not have access to the wider internet. This
location allows us to tell cartopy (and other similar tools) where to find the
files if they can not be downloaded at runtime. To reiterate, this setting is
not for model or observational datasets, rather it is for data files used in
plotting such as coastline descriptions and so on.


Tip: You choose your config.yml file at run time, so you could have several
available with different purposes. One for formalised run, one for debugging,
etc...

Developper config file
================================

This configuration file describes the file system structure for several
key projects (CMIP5, CMIP6) on several key machines (BADC, CP4CDS, DKRZ, ETHZ,
SMHI, BSC).

For instance, the code structure of the CMIP5 project is set up differently
at each site.

.. code-block:: yml

  CMIP5:
    input_dir:
      default: '/'
      BADC: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/latest/[short_name]'
      CP4CDS: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/[short_name]/latest/'
      DKRZ: '[institute]/[dataset]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/[latestversion]/[short_name]'
      ETHZ: '[exp]/[mip]/[short_name]/[dataset]/[ensemble]/'
      SMHI: '[dataset]/[ensemble]/[exp]/[frequency]'
      BSC: '[project]/[exp]/[dataset.lower]'

So, for instance, a typical file path on BADC would be:

.. code-block:: yml

    MOHC/HadGEM2-CC/rcp85/mon/ocean/Omon/r1i1p1/latest


To be more explicit

.. code-block:: yml

    [institute]/[dataset ]/[exp]/[frequency]/[modeling_realm]/[mip]/[ensemble]/latest/[short_name]
    MOHC       /HadGEM2-CC/rcp85/mon        /ocean           /Omon /r1i1p1    /latest/tos

.. warning::
    Section to be completed

References config file
================================

.. warning::
    Section to be completed

Environment config file
================================

This is the conda envronment which informs conda as to which packages
and which versions are needed in order to install ESMValTool.
It is called in the update and install commands:

.. code-block:: bash

    conda env install --name esmvaltool --file ESMValTool/environment.yml
    conda env update --name esmvaltool --file ESMValTool/environment.yml
