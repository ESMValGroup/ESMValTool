.. _config-user:

***********************
Configuring ESMValTool
***********************

The ``config-user.yml`` configuration file contains all the global level
information needed by ESMValTool. The configuration is passed to ESMValTool
as a command line argument (see :ref:`Running ESMValTool <running>`).

An example configuration file can be downloaded `here <https://github.com/ESMValGroup/ESMValTool/blob/version2_development/config-user-example.yml>`_
and tailored for your system using the explanation below.

Most of these settings are fairly self-explanatory, ie:

.. code-block:: yaml

  # Diagnostics create plots? [true]/false
  write_plots: true
  # Diagnositcs write NetCDF files? [true]/false
  write_netcdf: true

The ``write_plots`` setting is used to inform ESMValTool about your preference
for saving figures. Similarly, the ``write_netcdf`` setting is a boolean which
turns on or off the writing of netCDF files.

The ```rootpath`` specifies the directories where ESMValTool will look for input
data. Similarly, ``output_dir`` specifies where ESMValTool will store its
output, i.e. figures, data, logs, etc. Make sure to set appropriate paths.

.. code-block:: yaml

  # Auxiliary data directory (used for some additional datasets)
  auxiliary_data_dir: ~/auxiliary_data

The ``auxiliary_data_dir`` setting is the path to place any required
additional auxiliary data files. This method was necessary because certain
Python toolkits such as cartopy will attempt to download data files at run
time, typically geographic data files such as coastlines or land surface maps.
This can fail if the machine does not have access to the wider internet. This
location allows us to tell cartopy (and other similar tools) where to find the
files if they can not be downloaded at runtime. To reiterate, this setting is
not for model or observational datasets, rather it is for data files used in
plotting such as coastline descriptions and so on.


.. note::

   The ``config-user.yml`` file is specified as argument at run time, so it is
   possible to have several available with different purposes: one for
   formalised runs, one for debugging, etc...
