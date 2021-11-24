.. _recipe_monitor:

Monitor
#######

Overview
========

Available recipes and diagnostics
=================================

Recipes are stored in `recipes/``

  - recipe_monitor.yml

Diagnostics are stored in `diag_scripts/monitor/`

  - monitor.py:
    Plots preprocessor output directly from the preprocessor.
  - compute_eofs.py:
    A sample on how to use the monitor structure to show other metrics.
    Computes and plots the map of the first EOF and the associated PC timeseries.

User settings
=============

User setting files are stored in recipes/ and in a dedicted yaml config file

recipe_monitor.yml
------------------

  * plots:
    a dictionary containing the plots to make, with its own options.
  * cartopy_data_dir:
    Path to cartopy data dir. Defaults to None.
    See https://scitools.org.uk/cartopy/docs/latest/cartopy.html
  * plot_folder:
    Path to the folder to store the figures. It is defined as the
    input paths in ``config-developer.yml``. See
    https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#input-file-paths
    for more details. Defaults to ``~/plots/{dataset}/{exp}/{modeling_realm}/{real_name}``
  * plot_filename:
    Filename pattern for the plots. it is defined as the input
    files in in ``config-developer.yml``. See
    https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/configure.html#input-file-paths
    for more details. Defaults to ``{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}``
  * config_file:
    Path to the monitor config file. Defaults to
    ``monitor_config.yml`` in the same folder as the diagnostic script.

Plot specific options:
^^^^^^^^^^^^^^^^^^^^^^

- monclim:
   + maps:
     List of maps to plot, as defined in the config file. Defaults to ``[global]``.
   + months:
     Select only specific months. Defaults to ``None`` (do not select).
   + plot_size:
     Size of ech individual figure. Default's to ``(5, 4)``.
   + columns:
     Number of columns in the plot. Defaults to ``3``.
   + rows:
     Number of rows in the plot. Defaults to ``4``.
- seasonclim:
   + maps:
     List of maps to plot, as defined in the config file. Defaults to ``[global]``.
- clim:
   + maps:
     List of maps to plot, as defined in the config file. Defaults to ``[global]``.
- annual_cycle: No options
- timeseries: No options.

monitor_config.yml
------------------

A yaml file containing map and variable specific options.

Contains two dictionaries, ``maps`` and ``variables``.

Each entry in ``maps`` correspond to a map definitions. See below a sample with
comments to define each option

.. code-block:: yaml

   maps:
      global: # Map name, choose a meaningful one
         projection: PlateCarree # Cartopy projection to use
         projection_kwargs: # Dictionary with Cartopy's projection kwargs.
            central_longitude: 285
         smooth: true # If true, interpolate values to get smoother maps. If not, all points in a cells will get the exact same color
         lon: [-120, -60, 0, 60, 120, 180] # Set longitude ticks
         lat: [-90, -60, -30, 0, 30, 60, 90] # Set latitude ticks
         colorbar_location: bottom
         extent: null # If defined, restrict the projection to a region. Format [lon1, lon2, lat1, lat2]
         suptitle_pos: 0.87 # Title position in the figure.

Each entry in ``variable`` correspond to a variable definitions.
Use the default entry to apply generic options to al variables.
See below a sample with comments to define each option

.. code-block:: yaml

   variables:
      # Define default. Variable definitions completely override the default
      # not just the values defined. If you want to override only the defined
      # values, use yaml anchors as shown
      default: &default
         colors: RdYlBu_r # Matplotlib colormap to use for the co,orbar
         N: 20 # Number of map intervals to plot
         bad: [0.9, 0.9, 0.9] # Color to use when no data
      pr:
         <<: *default
         colors: gist_earth_r
         # Define bounds of the colorbar, as a list of
         bounds: 0-10.5,0.5 # Set colorbar bounds, as a list or in the format min-max,interval
         extend: max # Set extend parameter of mpl colorbar. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html
      sos:
         # If default is defined, entries are treated as map specific option.
         # Missing values in map definitionas are taken from variable's default
         # definition
         default:
            <<: *default
            bounds: 25-41,1
            extend: both
         arctic:
            bounds: 25-40,1
         antarctic:
            bounds: 30-40,0.5
      nao: &nao
         <<: *default
         extend: both
         # Variable definitions can override map parameters. Use with caution.
         bounds: [-0.03, -0.025, -0.02, -0.015, -0.01, -0.005, 0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
         projection: PlateCarree
         smooth: true
         lon: [-90, -60, -30, 0, 30]
         lat: [20, 40, 60, 80]
         colorbar_location: bottom
         suptitle_pos: 0.87
      sam:
         <<: *nao
         lat: [-90, -80, -70, -60, -50]
         projection: SouthPolarStereo
         projection_kwargs:
            central_longitude: 270
         smooth: true
         lon: [-120, -60, 0, 60, 120, 180]

Variables
=========

* Any, but the dimensionality should match the expected by each plot


Observations and reformat scripts
=================================

*None*

Example plots
=============

.. _fig_climglobal:
.. figure::  /recipes/figures/monitor/clim.png
   :align:   center
   :width:   14cm

Global climatology of tas.

.. _fig_seasonclimglobal:
.. figure::  /recipes/figures/monitor/seasonclim.png
   :align:   center
   :width:   14cm

Seasonal climatology of pr, with a custom colorbar

.. _fig_monthlyclimglobal:
.. figure::  /recipes/figures/monitor/monclim.png
   :align:   center
   :width:   14cm

Monthly climatology of sivol, only for March and September.

.. _fig_timeseries:
.. figure::  /recipes/figures/monitor/timeseries.png
   :align:   center
   :width:   14cm

Timeseries of Niño 3.4 index, computed directly with the preprocessor.

.. _fig_annual_cycle:
.. figure::  /recipes/figures/monitor/annualcycle.png
   :align:   center
   :width:   14cm

Annual cycle of tas.
