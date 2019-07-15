.. _recipes_quicklook:

Clouds
======

Overview
--------

The recipe recipe_quicklook.yml computes the global and zonal
timeseries of variables. If the ESMValTool is running in the quicklook mode
the concatinated files will be plotted.
The diagnostics read in lat-time fields.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_quicklook.yml

Diagnostics are stored in diag_scripts/quicklooks/

    * global_timeseries.py: global timeseries
    * zonal_timeseries.py: zonal timeseries


User settings in recipe
-----------------------

#. Script global_timeseries.py

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * time_int: min and max for time axis
   * y_min: min of y axis
   * y_max: max of y axis
   * multimodel_plot: if True: additional plot with all datasets
                 qicklook mode: all concatinated files
                 no quicklook mode: all dataset given in recipe

#. Script zonal_timeseries.py

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * time_int: min and max for time axis
   * lat_int: min and max latitude values
   * val_levs: values for contour levels


Example plots
-------------

.. _fig_global:
.. figure::  /recipes/figures/quicklooks/Model_HadGEM2-CC_tas_global_timeseries.png
   :align:   center

   Timeseries of global mean of variable tas for dataset HadGEM2-CC.

.. _fig_global_multi:
.. figure::  /recipes/figures/quicklooks/MultiModel_tas_global_timeseries.png
   :align:   center

   Timeseries of global mean of variable tas.

.. _fig_zonal:
.. figure::  /recipes/figures/quicklooks/Model_HadGEM2-CC_tas_zonal_timeseries.png
   :align:   center

   Timeseries of zonal mean of variable tas for dataset HadGEM2-CC.

