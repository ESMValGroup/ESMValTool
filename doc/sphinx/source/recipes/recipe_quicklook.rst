.. _recipes_quicklook:

Quicklooks
==========

Overview
--------

The recipe ``recipe_quicklook.yml`` computes the global and zonal timeseries of
variables. If the ESMValTool is running in the quicklook mode the concatenated
files will be plotted. The diagnostics can read any fields including the
necessary coordinates ('time' and in the case of zonal means 'latitude') and
automatically computes means over all other dimensions.


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

   * ``multi_dataset_plot``, *bool*, optional: Additionally plot all datasets.
     in one file.

   * ``time_range``, *list of float*, optional: Range for the ``time``
     coordinate in the plot.

   * ``y_range``, *list of float*, optional: Range for the ``y`` coordinate
     in the plot (= variable which is considered).

#. Script zonal_timeseries.py

   * ``latitude_range``, *list of float*, optional: Range for the ``latitude``
     coordinate in the plot.

   * ``levels``, *list of float*, optional: Values for contour levels.

   * ``time_range``, *list of float*, optional: Range for the ``time``
     coordinate in the plot.


Example plots
-------------

.. _fig_global:
.. figure::  /recipes/figures/quicklooks/tas_HadGEM2-CC_global_timeseries.png
   :align:   center

   Time series plot of global mean of variable tas for dataset HadGEM2-CC.

.. _fig_global_multi:
.. figure::  /recipes/figures/quicklooks/tas_global_timeseries.png
   :align:   center

   Time series plot of global mean of variable tas.

.. _fig_zonal:
.. figure::  /recipes/figures/quicklooks/tas_HadGEM2-CC_zonal_timeseries.png
   :align:   center

   Time series plot of zonal mean of variable tas for dataset HadGEM2-CC.

