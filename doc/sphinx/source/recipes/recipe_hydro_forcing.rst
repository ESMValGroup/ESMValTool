.. _recipes_hydro_forcing:

Hydro forcing comparison
========================

Overview
--------

This recipe can be used to assess the agreement between forcing datasets
(i.e. MSWEP, ERA5, ERA-Interim) for a certain catchment. The recipe can be used
to:

1. Plot a timeseries of the raw daily data
2. Plot monthly aggregrated data over a certain period
3. Plot the monthly climate statistics over a certain period


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * ``recipe_hydro_forcing.yml``

Diagnostics are stored in esmvaltool/diag_scripts/hydrology/

    * ``hydro_forcing.py``: Compares and plots precipitation for MSWEP / ERA5 / ERA-5 interim datasets


User settings in recipe
-----------------------

All hydrological recipes require a shapefile as an input to produce forcing data. This shapefile determines the shape of the basin for which the data will be cut out and processed. All recipes are tested with `the shapefiles <https://github.com/eWaterCycle/recipes_auxiliary_datasets/tree/master/>`_  that are used for the eWaterCycle project. In principle any shapefile can be used, for example, the freely available basin shapefiles from the `HydroSHEDS project <https://www.hydrosheds.org/>`_.

#. recipe ``hydrology/hydro_forcing.yml``

  *Optional preprocessor settings:*

    * ``extract_shape``: The region specified here should match the catchment

  *Required settings for script:*

    * ``entry_point``: Define which diagnostic function to run. Choices:

      * ``sample_year``: Plot a timeseries of the raw daily data
      * ``total_precipitation``: Plot monthly aggregrated data over a certain period
      * ``climatology``: Plot the monthly climate statistics over a certain period

Variables
---------

* pr (atmos, daily or monthly, longitude, latitude, time)


Observations
------------

All data can be used directly without any preprocessing.

*  ERA-Interim
*  ERA5
*  MSWEP

.. References
.. ----------

.. * xxx

Example plots
-------------

.. _fig_hydro_forcing_1:
.. figure::  /recipes/figures/hydrology/daily_precipitation_plot.png
  :align:   center

  Daily precipitation for 2015-01-01:2016-12-31.

.. _fig_hydro_forcing_2:
.. figure::  /recipes/figures/hydrology/monthly_total_precipitation_plot.png
  :align:   center

  Monthly total precipitation for 2015-01:2016-12.

.. _fig_hydro_forcing_3:
.. figure::  /recipes/figures/hydrology/precipitation_per_month_plot.png
  :align:   center

  Precipitation per month for 2015:2016.
