.. _recipes_sea_surface_salinity:

Sea Surface Salinity Evaluation
===============================

Overview
--------

This recipe compares the regional means of sea surface salinity with a
reference dataset (ESACCI-SEA-SURFACE-SALINITY v1 or v2 by default).
To do this, the recipe generate plots for the timeseries of each region and
a radar plot showing the correlation of dataset and reference timeseries for
each region during the time they both exists.

The recipe is created in a way that should make possible (although is not
tested) to use it for other variables and datasets, even for more that one at
a time.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_sea_surface_salinity

Diagnostics are stored in diag_scripts/sea_surface_salinity/

    * compare_salinity.py: plot timeseries for each region and generate radar
      plot.


User settings in recipe
-----------------------

#. compare_salinity.py

   *Required settings for script*

   none

   *Optional settings for script*

   none

   *Required settings for variables*

   * ref_model: name of reference data set

   *Optional settings for variables*

   none

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*


Variables
---------

* sos (ocean, monthly mean, time depth_id)


Observations and reformat scripts
---------------------------------

* ESACCI-SEA-SURFACE-SALINITY (sos)


References
----------

* Please, contact authors


Example plots
-------------

.. figure:: /recipes/figures/sea_surface_salinity/radar.png
   :align: center

   Radar plot showing correlation of average sea surface salinity for multiple
   regions with the observations
