.. _recipes_seaice_extents_sh:

Sea Ice area and extents
========================

Overview
--------

This recipe plots sea ice concentration from CICE (sea ice model) output 
around the southern polar region and compares it to the NSIDC CDR 
(National Snow and Ice Data Centre, Climate Data Record) dataset.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_seaice_extents_sh.yml

Diagnostics are stored in esmvaltool/diag_scripts/seaice_area_extents/

    * seaicearea_trends.py: plot minima and maxima sea ice area trends
    * seaice_mapextents.py: plot sea ice extent and differences to Observations


User settings in recipe
-----------------------

#. Script seaice_mapextents.py

   *Required settings for script*

   * `months`: months by month number which the mean are to be plotted


Variables
---------

* siconc (seaIce, monthly, longitude latitude time)
* areacello (fx)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* NSIDC CDR sh (siconc - esmvaltool/cmorizers/data/formatters/datasets/nsidc_g02202_sh.py)


References
----------

* COSIMA(Consortium for Ocean-Sea Ice Modelling in Australia) recipe: https://cosima-recipes.readthedocs.io/en/latest/DocumentedExamples/SeaIce_Obs_Model_Compare.html

Example plots
-------------

.. _trends:
.. figure::  /recipes/figures/seaice_extents_sh/min_trend.png
   :align:   center

   Minima trends of sea ice area with observation data. ACCESS OM model data years from 0, added 1652 years to model years for comparability.

.. _map extents:
.. figure::  /recipes/figures/seaice_extents_sh/map_difference.png
   :align:   center

   Difference and extents of models with observations for selected months.