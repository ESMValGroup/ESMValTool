.. _recipes_eady_growth_rate:

Eady growth rate
================

Overview
--------

This recipe computes the maximum Eady Growth Rate and performs the annual and seasonal means, storing 
the results for each dataset. 
For the seasonal means, the results are plotted over the North-Atlantic region for the selected
pressure levels.


Available recipes and diagnostics
---------------------------------

Recipes are stored in ``esmvaltool/recipes/``

    * ``recipe_eady_growth_rate.yml``

Diagnostics are stored in ``esmvaltool/diag_scripts/eady_growth_rate/``

    * ``eady_growth_rate.py``: Computes and stores the eady growth rate. 
      Plots can be produced for the seasonal mean over the North Atlantic region.


User settings in recipe
-----------------------

#. Script ``eady_growth_rate.py``

   *Required settings for script*

   * ``time_statistic``: Set to `'annual'` to compute the annual mean. Set to `'seasonal'` to compute the seasonal mean.

   *Optional settings for script*

   * ``plot_levels``: list of pressure levels to be plotted for the seasonal mean. If not specified, all levels will be plotted.


Variables
---------

* ta (atmos, monthly mean, longitude latitude level time)
* zg (atmos, monthly mean, longitude latitude level time)
* ua (atmos, monthly mean, longitude latitude level time) 

References
----------

Brian J Hoskins and Paul J Valdes. On the existence of storm-tracks. Journal of the atmospheric sciences, 47(15):1854–1864, 1990.

Example plots
-------------

.. _fig_eady_growth_rate:
.. figure::  /recipes/figures/eady_growth_rate/HadGEM3-GC31-LM_winter_eady_growth_rate_70000.png 
   :align:   center

   Eady Growth Rate values over the North-Atlantic region at 70000 Pa.
