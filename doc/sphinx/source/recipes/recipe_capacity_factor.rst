.. _yml_capacity_factor:

Capacity factor of wind power: Ratio of average estimated power to theoretical maximum power
============================================================================================

Overview
--------

The goal of this diagnostic is to compute the wind capacity factor,  taking as input the daily instantaneous surface wind speed, which is then extrapolated to obtain the  wind speed at a height of 100 m as described in Lledó (2017). 

The capacity factor is a normalized indicator of the suitability of wind speed conditions to produce electricity, irrespective of the size and number of installed turbines. This indicator is provided for three different classes of wind turbines (IEC, 2005) that are designed specifically for low, medium and high wind speed conditions. 

The user can select the region, temporal range and season of interest. 

The output of the recipe is a netcdf file containing the capacity factor for each of the three turbine classes.
.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_capacity_factor.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* capacity_factor.R: calculates the capacity factor for the three turbine classes.
* PC.r: calculates the power curves for the three turbine classes.


User settings
-------------

User setting files are stored in recipes/

#. recipe_capacity_factor.yml

   *Required settings for script*

   * power_curves: (should not be changed)

Variables
---------

* sfcWind (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* IEC. (2005). International Standard IEC 61400-1, third edition, International Electrotechnical Commission. https://webstore.iec.ch/preview/info_iec61400-1%7Bed3.0%7Den.pdf

* Lledó, L. (2017). Computing capacity factor. Technical note BSC-ESS-2017-001, Barcelona Supercomputing Center. Available online at https://earth.bsc.es/wiki/lib/exe/fetch.php?media=library:external:bsc-ess-2017-001-c4e_capacity_factor.pdf [last accessed 11 October 2018]

Example plots
-------------

.. _fig_capfactor1:
.. figure::  /recipes/figures/capacity_factor/capacity_factor_IPSL-CM5A-LR_1980-2005.png
   :align:   center
   :width:   14cm

