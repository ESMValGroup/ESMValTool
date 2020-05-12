.. _recipes_capacity_factor:

Capacity factor of wind power: Ratio of average estimated power to theoretical maximum power
============================================================================================

Overview
--------

The goal of this diagnostic is to compute the wind capacity factor,  taking as input the daily instantaneous surface wind speed, which is then extrapolated to obtain the  wind speed at a height of 100 m as described in Lledó (2017).

The capacity factor is a normalized indicator of the suitability of wind speed conditions to produce electricity, irrespective of the size and number of installed turbines. This indicator is provided for three different classes of wind turbines (IEC, 2005) that are designed specifically for low, medium and high wind speed conditions.

The user can select the region, temporal range and season of interest.

The output of the recipe is a netcdf file containing the capacity factor for each of the three turbine classes.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_capacity_factor.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* capacity_factor.R: calculates the capacity factor for the three turbine classes.
* PC.R: calculates the power curves for the three turbine classes.


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

Main features of the selected turbines:

=================  ==================  ================  ==================  =================  ===================
Turbine name       Rotor diameter (m)  Rated power (MW)  Cut-in speed (m/s)  Rated speed (m/s)  Cut-out speed (m/s)

-----------------  ------------------  ----------------  ------------------  -----------------  -------------------
Enercon E70 2.3MW  70                  2.3               2.0                 16.0               25.0
Gamesa G80 2.0MW   80                  2.0               4.0                 17.0               25.0
Gamesa G87 2.0MW   87                  2.0               4.0                 16.0               25.0
Vestas V100 2.0MW  100                 2.0               3.0                 15.0               20.0
Vestas V110 2.0MW  110                 2.0               3.0                 11.5               20.0
=================  ==================  ================  ==================  =================  ===================

References
----------

* IEC. (2005). International Standard IEC 61400-1, third edition, International Electrotechnical Commission. https://webstore.iec.ch/preview/info_iec61400-1%7Bed3.0%7Den.pdf

* Lledó, L. (2017). Computing capacity factor. Technical note BSC-ESS-2017-001, Barcelona Supercomputing Center. Available online at https://earth.bsc.es/wiki/lib/exe/fetch.php?media=library:external:bsc-ess-2017-001-c4e_capacity_factor.pdf [last accessed 11 October 2018]

Example plots
-------------

.. _fig_capfactor1:
.. figure::  /recipes/figures/capacity_factor/capacity_factor_IPSL-CM5A-MR_2021-2050.png
   :align:   center
   :width:   14cm

Wind capacity factor for five turbines: Enercon E70 (top-left), Gamesa G80 (middle-top), Gamesa G87 (top-right), Vestas V100 (bottom-left) and Vestas V110 (middle-bottom) using the IPSL-CM5A-MR simulations for the r1p1i1 ensemble for the rcp8.5 scenario during the period 2021-2050.
