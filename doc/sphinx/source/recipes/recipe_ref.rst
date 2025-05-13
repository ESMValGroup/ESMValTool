.. _recipes_REF:

CMIP Rapid Evaluation Framework (REF)
======================================

Overview
--------

Here ESMValTool recipes are collected which will be used in the CMIP 
`Rapid Evaluation Framework (REF) <https://wcrp-cmip.org/cmip7/rapid-evaluation-framework/>`__.


Available recipes 
-----------------

Recipes are stored in recipes/ref

* recipe_ref_cre.yml:
    Maps and zonal means of longwave and shortwave cloud radiative effect
* recipe_fire.yml
    Maps of burnt area fraction, fire weather control, and fuel load continuity control.


Example plots:
-----------------

.. _fig_ref_1:
.. figure::  /recipes/figures/ref/map_lwcre_MPI-ESM1-2-LR_Amon.png
   :align:   center
   
   Geographical map of the climatological mean longwave cloud radiative 
   effect from MPI-ESM1-2-LR and CERES-EBAF Ed4.2 and their difference.

.. _fig_ref_2:
.. figure::  /recipes/figures/ref/variable_vs_lat_lwcre_ambiguous_dataset_Amon.png
   :align:   center

   Zonal averages of the climatological mean longwave cloud radiative
   effect from CERES-EBAF Ed4.2 (solid black), ESACCI-CLOUD (dashed black), 
   ISCCP-FH (dotted black) and the MPI-ESM1-2-LR model (blue).

.. _fig_ref_3:
.. figure::  /recipes/figures/ref/burnt_fraction_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center
   
   Burnt area fraction for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.

.. _fig_ref_4:
.. figure::  /recipes/figures/ref/fire_weather_control_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center
   
   Fire weather control for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.

.. _fig_ref_fire_5:
.. figure::  /recipes/figures/ref/fuel_load_continuity_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center
   
   Fuel load continuity control for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.