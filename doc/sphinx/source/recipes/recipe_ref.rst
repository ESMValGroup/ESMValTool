.. _recipes_REF:

CMIP Rapid Evaluation Framework (REF)
======================================

Overview
--------

Here ESMValTool recipes are collected which will be used in the CMIP
`Rapid Evaluation Framework (REF) <https://wcrp-cmip.org/cmip7/rapid-evaluation-framework/>`__.


Available recipes
-----------------

Recipes are stored in recipes

* :ref:`recipe_ecs.yml <recipes_ecs>`:
  Calculate equilibrium climate sensitivity (ECS)
* :ref:`recipe_tcr.yml <recipes_tcr>`:
  Calculate transient climate response (TCR)
* :ref:`recipe_tcre.yml <recipes_tcre>`:
  Calculate transient climate response to cumulative CO2 emissions (TCRE)
* ref/recipe_ref_cre.yml:
  Maps and zonal means of longwave and shortwave cloud radiative effect
* ref/recipe_ref_sea_ice_seasonal.yml:
  Seasonal cycle of Arctic (NH) and Antarctic (SH) sea ice area and extent


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
.. figure::  /recipes/figures/ref/annual_cycle_siextent_nh_ambiguous_dataset_SImon_historical_r1i1p1f1.png
   :align:   center
   :width:   8cm

   Average seasonal cycle of the Arctic (NH) sea ice extent from MPI-ESM1-2-LR 
   (red line) compared with ESACCI-SEAICE (black line).
