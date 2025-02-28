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
* :ref:`recipe_ref_cre.yml <recipe_ref_cre>`:
  Maps and zonal means of longwave and shortwave cloud radiative effect
* :ref:`recipe_model_benchmarking_timeseries_region.yml <recipes_trecipe_model_benchmarking_timeseries_region>`:
  Time series in comparison with refernce data for a region defined through a shape file, based on model_evaluation/recipe_model_benchmarking_timeseries.yml
* :ref:`recipe_model_benchmarking_boxplots_region_trend.yml <recipes_model_benchmarking_boxplots_region_trend>`:
  Benchmarking plot with the linear trend for a region defined through a shape file
* :ref:`recipe_monitor_anncyc_regions.yml <recipe_monitor_anncyc_regions>`:
  Annual cycle for several regions as multi panel plot, based on monitor/recipe_monitor.yml
* :ref:`recipe_portrait_regions.yml <recipe_portrait_regions>`:
    Portrait plot for several regions and variables, based on recipe_portrait.yml


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
.. figure::  /recipes/figures/ref/benchmarking_boxplot_tas_MIROC6_Amon_historical_r1i1p1f1.png
   :align:   center

   Comparing the linear trend for one models (MIROC6 as examlple) to observational data set
   (tas_land: HadCRUT5 and ERA5; sst: HadISST and ERA5; pr: GPCP-SG and ERA5) for the period 2001 to 2014.
   Each box indicates the range from the first quartile to the third quartile, the vertical lines show the median, 
   and the whiskers present the minimum and maximum values.



