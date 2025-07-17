.. _recipes_REF:

CMIP Rapid Evaluation Framework (REF)
======================================

Overview
--------

Here ESMValTool recipes are collected which will be used in the CMIP
`Rapid Evaluation Framework (REF) <https://wcrp-cmip.org/cmip7/rapid-evaluation-framework/>`__.


Available recipes
-----------------

Recipes are stored in `recipes/`

* :ref:`recipe_ecs.yml <recipes_ecs>`:
  Calculate equilibrium climate sensitivity (ECS)
* :ref:`recipe_tcr.yml <recipes_tcr>`:
  Calculate transient climate response (TCR)
* :ref:`recipe_tcre.yml <recipes_tcre>`:
  Calculate transient climate response to cumulative CO2 emissions (TCRE)
* ref/recipe_ref_cre.yml:
  Maps and zonal means of longwave and shortwave cloud radiative effect
* ref/recipe_ref_timeseries_region.yml:
  Time series in comparison with reference data for a selected IPCC region defined through a shape file, based on :ref:`recipe_ref_timeseries.yml <recipe_benchmarking>`
* ref/recipe_ref_annual_cycle_region.yml:
  Annual cycle in comparison with reference data for a selected IPCC region defined through a shape file, based on :ref:`recipe_ref_annual_cycle.yml <recipe_benchmarking>`
* ref/recipe_ref_trend_regions.yml:
  Linear Trends for all IPCC land regions compared with reference data, based on :ref:`recipe_ref_trend_regions.yml <recipes_seaborn_diag>`
* ref/recipe_ref_scatterplot.yml:
  2D histograms with focus on clouds
* ref/recipe_ref_sea_ice_area_basic.yml:
  Seasonal cycle of Arctic (NH) and Antarctic (SH) sea ice area, time series
  of Arctic September (NH) and Antarctic February (SH) sea ice area
* ref/recipe_ref_ozone.yml:
  NH/SH polar cap (60 degrees to 90 degrees) March/September total column ozone
  time series, zonal mean total column ozone vs. time map plot, zonal mean
  total column ozone vs. annual cycle map plot, altitude vs. zonal mean
  ozone profile climatology map plot


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

.. figure::  /recipes/figures/ref/seaborn_jointplot.png
   :align:   center

   2D histogram of total cloud fraction (ESACCI-CLOUD) and shortwave cloud radiative
   effect (CERES-EBAF) for the years 2001-2016 with 1D histograms attached.

.. _fig_ref_4:
.. figure::  /recipes/figures/ref/annual_cycle_sea_ice_area_nh_ambiguous_dataset_ambiguous_mip_historical_r1i1p1f1.png
   :align:   center
   :width:   8cm

   Average seasonal cycle of the Arctic (NH) sea ice area from MPI-ESM1-2-LR
   (red line) compared with OSISAF/CCI (blue line). Created with recipe_ref_sea_ice_area_basic.yml.

.. _fig_ref_5:
.. figure::  /recipes/figures/ref/timeseries_sea_ice_area_nh_sep_ambiguous_dataset_ambiguous_mip_historical_r1i1p1f1.png
   :align:   center
   :width:   8cm

   Time series of Arctic (NH) September (NH) sea ice area from MPI-ESM1-2-LR
   (red line) compared with OSISAF/CCI (blue line). Created with recipe_ref_sea_ice_area_basic.yml.

.. _fig_ref_6:
.. figure::  /recipes/figures/ref/zonal_mean_profile_o3_CNRM-ESM2-1_historical.png
   :align:   center
   :width:   8cm

   Zonal mean vertically resolved ozone climatology from CNRM-ESM2-1 compared with ESACCI-OZONE
   for the years 1990 to 2000. Created with recipe_ref_ozone.yml.

.. _fig_ref_7:
.. figure::  /recipes/figures/ref/timeseries_tas_ambiguous_dataset_Amon_historical_r1i1p1f1.png
   :align:   center
   :width:   8cm

   Time series of near-surface air temperature anomalies from MIROC6 compared with HadCRUT5
   for N.Europe for the years 1980 to 2014 (reference period 1980 to 2009). Created with recipe_ref_timeseries_region.yml.
