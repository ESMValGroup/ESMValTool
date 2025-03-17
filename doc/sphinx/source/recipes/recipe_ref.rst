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
* ref/recipe_ref_scatterplot.yml:
  2D histograms with focus on clouds


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
.. figure::  /recipes/figures/ref/seaborn_jointplot.png
   :align:   center

   2D histogram of total cloud fraction (ESACCI-CLOUD) and shortwave cloud radiative
   effect (CERES-EBAF) for the years 2001-2016 with 1D histograms attached.
