.. _recipe_eyring13jgr:

Ozone and associated climate impacts
====================================

Overview
--------

This recipe is implemented into the ESMValTool to evaluate atmospheric chemistry and the climate impact of stratospheric ozone changes. It reproduces selected plots from Eyring et al. (2013).

The following plots are reproduced:

* Zonal mean of long-term zonal wind with linear trend

.. _`Eyring et al. (2013)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrd.50316

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_eyring13jgr_12.yml

Diagnostics are stored in esmvaltool/diag_scripts/eyring13jgr/

* eyring13jgr_fig12.ncl

User settings in recipe
-----------------------
#. Preprocessor

   * ``zonal`` : Regridding and zonal mean used by eyring13jgr_fig12

#. Script <eyring13jgr_fig12.ncl>

   *Required settings for script*

   * ``e13fig12_exp_MMM``: name of the experiments for the MMM

   *Optional settings for script*

   * ``e13fig12_start_year``: year when to start the climatology calculation
   * ``e13fig12_end_year``: year when to end the climatology calculation
   * ``e13fig12_multimean``: calculate multimodel mean (default: False)
   * ``e13fig12_season``: season (default: ANN (annual))

   *Required settings for variables*
   
   * ``preprocessor``: zonal
   * ``reference_dataset``: name of the reference model or observation for regridding and bias calculation (e.g. ERA5).
   *  ``mip``:  Amon.

Variables
---------

*  ua (atmos, monthly mean, longitude latitude level time)

Observations and reformat scripts
---------------------------------

* ERA5
  *Reformatting with:* recipes/cmorizers/recipe_era5.yml


Example plots
-------------

.. _fig_eyring13jgr_12:
.. figure::  /recipes/figures/eyring13jgr/fig_eyr13jgr_12.png
   :align:   center

   Long-term mean (thin black contour) and linear trend (colour) of zonal mean DJF zonal winds for the multi-model mean CMIP6 over 1995-2014
