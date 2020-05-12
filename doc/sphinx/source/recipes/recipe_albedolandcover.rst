.. _recipes_albedolandcover:

Landcover - Albedo
==================


Overview
--------

The diagnostic determines the coefficients of multiple linear regressions fitted between the albedo values and the tree, shrub, short vegetation (crops and grasses) fractions of each grid cell within spatially moving windows encompassing 5x5 model grid cells. Solving these regressions provides the albedo values for trees, shrubs and short vegetation (crops and grasses) from which the albedo changes associated with transitions between these three landcover types are derived. The diagnostic distinguishes between snow-free and snow-covered grid cells.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_albedolandcover.yml

Diagnostics are stored in diag_scripts/landcover/

    * albedolandcover.py


User settings
-------------

Several parameters can be set in the recipe


Variables
---------

* rsus (atmos, monthly mean, time latitude longitude)
* rsds (atmos, monthly mean, time latitude longitude)
* snc (landice, monthly mean, time latitude longitude)
* grassFrac (land, monthly mean, time latitude longitude)
* treeFrac (land, monthly mean, time latitude longitude)
* shrubFrac (land, monthly mean, time latitude longitude)
* cropFrac (land, monthly mean, time latitude longitude)
* pastureFrac (land, monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

A reformatting script for observational data is available here:
    * cmorize_obs_duveiller2018.py


References
----------

* Duveiller, G., Hooker, J. and Cescatti, A., 2018a. A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5: 180014.

* Duveiller, G., Hooker, J. and Cescatti, A., 2018b. The mark of vegetation change on Earth’s surface energy balance. Nature communications, 9(1): 679.

Example plots
-------------

.. _fig_landcoveralbedo_CMIP5_MPI-ESM-LR:
.. figure::  /recipes/figures/albedolandcover/MPI-ESM-LR_albedo_change_from_tree_to_crop-grass.png
   :align:   center
   :width:   14cm

   Example of albedo change from tree to crop and grass for the CMIP5 model MPI-ESM-LR derived for the month of July and averaged over the years 2000 to 2004.
