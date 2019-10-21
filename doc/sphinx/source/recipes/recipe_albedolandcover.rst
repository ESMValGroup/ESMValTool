.. _recipes_albedolandcover:

Landcover - Albedo
==================


Overview
--------

The diagnostic computes the relationship between landcover and albedo using a multiple linear regression technique.


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
* grassFrac    (land, monthly mean, time latitude longitude)
* treeFrac     (land, monthly mean, time latitude longitude)
* shrubFrac    (land, monthly mean, time latitude longitude)
* cropFrac     (land, monthly mean, time latitude longitude)
* pastureFrac  (land, monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

A reformatting script for observational data is available here:
    * cmorize_obs_duveiller2018.py


References
----------

* Duveiller, G., Hooker, J. and Cescatti, A., 2018a. A dataset mapping the potential biophysical effects of vegetation cover change. Scientific Data, 5: 180014.

* Duveiller, G., Hooker, J. and Cescatti, A., 2018b. The mark of vegetation change on Earth’s surface energy balance. Nature communications, 9(1): 679.
