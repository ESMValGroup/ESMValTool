.. _recipes_steric_patterns:

Calculate thermosteric sea-level patterns from CMIP6 models
===========================================================

Overview
--------

The recipe recipe_steric_patterns calculates thermosteric sea-level change
patterns from CMIP6 model datasets.

.. note::
  The regrid setting in the recipe is set to a (180, 360) grid to put 
  all models on the same grid.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_steric_patterns.yml

Diagnostics are stored in esmvaltool/diag_scripts/steric_patterns/

* steric_patterns.py: generates thermosteric patterns from input datasets
* sub_funcs.py: set of sub functions to assist with the main script


User settings in recipe
-----------------------

#. Script steric_patterns.py

   *Required settings for script*

   None

   *Optional settings for script*

   None

   *Required settings for variables*

   * short_name
   * additional_datasets
   * exp

   *Optional settings for variables*

   None

   *Required settings for preprocessor*

   * regular_grid: regrids data to a (180, 360) grid

   *Optional settings for preprocessor*

   None


Variables
---------

#. Script steric_patterns.py

* zostoga (ocean, monthly, time)
* zos (ocean, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

None


References
----------

* Palmer, M. D., Gregory, J. M., Bagge, M., Calvert, D., Hagedoorn, J. M., 
  & Howard, T., et al. (2020). Exploring the drivers of global and local 
  sea-level change over the 21st century and beyond. Earth's Future, 8, 
  e2019EF001413. https://doi.org/10.1029/2019EF001413

* Perks, R., & Weeks, J. (2023). MetOffice/ProFSea-tool: v1.0.0 (v1.0.0). 
  Zenodo. https://doi.org/10.5281/zenodo.10255468


Example plots
-------------

.. _fig_steric_patterns_1:
.. figure::  /recipes/figures/steric_patterns/detrended.png
   :align:   center
   :width: 80%

   Detrended zostoga, correcting for model drift using the pre-industrial
   (PiControl) experiment. This is done for each model.

.. _fig_steric_patterns_2:
.. figure::  /recipes/figures/steric_patterns/regressions.png
   :align:   center
   :width: 80%

   Example of the regressions between the global thermal expansion and 
   local dynamic sea-level height for three random grid-cells.

.. _fig_steric_patterns_3:
.. figure::  /recipes/figures/steric_patterns/predictions.png
   :align:   center
   :width: 80%

   Example predictions from the UKESM1-0-LL model patterns for each SSP, 
   as well as a timeseries of globally-averaged mean-squared error.
