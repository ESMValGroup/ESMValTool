.. _recipes_climate_patterns:

Generating Climate Patterns from CMIP6 Models
=============================================

Overview
--------

The recipe recipe_climate_patterns generates climate patterns from CMIP6 model
datasets. It also generates a set of parameters to tune the energy-balance
model in IMOGEN-JULES.

.. note::
  The regrid setting in the recipe is set to a 2.5x3.75 grid. This is done to
  match the current resolution in the IMOGEN-JULES framework, but can be
  adjusted with no issues for a finer/coarser patterns grid.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_climate_patterns.yml

Diagnostics are stored in esmvaltool/diag_scripts/climate_patterns/

    * climate_patterns.py: generates climate patterns from input datasets
    * ebm_parameters.py: outputs a set of ebm parameters for IMOGEN-JULES
    * rename_variables.py: renames variables depending on user specifications
    * sub_functions.py: set of sub functions to assist with driving scripts
    * plotting.py: contains all plotting functions for driving scripts


User settings in recipe
-----------------------

#. Script climate_patterns.py

   *Required settings for script*

   None

   *Optional settings for script*

   * grid: whether you want to remove Antarctic latitudes or not
   * imogen_mode: output imogen-specific var names + .nc files
   * output_r2_scores: output measures of pattern robustness (adds runtime)
   * parallelise: parallelise over models or not
   * parallel_threads: if you want to paralellise, how many threads you want

   *Required settings for variables*

   * short_name
   * additional_datasets

   *Optional settings for variables*

   None

   *Required settings for preprocessor*

   * monthly_statistics: converts data to mean monthly data

   *Optional settings for preprocessor*

   * regrid: regrids data

#. Script ebm_parameters.py

   *Required settings for script*

   None

   *Optional settings for script*

   * include_params: includes input parameters if known
   * parallelise: parallelise over models or not
   * parallel_threads: if you want to paralellise, how many threads you want

   *Required settings for variables*

   * short_name
   * mip (only for land_frac)
   * project (only for land_frac)
   * additional_datasets

   *Optional settings for variables*

   None

   *Required settings for preprocessor*

   * annual_statistics: converts data to mean annual data

   *Optional settings for preprocessor*

   None


Variables
---------

#. Script climate_patterns.py

* tasmax (atmos, monthly, longitude latitude time)
* tasmin (atmos, monthly, longitude latitude time)
* tas (atmos, monthly, longitude latitude time)
* hurs (atmos, monthly, longitude latitude time)
* huss (atmos, monthly, longitude latitude time)
* pr (atmos, monthly, longitude latitude time)
* sfcWind (atmos, monthly, longitude latitude time)
* ps (atmos, monthly, longitude latitude time)
* rsds (atmos, monthly, longitude latitude time)
* rlds (atmos, monthly, longitude latitude time)

#. Script ebm_parameters.py

* land_frac (land, longitude latitude)
* tas (atmos, annual, longitude latitude time)
* rlut (atmos, annual, longitude latitude time)
* rsut (atmos, annual, longitude latitude time)
* rsdt (atmos, annual, longitude latitude time)


Observations and reformat scripts
---------------------------------

None

References
----------

* Huntingford, C., Cox, P. An analogue model to derive additional climate
  change scenarios from existing GCM simulations.
  Climate Dynamics 16, 575–586 (2000). https://doi.org/10.1007/s003820000067

Example plots
-------------

.. _fig_climate_patterns_1:
.. figure::  /recipes/figures/climate_patterns/ebm_plots.png
   :align:   center
   :width: 80%

   Linear regression between tas and rtmt for the 4x abrupt CO2 expriment (top),
   derived model total radiative forcing SSP1-26 and SSP5-85 (middle), EBM's pediction
   of global surface temperature vs model output. (bottom)

.. _fig_climate_patterns_2:
.. figure::  /recipes/figures/climate_patterns/patterns.png
   :align:   center
   :width: 80%

   Patterns generated for CMIP6 models, gridded view. Patterns are shown per
   variable, for the month of January.

.. _fig_climate_patterns_3:
.. figure::  /recipes/figures/climate_patterns/score_timeseries.png
   :align:   center
   :width: 80%

   R2 scores of patterns fitting per variable. Diversity of scores sits in the
   literatures' range: with temperature, specific humidity and longwave
   downwelling radiation being the most robust fits.
