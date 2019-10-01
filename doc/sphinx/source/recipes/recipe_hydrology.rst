.. _recipes_hydrology:

Data pre-processing for hydrological models
============================================

Overview
--------

We provide a collection of scripts that pre-processes environmental data for use in several hydrological models:

PCR-GLOBWB
**********
PCR-GLOBWB (PCRaster Global Water Balance) is a large-scale hydrological model intended for global to regional studies and developed at the Department of Physical Geography, Utrecht University (Netherlands). The recipe pre-processes ERA-Interim reanalyses data for use in the PCR-GLOBWB.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_pcrglobwb.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * pcrglobwb.py


User settings in recipe
-----------------------

#. recipe_pcrglobwb.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979

Variables
---------

#. recipe_pcrglobwb.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)

Observations and reformat scripts
---------------------------------
*Note: see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for download instructions.*

*  ERA-Interim (tas, pr - esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)

References
----------

* Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P.: PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018, 2018.

HYPE
****

The hydrological catchment model HYPE simulates water flow and substances on their way from precipitation through soil, river and lakes to the river outlet.
HYPE is developed at the Swedish Meteorological and Hydrological Institute. The recipe pre-processes ERA-Interim reanalyses data for use in the HYPE.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_hype.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * hype.py


User settings in recipe
-----------------------

#. recipe_hype.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979

Variables
---------

#. recipe_hype.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)

Observations and reformat scripts
---------------------------------
*Note: see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for download instructions.*

*  ERA-Interim (tas, pr - esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)

References
----------

* Arheimer, B., Lindström, G., Pers, C., Rosberg, J. och J. Strömqvist, 2008. Development and test of a new Swedish water quality model for small-scale and large-scale applications. XXV Nordic Hydrological Conference, Reykjavik, August 11-13, 2008. NHP Report No. 50, pp. 483-492.
* Lindström, G., Pers, C.P., Rosberg, R., Strömqvist, J., Arheimer, B. 2010. Development and test of the HYPE (Hydrological Predictions for the Environment) model – A water quality model for different spatial scales. Hydrology Research 41.3-4:295-319.
