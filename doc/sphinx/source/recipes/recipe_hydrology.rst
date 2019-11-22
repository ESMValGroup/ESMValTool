.. _recipes_hydrology:

Hydrological models - data pre-processing
=========================================

Overview
--------

We provide a collection of scripts that pre-processes environmental data for use in several hydrological models:

PCR-GLOBWB
**********
PCR-GLOBWB (PCRaster Global Water Balance) is a large-scale hydrological model intended for global to regional studies and developed at the Department of Physical Geography, Utrecht University (Netherlands). The recipe pre-processes ERA-Interim reanalyses data for use in the PCR-GLOBWB.

wflow_sbm and wflow_topoflex
****************************
Forcing data for the `wflow_sbm <https://wflow.readthedocs.io/en/latest/wflow_sbm.html>`_
and `wflow_topoflex <https://wflow.readthedocs.io/en/latest/wflow_topoflex.html>`_
hydrological models can be prepared using recipe_wflow.yml.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_pcrglobwb.yml
    * recipe_wflow.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * pcrglobwb.py
    * wflow.py

User settings in recipe
-----------------------

#. recipe_pcrglobwb.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979

#. recipe_wflow.yml

   *Required preprocessor settings:*

      * extract_region: the region specified here should match the catchment
      * daily_statistics: if the frequency of the input data is not daily, it
        should be converted to daily using the preprocessor function
        daily_statistics with ``operator: mean``.

   *Required diagnostic script settings:*

	    * basin_name: name of the catchment
	    * dem_file: netcdf file containing a digital elevation model with
	      elevation in meters and coordinates latitude and longitude.

Variables
---------

#. recipe_pcrglobwb.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)

#. recipe_wflow.yml

   * orog (fx, longitude, latitude)
   * pr (atmos, daily or hourly mean, longitude, latitude, time)
   * psl (atmos, daily or hourly mean, longitude, latitude, time)
   * rsds (atmos, daily or hourly mean, longitude, latitude, time)
   * rsdt (atmos, daily or hourly mean, longitude, latitude, time)
   * tas (atmos, daily or hourly mean, longitude, latitude, time)

Observations and reformat scripts
---------------------------------
*Note: see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for download instructions.*

*  ERA-Interim (esmvaltool/cmorizers/obs/cmorize_obs_era_interim.py)
*  ERA5 (esmvaltool/cmorizers/obs/cmorize_obs_era5.py)

References
----------

* Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P.: PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018, 2018.
