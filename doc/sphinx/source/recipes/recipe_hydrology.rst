.. _recipes_hydrology:

Hydrological models - data pre-processing
=========================================

Overview
--------

We provide a collection of scripts that pre-processes environmental data for use in several hydrological models:

PCR-GLOBWB
**********
PCR-GLOBWB (PCRaster Global Water Balance) is a large-scale hydrological model intended for global to regional studies and developed at the Department of Physical Geography, Utrecht University (Netherlands). The recipe pre-processes ERA-Interim reanalyses data for use in the PCR-GLOBWB.

MARRMoT
**********
MARRMoT (Modular Assessment of Rainfall-Runoff Models Toolbox) is a rainfall-runoff model comparison framework that allows objective comparison between different conceptual hydrological model structures https://github.com/wknoben/MARRMoT. The recipe pre-processes ERA-Interim and ERA5 reanalyses data for use in the MARRMoT.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_pcrglobwb.yml
    * recipe_marrmot.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * pcrglobwb.py
    * marrmot.py


User settings in recipe
-----------------------

#. recipe_pcrglobwb.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979



#. recipe_marrmot.yml

   There are two diagnostics, one for daily and one for hourly data.

   *Required preprocessor settings:*

      The settings below should not be changed.

      *extract_shape:*

         * shapefile: meuse_hydrosheds.shp (MARRMoT is a hydrological Lumped model that needs catchment-aggregated forcing data. The catchment is provided as a shapefile, the path can be relative to ``auxiliary_data_dir`` as defined in config-user.yml.).
         * method: contains
         * crop: true

      *daily_statistics:*

         * operator: mean (MARRMoT needs daily forcing data. Hourly forcing data are converted to daily values by mean operator).

   *Required diagnostic script settings:*

      * basin: Name of the catchment


Variables
---------

#. recipe_pcrglobwb.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)


#. recipe_marrmot.yml

   * tas ( longitude, latitude, time)
   * pr (longitude, latitude, time)

  MARRMoT requires potential evapotranspiration (PET). The variable PET is not available in ERA-Interim archive. Thus, we use the debruin function (De Bruin et al. 2016) to obtain PET using both ERA-Interim and ERA5. This function needs the psl, rsds, and rsdt variables.

   * psl (longitude, latitude, time)
   * rsds (longitude, latitude, time)
   * rsdt (longitude, latitude, time)


Observations and reformat scripts
---------------------------------
*Note: see headers of cmorization scripts (in esmvaltool/cmorizers/obs) for download instructions.*

*  ERA-Interim (esmvaltool/cmorizers/obs/cmorize_obs_ERA-Interim.py)
*  ERA5 (esmvaltool/cmorizers/obs/cmorize_obs_ERA5.py)


Output
---------

#. recipe_pcrglobwb.yml


#. recipe_marrmot.yml

    The forcing data, the start and end times of the forcing data, the latitude and longitude of the catchment are saved in a .mat file as a data structure readable by MATLAB or Octave.


References
----------

* Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P.: PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018, 2018.
* De Bruin, H. A. R., Trigo, I. F., Bosveld, F. C., Meirink, J. F.: A Thermodynamically Based Model for Actual Evapotranspiration of an Extensive Grass Field Close to FAO Reference, Suitable for Remote Sensing Application, American Meteorological Society, 17, 1373-1382, DOI: 10.1175/JHM-D-15-0006.1, 2016.
