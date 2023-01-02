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

MARRMoT requires potential evapotranspiration (evspsblpot). The variable evspsblpot is not available in ERA-Interim. Thus, we use the debruin function (De Bruin et al. 2016) to obtain evspsblpot using both ERA-Interim and ERA5. This function needs the variables tas, psl, rsds, and rsdt as input.

wflow_sbm and wflow_topoflex
****************************
Forcing data for the `wflow_sbm <https://wflow.readthedocs.io/en/latest/wflow_sbm.html>`_
and `wflow_topoflex <https://wflow.readthedocs.io/en/latest/wflow_topoflex.html>`_
hydrological models can be prepared using recipe_wflow.yml.
If PET is not available from the source data (e.g. ERA-Interim), then it can be derived from psl, rsds and rsdt using De Bruin's 2016 formula (De Bruin et al. 2016). For daily ERA5 data, the time points of these variables are shifted 30 minutes with respect to one another. This is because in ERA5, accumulated variables are recorded over the past hour, and in the process of cmorization, we shift the time coordinates to the middle of the interval over which is accumulated. However, computing daily statistics then averages the times, which results in 12:00 UTC for accumulated variables and 11:30 UTC for instantaneous variables. Therefore, in this diagnostic, the time coordinates of the daily instantaneous variables are shifted 30 minutes forward in time.

LISFLOOD
********
`LISFLOOD <https://ec-jrc.github.io/lisflood-model/>`_ is a spatially distributed water resources model, developed by the Joint Research Centre (JRC) of the European Commission since 1997. We provide a recipe to produce meteorological forcing data for the Python 3 version of LISFLOOD.

LISFLOOD has a separate preprocessor LISVAP that derives some additional variables. We don't replace LISVAP. Rather, we provide input files that can readily be passed to LISVAP and then to LISFLOOD.


HYPE
****

The hydrological catchment model HYPE simulates water flow and substances on their way from precipitation through soil, river and lakes to the river outlet.
HYPE is developed at the Swedish Meteorological and Hydrological Institute. The recipe pre-processes ERA-Interim and ERA5 data for use in HYPE.

GlobWat
*******
GlobWat is a soil water balance model that has been provided by the Food and Agriculture Organization (FAO) to assess water use in irrigated agriculture (http://www.fao.org/nr/water/aquamaps). The recipe pre-processes ERA-Interim and ERA5 reanalyses data for use in the GlobWat model. GlobWat requires potential evapotranspiration (evspsblpot) as input. The variable evspsblpot is not available in ERA-Interim. Thus, we use debruin function (De Bruin et al. 2016) or the langbein method (Langbein et al. 1949) to obtain evspsblpot using both ERA-Interim and ERA5. The Langbein function needs a variable tas and the debruin function besides that needs the variables psl, rsds, and rsdt as input. In order to calculate monthly/daily pet with Langbein method we assumed that tas is constant over time and the average value is equal to the annual average.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/hydrology

    * recipe_pcrglobwb.yml
    * recipe_marrmot.yml
    * recipe_wflow.yml
    * recipe_lisflood.yml
    * recipe_hype.yml
    * recipe_globwat.yml

Diagnostics are stored in esmvaltool/diag_scripts/hydrology

    * pcrglobwb.py
    * marrmot.py
    * wflow.py
    * lisflood.py
    * hype.py
    * globwat.py 


User settings in recipe
-----------------------

All hydrological recipes require a shapefile as an input to produce forcing data. This shapefile determines the shape of the basin for which the data will be cut out and processed. All recipes are tested with `the shapefiles <https://github.com/eWaterCycle/recipes_auxiliary_datasets/tree/main/>`_  that are used for the eWaterCycle project. In principle any shapefile can be used, for example, the freely available basin shapefiles from the `HydroSHEDS project <https://www.hydrosheds.org/>`_. 

#. recipe_pcrglobwb.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979

#. recipe_marrmot.yml

   There is one diagnostic ``diagnostic_daily`` for using daily data.

   *Required preprocessor settings:*

      The settings below should not be changed.

      *extract_shape:*

         * shapefile: Meuse.shp (MARRMoT is a hydrological Lumped model that needs catchment-aggregated forcing data. The catchment is provided as a shapefile, the path can be relative to ``auxiliary_data_dir`` as defined in config-user.yml.).
         * method: contains
         * crop: true

   *Required diagnostic script settings:*

      * basin: Name of the catchment

#. recipe_wflow.yml

   *Optional preprocessor settings:*

      * extract_region: the region specified here should match the catchment

   *Required diagnostic script settings:*

	    * basin: name of the catchment
	    * dem_file: netcdf file containing a digital elevation model with
	      elevation in meters and coordinates latitude and longitude.
              A wflow example dataset is available at: https://github.com/openstreams/wflow/tree/master/examples/wflow_rhine_sbm
              The example dem_file can be obtained from https://github.com/openstreams/wflow/blob/master/examples/wflow_rhine_sbm/staticmaps/wflow_dem.map 
	    * regrid: the regridding scheme for regridding to the digital elevation model. Choose ``area_weighted`` (slow) or ``linear``.

#. recipe_lisflood.yml

   *Required preprocessor settings:*

      * extract_region: A region bounding box slightly larger than the shapefile. This is run prior to regridding, to save memory.
      * extract_shape:*

         * shapefile: A shapefile that specifies the extents of the catchment.

         These settings should not be changed

         * method: contains
         * crop: true

      * regrid:*

         * target_grid: Grid of LISFLOOD input files

         These settings should not be changed

         * lon_offset: true
         * lat_offset: true
         * scheme: linear

   There is one diagnostic ``diagnostic_daily`` for using daily data.

   *Required diagnostic script settings:*

      * catchment: Name of the catchment, used in output filenames

#. recipe_hype.yml

   *Required preprocessor settings:*

   * start_year: 1979
   * end_year: 1979
   * shapefile: Meuse_HYPE.shp (expects shapefile with subcatchments)

   These settings should not be changed

   * method: contains
   * decomposed: true

#. recipe_globwat.yml

   *Required preprocessor settings:*

   * start_year: 2004
   * end_year: 2004
   * target_grid_file: grid of globwat input files. A target file has been generated from one of the GlobWat models sample files (prc01wb.asc) for regridding ERA5 and ERA-Interim datasets. The ASCII file can be found at: https://storage.googleapis.com/fao-maps-catalog-data/geonetwork/aquamaps/GlobWat-InputP1_prec.zip. You can use the GDAL translator to convert the file from ASCII format to NetCDF format by entering the following command into the terminal: gdal_translate -of netCDF prc01wb.asc globwat_target_grid.nc

   *Optional preprocessor settings:*

   * area_selection: A region bounding box to extract the data for a specific region. The area selection preprocessor can be used by users to process the data for their desired region. The data will be processed at the global scale if the preprocessor in the recipe is commented.
   * regrid_scheme: The area-weighted regridding scheme is used as a default regridding scheme to ensure that the total volume of water is consistent before and after regridding.
   * langbein_pet: Can be set to True to use langbein function for calculating evspsblpot (default is de bruin method)


Variables
---------

#. recipe_pcrglobwb.yml

   * tas (atmos, daily, longitude, latitude, time)
   * pr (atmos, daily, longitude, latitude, time)

#. recipe_marrmot.yml

   * pr (atmos, daily or hourly mean, longitude, latitude, time)
   * psl (atmos, daily or hourly mean, longitude, latitude, time)
   * rsds (atmos, daily or hourly mean, longitude, latitude, time)
   * rsdt (atmos, daily or hourly mean, longitude, latitude, time)
   * tas (atmos, daily or hourly mean, longitude, latitude, time)

#. recipe_wflow.yml

   * orog (fx, longitude, latitude)
   * pr (atmos, daily or hourly mean, longitude, latitude, time)
   * tas (atmos, daily or hourly mean, longitude, latitude, time)

   Either potential evapotranspiration can be provided:

   * evspsblpot(atmos, daily or hourly mean, longitude, latitude, time)

   or it can be derived from tas, psl, rsds, and rsdt using the De Bruin formula, in that case the following variables need to be provided:

   * psl (atmos, daily or hourly mean, longitude, latitude, time)
   * rsds (atmos, daily or hourly mean, longitude, latitude, time)
   * rsdt (atmos, daily or hourly mean, longitude, latitude, time)

#. recipe_lisflood.yml

   * pr (atmos, daily, longitude, latitude, time)
   * tas (atmos, daily, longitude, latitude, time)
   * tasmax (atmos, daily, longitude, latitude, time)
   * tasmin (atmos, daily, longitude, latitude, time)
   * tdps (atmos, daily, longitude, latitude, time)
   * uas (atmos, daily, longitude, latitude, time)
   * vas (atmos, daily, longitude, latitude, time)
   * rsds (atmos, daily, longitude, latitude, time)

#. recipe_hype.yml

   * tas (atmos, daily or hourly, longitude, latitude, time)
   * tasmin (atmos, daily or hourly, longitude, latitude, time)
   * tasmax (atmos, daily or hourly, longitude, latitude, time)
   * pr (atmos, daily or hourly, longitude, latitude, time)

#. recipe_globwat.yml

   * pr (atmos, daily or monthly, longitude, latitude, time)
   * tas (atmos, daily or monthly, longitude, latitude, time)
   * psl (atmos, daily or monthly, longitude, latitude, time)
   * rsds (atmos, daily or monthly, longitude, latitude, time)
   * rsdt (atmos, daily or monthly , longitude, latitude, time)

Observations and reformat scripts
---------------------------------
*Note: download instructions can be obtained with `esmvaltool data info DATASET` or in headers of cmorization scripts.*

*  ERA-Interim (esmvaltool/cmorizers/data/formatters/datasets/era_interim.py)
*  ERA5 (esmvaltool/diag_scripts/cmorizers/era5.py)

Output
---------

#. recipe_pcrglobwb.yml

#. recipe_marrmot.yml

    The forcing data, the start and end times of the forcing data, the latitude and longitude of the catchment are saved in a .mat file as a data structure readable by MATLAB or Octave.

#. recipe_wflow.yml

	The forcing data, stored in a single NetCDF file.

#. recipe_lisflood.yml

   The forcing data, stored in separate files per variable.

#. recipe_globwat.yml

   The forcing data, stored in separate files per timestep and variable.

References
----------

* Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P.: PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018, 2018.
* De Bruin, H. A. R., Trigo, I. F., Bosveld, F. C., Meirink, J. F.: A Thermodynamically Based Model for Actual Evapotranspiration of an Extensive Grass Field Close to FAO Reference, Suitable for Remote Sensing Application, American Meteorological Society, 17, 1373-1382, DOI: 10.1175/JHM-D-15-0006.1, 2016.
* Arheimer, B., Lindström, G., Pers, C., Rosberg, J. och J. Strömqvist, 2008. Development and test of a new Swedish water quality model for small-scale and large-scale applications. XXV Nordic Hydrological Conference, Reykjavik, August 11-13, 2008. NHP Report No. 50, pp. 483-492.
* Lindström, G., Pers, C.P., Rosberg, R., Strömqvist, J., Arheimer, B. 2010. Development and test of the HYPE (Hydrological Predictions for the Environment) model – A water quality model for different spatial scales. Hydrology Research 41.3-4:295-319.
* van der Knijff, J. M., Younis, J. and de Roo, A. P. J.: LISFLOOD: A GIS-based distributed model for river basin scale water balance and flood simulation, Int. J. Geogr. Inf. Sci., 24(2), 189–212, 2010.
* Hoogeveen, J., Faurès, J. M., Peiser, L., Burke, J., de Giesen, N. V.: GlobWat--a global water balance model to assess water use in irrigated agriculture, Hydrology & Earth System Sciences Discussions, 2015 Jan 1;12(1), Doi:10.5194/hess-19-3829-2015.
* Langbein, W.B., 1949. Annual runoff in the United States. US Geol. Surv.(https://pubs.usgs.gov/circ/1949/0052/report.pdf)
