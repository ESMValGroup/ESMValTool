.. _ipccwg1ar6ch3_santer20jclim:

Trends in Water Vapour Path
===========================

Overview
--------

This recipe plots the probability density function (PDF) for the linear trend of the Water Vapour Path. It displays the result separately for CMIP5 and CMI6 model results and ensemble members of single models compared to observation and reanalysis data. A filter to use the sampling of observational data on the model results can be applied. It is based on the methods used in Santer et al. 2021.

The following plots are reproduced:

* PDF for the linear trends of the water vapour path for CMIP5 and CMIP6 compared to observations and reanalysis data (tested for RSS and ERA5)

* PDF for the linear trends of the water vapour path for given ensemble members of chosen models compared to observations and reanalysis data (tested for RSS and ERA5)

.. _`Santer et al. 2021`: submitted, but not published, yet

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_ipccwg1ar6ch3_santer20jclim.yml

Diagnostics are stored in esmvaltool/diag_scripts/santer20jclim/

* santer20jclim/santer20jclimfig.py

User settings in recipe
-----------------------
#. Preprocessor

  * ``tropical_ocean`` : Mask out land, regridding on RSS grid and extract region between 20S and 20N


#. Script <santer20jclimfig.py>

   *Required settings for script*


   *Optional settings for script*

   * ``sample_obs``: Path to the mask to filter the model data for a given sampling of observational data. All data sets must be interpolated to the same grid as the mask. The mask must cover at least the time period used for the data. If the keyword is not provided no mask is applied and interpolation to a common grid is not necessary.
   * ``add_model_dist``: List datasets with several ensemble members to plot their PDF. If the keyword is provided an extra plot is created contain their PDFs else only the PDF for CMIP5 and CMIP6 is plotted.

   *Required settings for variables*

   * ``preprocessor``: tropical_ocean, with the optional setting "sample_obs" the regridding to a grid fitting the mask is require, the other preprocessor settings are chosen to match the methods of Santer et al. 2021 but the recipe would accept other settings.
   * ``reference_dataset``: name of the reference data set for regridding, this must be RSS if a mask based on RSS should be applied (given at "sample_obs").

Variables
---------

*  prw (atmos, monthly mean, longitude latitude time)

Observations and reformat scripts
---------------------------------

* ERA5
  *Reformatting with:* on the flight cmorizer as fix (native6)

* RSS
  *Reformatting with:* esmvaltool/cmorizers/obs/cmorize_obs_rss.ncl

* Mask based on RSS
  *Reformatting with:* esmvaltool/cmorizers/obs/cmorize_obs_rssanom.ncl


Example plots
-------------

.. _fig1:
.. figure::  /recipes/figures/santer20jclim/fig1.png
   :align:   center

   PDF of the linear trends (1988-2019) of the water vapour path for CMIP5 and CMIP6 compared to RSS and ERA5. For model data historical and future scenario (RCP8.5 for CMIP5 and SSP585 for CMIP6 are combined). A masked derived from the RSS anomalies and climatology is applied to all data sets. The fit is calculated as kernel-density estimate with Gaussian kernels using Scott’s Rule. If several ensemble members are included for one model they are weighted with the inverse of the number of ensemble members.
