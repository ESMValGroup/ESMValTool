.. _ipccwg1ar6ch3_santer21jclim:

Trends in Water Vapour Path
===========================

Overview
--------

This recipe reproduces Figure 3.12 from IPCC AR6 WGI Chapter 3: Human Influence on the Climate System (`Eyring et al., 2021`_).
The plot is based on the methods used in `Santer et al., 2007`_ and `Santer et al., 2021`_.
It shows the probability density function (PDF) for the linear trend of the Water Vapour Path.
The results for CMIP5 and CMIP6 models are displayed separately and are compared to observation and reanalysis data.
The probability density function is estimated by using Gaussian kernel density estimation.
The models are weighted by their respective number of realisations.
An additional plot allows to show ensemble members of a few select models.
A filter to use the sampling of observational data on the model results can be applied.

The following plots are produced:

* PDF for the linear trends of the water vapour path for CMIP5 and CMIP6 compared to observations and reanalysis data (tested for RSS and ERA5)

* PDF for the linear trends of the water vapour path for given ensemble members of chosen models compared to observations and reanalysis data (tested for RSS and ERA5)

.. _Eyring et al., 2021: https://www.ipcc.ch/report/ar6/wg1/chapter/chapter-3/
.. _Santer et al., 2007: https://www.pnas.org/doi/full/10.1073/pnas.0702872104
.. _Santer et al., 2021: https://journals.ametsoc.org/view/journals/clim/34/15/JCLI-D-20-0768.1.xml

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/ipccwg1ar6ch3/

* recipe_ipccwg1ar6ch3_santer21jclim.yml

Diagnostics are stored in esmvaltool/diag_scripts/ipcc_ar6/

* santer21jclimfig.py

User settings in recipe
-----------------------
#. Preprocessor

   * ``tropical_ocean`` : Mask out land, regridding on RSS grid and extract region between 50S and 50N


#. Script <santer21jclimfig.py>

   *Required settings for script*


   *Optional settings for script*

   * ``sample_obs``: Path to the mask to filter the model data for a given sampling of observational data.
     All data sets must be interpolated to the same grid as the mask. The mask must cover at least the time period used for the data.
     If the keyword is not provided no mask is applied and interpolation to a common grid is not necessary.
   * ``add_model_dist``: List datasets with several ensemble members to plot their PDF.
     If the keyword is provided an extra plot is created containing their PDFs, else only the PDFs for CMIP5 and CMIP6 are plotted.
   * ``histmin`` and ``histmax``: The minimum and maximum value of the histograms, also acting as x-axis limits.
     If no values are provided, they are automatically set to the minimum and maximum values of the model data.
   * ``ymax``: The upper y-axis boundary (lower is set to zero). If no value is provided, this is also set automatically.

   *Required settings for variables*

   * ``preprocessor``: tropical_ocean, with the optional setting "sample_obs"
     the regridding to a grid fitting the mask is required, the other
     preprocessor settings are chosen to match the methods of `Santer et al., 2021`_
     but the recipe would accept other settings.
   * ``reference_dataset``: name of the reference data set for regridding,
     this must be RSS if a mask based on RSS should be applied
     (given at "sample_obs").

Variables
---------
Tested for:

*  prw (atmos, monthly mean, longitude latitude time)
*  tas (atmos, monthly mean, longitude latitude time)
*  pr (atmos, monthly mean, longitude latitude time)

Other variables should be possible.

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
.. figure::  /recipes/figures/santer21jclim/fig1.png
   :align:   center

   PDF of the linear trends (1988-2019) of the water vapour path for CMIP5 and CMIP6 compared to RSS and ERA5.
   For model data historical and future scenario (RCP8.5 for CMIP5 and SSP585 for CMIP6 are combined).
   A masked derived from the RSS anomalies and climatology is applied to all data sets.
   The fit is calculated as kernel-density estimate with Gaussian kernels using Scott’s Rule.
   If several ensemble members are included for one model they are weighted with the inverse of the number of ensemble members.

.. _fig2:
.. figure::  /recipes/figures/santer21jclim/fig2.png
   :align:   center

   PDF of the linear trends (1988-2019) of the water vapour path for the ensemble members
   of 'MPI-ESM1-2-LR' and 'CanESM5' models compared to ERA5.
