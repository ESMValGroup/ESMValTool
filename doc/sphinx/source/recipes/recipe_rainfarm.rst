.. _recipes_rainfarm:

RainFARM stochastic downscaling
===============================


Overview
--------

Precipitation extremes and small-scale variability are essential drivers in many climate change impact studies. However, the spatial resolution currently achieved by global and regional climate models is still insufficient to correctly identify the fine structure of precipitation intensity fields. In the absence of a proper physically based representation, this scale gap can be at least temporarily bridged by adopting a stochastic rainfall downscaling technique (Rebora et al, 2006). With this aim, the Rainfall Filtered Autoregressive Model (RainFARM) was developed to apply the stochastic precipitation downscaling method to climate models. The RainFARM Julia library and command-line tool version (https://github.com/jhardenberg/RainFARM.jl) was implemented as recipe. The stochastic method allows to predict climate variables at local scale from information simulated by climate models at regional scale: It first evaluates the statistical distribution of precipitation fields at regional scale and then applies the relationship to the boundary conditions of the climate model to produce synthetic fields at the requested higher resolution. RainFARM exploits the nonlinear transformation of a Gaussian random precipitation field, conserving the information present in the fields at larger scale (Rebora et al., 2006; D’Onofrio et al., 2014).


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_rainfarm.yml

Diagnostics are stored in diag_scripts/rainfarm/

* rainfarm.jl


User settings
-------------

*Required settings for script*

* slope: spatial spectral slope (set to 0 to compute automatically from large scales)
* nens: number of ensemble members to be calculated
* nf: number of subdivisions for downscaling (e.g. 8 will produce output fields with linear resolution increased by a factor 8)
* conserv_glob: logical, if to conserve precipitation over full domain
* conserv_smooth: logical, if to conserve precipitation using convolution (if neither conserv_glob or conserv_smooth is chosen, box conservation is used)
* weights_climo: set to false or omit if no orographic weights are to be used, else set it to the path to a fine-scale precipitation climatology file. If a relative file path is used, `auxiliary_data_dir` will be searched for this file. The file is expected to be in NetCDF format and should contain at least one precipitation field. If several fields at different times are provided, a climatology is derived by time averaging. Suitable climatology files could be for example a fine-scale precipitation climatology from a high-resolution regional climate model (see e.g. Terzago et al. 2018), a local high-resolution gridded climatology from observations, or a reconstruction such as those which can be downloaded from the WORLDCLIM (http://www.worldclim.org) or CHELSA (http://chelsa-climate.org) websites. The latter data will need to be converted to NetCDF format before being used (see for example the GDAL tools (https://www.gdal.org).


Variables
---------

* pr (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Terzago et al. 2018, Nat. Hazards Earth Syst. Sci., 18, 2825-2840
* D'Onofrio et al. 2014, J of Hydrometeorology 15, 830-843
* Rebora et. al 2006, JHM 7, 724

Example plots
-------------

.. figure:: /recipes/figures/rainfarm/rainfarm.png
   :width: 14cm

   Example of daily cumulated precipitation from the CMIP5 EC-EARTH model on a specific day, downscaled using RainFARM from its original resolution (1.125°) (left panel), increasing spatial resolution by a factor of 8 to 0.14°; Two stochastic realizations are shown (central and right panel). A fixed spectral slope of s=1.7 was used. Notice how the downscaled fields introduce fine scale precipitation structures, while still maintaining on average the original coarse-resolution precipitation. Different stochastic realizations are shown to demonstrate how an ensemble of realizations can be used to reproduce unresolved subgrid variability. (N.B.: this plot was not produced by ESMValTool - the recipe output is netcdf only). 
