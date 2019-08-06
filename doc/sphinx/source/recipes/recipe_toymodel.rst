.. _recipes_toymodel:

Toymodel
====================================================

Overview
--------

The goal of this diagnostic is to simulate single-model ensembles from an observational dataset to investigate the effect of observational uncertain.  For further discussion of this synthetic value generator, its general application to forecasts and its limitations, see Weigel et al. (2008). The output is a netcdf file containing the synthetic observations. Due to the sampling of the perturbations from a Gaussian distribution, running the recipe multiple times, with the same observation dataset and input parameters, will result in different outputs.


Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_toymodel.yml


Diagnostics are stored in diag_scripts/magic_bsc/

* toymodel.R: generates a single model ensemble of synthetic observations




User settings
-------------

User setting files are stored in recipes/

#.	recipe_toymodel.yml

   *Required settings for preprocessor*

	extract_region:

   * start_longitude: minimum longitude
   * end_longitude: maximum longitude
   * start_latitude: minimum longitude
   * end_latitude: maximum latitude
  
  	extract_levels: (for 3D variables)

   * levels: [50000] # e.g. for 500 hPa level
   

   *Required settings for script*

   * number_of_members: integer specifying the number of members to be generated
   * beta: the user defined underdispersion (beta >= 0)


Variables
---------

* any variable (atmos/ocean, daily-monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Bellprat, O., Massonnet, F., Siegert, S., Prodhomme, C., Macias-Gómez, D., Guemas, V., & Doblas-Reyes, F. (2017). Uncertainty propagation in observational references to climate model scales. Remote Sensing of Environment, 203, 101-108.

* Massonet, F., Bellprat, O. Guemas, V., & Doblas-Reyes, F. J. (2016). Using climate models to estimate the quality of global observational data sets. Science, aaf6369.

* Weigel, A. P., Liniger, M. A., & Appenzeller, C. (2008). Can multi-model combinations really enhance the prediction skill of probabilistic ensemble forecasts? Quarterly Journal of the Royal Meteorological Society, 134(630), 241-260.


Example plots
-------------

.. _fig_toymodel:
.. figure::  /recipes/figures/toymodel/synthetic_CMIP5_bcc-csm1-1_Amon_rcp45_r1i1p1_psl_2051-2060.jpg




