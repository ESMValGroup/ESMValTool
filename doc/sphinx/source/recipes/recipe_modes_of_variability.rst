.. _recipes_modes_of_variability:

Modes of variability
====================================================

Overview
--------

The goal of this recipe is to compute modes of variability from a reference/observational dataset and a set of climate projections and calculate the root-mean-square error between the mean anomalies obtained for the clusters from the reference and projection data sets. This is done through K-means clustering applied either directly to the spatial data or after computing the EOFs. The user can specify the number of clusters to be computed. The recipe output consist of netcdf files of the time series of the cluster occurrences, the mean anomaly corresponding to each cluster at each location and the corresponding p-value, for both the observed and projected weather regimes and the RMSE between them. 
 

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_modes_of_variability_wp4.yml


Diagnostics are stored in diag_scripts/magic_bsc/

* WeatherRegime.R - function for computing the EOFs and k-means clusters.

* weather_regime.r: applies the above weather regimes function to the datasets 



User settings
-------------

User setting files are stored in recipes/

#. recipe_modes_of_variability_wp4.yml

   *Required settings for script*

   * start_historical: start date of the reference dataset to be used (please make sure this matches the available data)
   * end_historical: end date of the reference dataset to be used (please make sure this matches the available data)
   * start_projection: start date of the projection dataset to be used
   * end_projection: end date of the projection dataset to be used
   * ncenters: number of centers to be computed by the k-means clustering algorithm
   * detrend_order: the order of the polynomial detrending to be applied
   * EOFs: logical indicating wether the k-means clustering algorithm is applied directly to the spatial data ('FALSE') or to the EOFs ('TRUE')
   * frequency: select the month or season for the diagnostic to be computed for


Variables
---------

* psl or sic (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Dawson, A., T. N. Palmer, and S. Corti, 2012: Simulating regime structures in weather and climate prediction models. Geophysical Research Letters, 39 (21), doi: 10.1029/2012GL053284.

* Ferranti, L., S. Corti, and M. Janousek, 2015: Flow-dependent verification of the ECMWF ensemble over the Euro-Atlantic sector. Quarterly Journal of the Royal Meteorological Society, 141 (688), 916-924, doi: 10.1002/qj.2411.

* Grams, C. M., Beerli, R., Pfenninger, S., Staffell, I., & Wernli, H. (2017). Balancing Europe's wind-power output through spatial deployment informed by weather regimes. Nature climate change, 7(8), 557.

* Hannachi, A., D. M. Straus, C. L. E. Franzke, S. Corti, and T. Woollings, 2017: Low Frequency Nonlinearity and Regime Behavior in the Northern Hemisphere Extra-Tropical Atmosphere. Reviews of Geophysics, doi: 10.1002/2015RG000509.

* Michelangeli, P.-A., R. Vautard, and B. Legras, 1995: Weather regimes: Recurrence and quasi stationarity. Journal of the atmospheric sciences, 52 (8), 1237-1256, doi: 10.1175/1520-0469(1995)052<1237:WRRAQS>2.0.CO

* Vautard, R., 1990: Multiple weather regimes over the North Atlantic: Analysis of precursors and successors. Monthly weather review, 118 (10), 2056-2081, doi: 10.1175/1520-0493(1990)118<2056:MWROTN>2.0.CO;2.

* Yiou, P., K. Goubanova, Z. X. Li, and M. Nogaj, 2008: Weather regime dependence of extreme value statistics for summer temperature and precipitation. Nonlinear Processes in Geophysics, 15 (3), 365-378, doi: 10.5194/npg-15-365-2008.




Example plots
-------------

.. _fig_modesofvar:
.. figure::  /recipes/figures/modes_of_variability/DJF-psl_observed_regimes.png
   :align:   center
   :width:   14cm




