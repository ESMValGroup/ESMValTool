.. _recipes_insurance_risk_index_wp7:

Combined Climate Extreme Index
====================================================

Overview
--------

The goal of this diagnostic is to compute time series of a number of extreme events: heatwave, coldwave, heavy precipitation, drought and high wind. Then, the user can combine these different components (with or without weights). The result is an index similar to the Climate Extremes Index (CEI; Karl et al., 1996), the modified CEI (mCEI; Gleason et al., 2008) or the Actuaries Climate Index (ACI; American Academy of Actuaries, 2018). The output consists of a netcdf file containing the area-weighted and multi-model multi-metric index. This recipe can be applied to data with any temporal resolution, and the running average is computed based on the user-defined window length (e.g. a window length of 5 would compute the 5-day running mean when applied to monthly data, or 5-month running mean when applied to monthly data).

In recipe_extreme_index_wp7.yml, after defining the area and reference and projection period, the metric indicating the extreme index is selected. The options are
* t90p to compute the number of days when the maximum temperature exceeds the 90th percentile,
* t10p to compute the number of days when the minimum temperature falls below the 10th percentile,
* Wx to compute the number of days when wind power (third power of wind speed) exceeds the 90th percentile,
* cdd to compute the maximum length of a dry spell, defined as the maximum number of consecutive days when the daily precipitation is lower than 1 mm, and
* rx5day to compute the maximum precipitation accumulated during 5 consecutive days.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_combined_indices_wp6.yml

* recipe_extreme_index_wp7.yml

Diagnostics are stored in diag_scripts/magic_bsc/

* combined_indices_wp6.r : calculates the area-weighted means and multi-model means, with or without weights

* risk_index.r



User settings
-------------

User setting files are stored in recipes/

#. recipe_combined_indices_wp6.yml

   *Required settings for script*

   * weights: either ‘equal’, for equal weights, ‘null’ for no weights, or a vector of integers the same length as the number of input datasets.
   * running_mean: an integer specifying the length of the window to be used for computing the running mean (does not work yet).
   * moninf: instead of running_mean an integer can be given to determine the first month of the seasonal mean to be computed (does not work yet).
   * monsup: an integer specifying the last month to be computed (does not work yet).
   * Multi_year_average: ‘true’ or ‘false’ to specify whether to compute the mean across all input years (does not work yet).

#. recipe_extreme_index_wp7.yml

   *Required settings for script*

   * metric: the metric to be computed, t90p, t10p, Wx, cdd, rx5day. See overview for a description of the different metrics (cdd does not work yet).


Variables
---------

* tasmax, tasmin, pr or sfcWind (atmos, daily, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*

References
----------

* Alexander L.V.  and Coauthors (2006). Global observed changes in daily climate extremes of temperature and precipitation. J. Geophys. Res., 111, D05109, doi:10.1029/2005JD006290.

* American Academy of Actuaries, Canadian Institute of Actuaries, Casualty Actuarial Society and Society of Actuaries. Actuaries Climate Index. http://actuariesclimateindex.org (2018-10-06).

* Donat, M., and Coauthors (2013). Updated analyses of temperature and precipitation extreme indices since the beginning of the twentieth century: The HadEX2 dataset. J.  Geophys. Res., 118, 2098–2118, doi:10.1002/jgrd.50150.

* Fouillet, A., Rey, G., Laurent, F., Pavillon, G. Bellec, S., Guihenneuc-Jouyaux, C., Clavel J., Jougla, E. and Hémon, D. (2006) Excess mortality related to the August 2003 heat wave in France. Int. Arch. Occup. Environ. Health, 80, 16–24.

* Gleason, K.L., J.H. Lawrimore, D.H. Levinson, T.R. Karl, and D.J. Karoly (2008). A Revised U.S. Climate Extremes Index. J. Climate, 21, 2124-2137

* Meehl, G. A., and Coauthors (2000). An introduction to trends inextreme weather and climate events: Observations, socio-economic impacts, terrestrial ecological impacts, and model projections. Bull. Amer. Meteor. Soc., 81, 413–416.

* Whitman, S., G. Good, E. R. Donoghue, N. Benbow, W. Y. Shou and S. X. Mou (1997). Mortality in Chicago attributed to the July 1995 heat wave. Amer. J. Public Health, 87, 1515–1518.

* Zhang, Y., M. Nitschke, and P. Bi (2013). Risk factors for direct heat-related hospitalization during the 2009 Adelaide heat-wave: A case crossover study. Sci. Total Environ., 442, 1–5.

* Zhang, X. , Alexander, L. , Hegerl, G. C., Jones, P. , Tank, A. K.,  Peterson, T. C., Trewin, B.  and Zwiers, F. W. (2011). Indices for  monitoring changes in extremes based on daily temperature and  precipitation data. WIREs Clim Change, 2: 851-870. doi:10.1002/wcc.147



Example plots
-------------

.. _fig_combinedindices1:
.. figure::  /recipes/figures/combined_climate_extreme_index/t90p_IPSL-CM5A-LR_rcp85_2020_2040.png
   :align:   center
   :width:   14cm



