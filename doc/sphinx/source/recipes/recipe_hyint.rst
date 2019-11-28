.. _recipes_hyint:

Hydroclimatic intensity and extremes (HyInt)
============================================


Overview
--------
The Earth’s hydrological cycle is of key importance both for the climate system and society. For example, the intensity and distribution of precipitation determines the availability or scarcity of fresh water in a certain region, and it is also related to the severity of hazardous events such as flooding or droughts. The simple investigation of average precipitation quantities can clearly hide some of the most relevant aspects of the hydrological cycle and its extremes (e.g., Giorgi et al., 2014). More in general, temperature and precipitation extremes have been the focus of recent climate studies attempting to capture the most relevant component of climate variability and impact on society in a changing climate (e.g., Alexander, 2016. A particular effort has been dedicated to developing and standardising indices that can be adopted for investigation studies with observations and climate models. This tool was developed to calculate a number of hydroclimatic and climate extremes indices and allow a multi-index evaluation of climate models. The tool firstly computes a set of 6 indices that allow to evaluate the response of the hydrological cycle to global warming with a joint view of both wet and dry extremes. The indices were selected following Giorgi et al. (2014) and include the simple precipitation intensity index (SDII), the maximum dry spell length (DSL) and wet spell length (WSL), the hydroclimatic intensity index (HY-INT), which is a measure of the overall behaviour of the hydroclimatic cycle (Giorgi et al., 2011), and the precipitation area (PA), i.e. the area over which at any given day precipitation occurs, (Giorgi et al., 2014). Secondly, also a selection of the 27 temperature and precipitation -based indices of extremes from the Expert Team on Climate Change Detection and Indices (ETCCDI) produced by the climdex (https://www.climdex.org) library can be ingested to produce a multi-index analysis. The tool allows then to perform a subsequent analysis of the selected indices calculating timeseries and trends over predefined continental areas, normalized to a reference period. Trends are calculated using the R `lm` function and significance testing performed with a Student T test on non-null coefficients hypothesis. Trend coefficients are stored together with their statistics which include standard error, t value and Pr(>|t|). The tool can then produce a variety of types of plots including global and regional maps, maps of comparison between models and a reference dataset, timeseries with their spread, trend lines and summary plots of trend coefficients.

The hydroclimatic indices calculated by the diagnostic and included in the output are defined as follows:

* PRY = mean annual precipitation
* INT = mean annual precipitation intensity (intensity during wet days, or simple precipitation intensity index SDII)
* WSL = mean annual wet spell length (number of consecutive days during each wet spell)
* DSL = mean annual dry spell lenght (number of consecutive days during each dry spell)
* PA  = precipitation area (area over which of any given day precipitation occurs)
* R95 = heavy precipitation index (percent of total precipitation above the 95% percentile of the reference distribution)
* HY-INT = hydroclimatic intensity. HY-INT = normalized(INT) x normalized(DSL).



Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_hyint.yml (evaluating the 6 hydroclimatic indices)

Diagnostics are stored in diag_scripts/hyint/

* hyint.R

and subroutines

* hyint_diagnostic.R
* hyint_functions.R
* hyint_parameters.R
* hyint_plot_trends.R
* hyint_etccdi_preproc.R
* hyint_metadata.R
* hyint_plot_maps.R
* hyint_preproc.R
* hyint_trends.R


User settings
-------------

*Required settings for script*


* norm_years: first and last year of reference normalization period to be used for normalized indices

* select_indices: indices to be analysed and plotted. Select one or more fields from the following list (order-sensitive): "pa_norm", "hyint",  "int_norm", "r95_norm", "wsl_norm", "dsl_norm", "int", "dsl", "wsl"

* select_regions: Select regions for timeseries and maps from the following list: GL=Globe, GL60=Global 60S/60N, TR=Tropics (30S/30N), SA=South America, AF=Africa, NA=North America, IN=India, EU=Europe, EA=East-Asia, AU=Australia

* plot_type: type of figures to be plotted. Select one or more from: 1=lon/lat maps per individual field/exp/multi-year mean, 2=lon/lat maps per individual field exp-ref-diff/multi-year mean, 3=lon/lat maps multi-field/exp-ref-diff/multi-year mean, 11=timeseries over required individual region/exp, 12=timeseries over multiple regions/exp, 13=timeseries with multiple models, 14=summary trend coefficients multiple regions, 15=summary trend coefficients multiple models

*Optional settings for script (with default setting)*

#. Data

   * rgrid (false): Define whether model data should be regridded. (a) false to keep original resolution; (b) set desired regridding resolution in cdo format e.g., "r320x160"; (c) "REF" to use resolution of reference model

#. Plotting

   * npancol (2): number of columns in timeseries/trends multipanel figures
   * npanrow (3): number of rows in timeseries/trends multipanel figures
   * autolevels (true): select automated (true) or pre-set (false) range of values in plots
   * autolevels_scale (1): factor multiplying automated range for maps and timeseries
   * autolevels_scale_t (1.5): factor multiplying automated range for trend coefficients

#. Maps

   * oplot_grid (false): plot grid points over maps
   * boxregion (false): !=0 plot region boxes over global maps with thickness = abs(boxregion); white (>0) or grey (<0).
   * removedesert (false) remove (flag as NA) grid points with mean annual pr < 0.5 mm/day (deserts, Giorgi2014). This affects timeseries and trends calculations too.

#. Timeseries and trends

   * weight_tseries (true): adopt area weights in timeseries
   * trend_years (false): (a) false = apply trend to all years in dataset; (b) [year1, year2] to apply trend calculation and plotting only to a limited time interval
   * add_trend (true): add linear trend to plot
   * add_trend_sd (false): add dashed lines of stdev range to timeseries
   * add_trend_sd_shade (false): add shade of stdev range to timeseries
   * add_tseries_lines (true): plot lines connecting timeseries points
   * add_zeroline (true): plot a dashed line at y=0
   * trend_years_only (false): limit timeseries plotting to the time interval adopted for trend calculation (excluding the normalization period)
   * scale100years (true): plot trends scaled as 1/100 years
   * scalepercent (false): plot trends as percent change


Variables
---------

* pr (atmos, daily mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None.


References
----------

* Giorgi et al., 2014, J. Geophys. Res. Atmos., 119, 11,695–11,708, doi:10.1002/ 2014JD022238
* Giorgi et al., 2011, J. Climate 24, 5309-5324, doi:10.1175/2011JCLI3979.1


Example plots
-------------

.. figure:: /recipes/figures/hyint/hyint_maps.png
   :width: 10cm

   Mean hydroclimatic intensity for the EC-EARTH model, for the historical + RCP8.5 projection in the period 1976-2099

.. figure:: /recipes/figures/hyint/hyint_timeseries.png
   :width: 12cm

   Timeseries for multiple indices and regions for the ACCESS1-0 model, for the historical + RCP8.5 projection in the period 1976-2099, normalized to the 1976-2005 historical period.

.. figure:: /recipes/figures/hyint/hyint_trends.png
   :width: 12cm

   Multi-model trend coefficients over selected indices for CMIP5 models in the RCP8.5 2006-2099 projection, normalized to the 1976-2005 historical period.
