.. _recipes_recipe_calculate_gwl_exceedance_stats:

Physical Climate at Global Warming Levels (GWLs)
==================================================

Overview
--------

This recipe calculates years of Global Warming Level (GWL) exceedances in CMIP
models as described in `Swaminathan et al (2022)`. Time series of the
anomalies in annual global mean surface air temperature (GSAT)
are calculated with respect to the 1850-1900 time-mean of each
individual time series. To limit the influence of short-term variability,
a 21-year centered running mean is applied to the time series. The year at
which the time series exceeds warming levels or temperatures such as 1.5C
is then recorded for the specific model ensemble member and future scenario.
Once the years of exceedance are calculated, the time averaged global
mean and standard deviation for the multimodel ensemble
over the 21-year period around the year of exceedance are plotted.
By selecting specific scenarios, the multimodel mean and spread for a single
future scenario can be plotted as shown in the examples below.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_calculate_gwl_exceedance_stats.yml

Diagnostics are stored in esmvaltool/diag_scripts/

* gwls/calculate_gwl_exceedance_years.py
* gwls/plot_gwl_exceedance_mm_stats.py



User settings in recipe
-----------------------

* Preprocessors

   * ``calculate_anomalies``

      * ``custom_order`` (*bool: True*) : Follow order of preprocessors as given.
      * ``area_statistics`` (*operator: mean*) : Calculate area averaged means.
      * ``annual_statistics`` (*operator: mean*): Calculate annual means.
      * ``anomalies``: Calculate anomalies for the full time period and with
         respect to the reference period.
      * ``extract_time`` : Extracting time period to calculate time series for.

   * ``multi_model_gwl_stats``

      * ``custom_order`` (*bool: True*) : Follow order of preprocessors as given.
      * ``extract_time`` : Extract time for the period of the time series calculation.
      * ``annual_statistics``  (*operator: mean*): Calculate annual means.
      * ``regrid``: Regrid to a common resolution so multimodel means can be calculated.


* Script calculate_gwl_exceedance_years.py

   * ``window_size`` : Number of years to average over to smooth the time series.
   * ``gwls``: Global warming levels for which years of exceedances are to be calculated.

* Script plot_gwl_exceedance_mm_stats.py

   * ``ancestors`` : Output file from the GWL exceedance calculation step and
      preprocessed output for the physical variable.
   * ``quickplot`` : Plotting options.

      * ``plot_type`` : For recording provenance information on type of plot.
      * ``cmap_mean``: Colormap for the mean map contour plot.
      * ``cmap_stdev``: Colormap for the standard deviation map contour plot.
      * ``title_var`` : Variable name for the plot title.
      * ``mean_level_params``: Start, end and step size for the levels in
        the contour map plot. Values are floats and used by the numpy arange
        function.
      * ``stdev_level_params``: As above but for the standard deviation plots.


Variables
---------

* tas (atmos, monthly mean, latitude, longitude, time)
* pr (atmos, monthly mean, latitude, longitude, time)

Observations and reformat scripts
---------------------------------

* None used but can be included.


References
----------

* Swaminathan, R., R. J. Parker, C. G. Jones, R. P. Allan, T. Quaife,
  D. I. Kelley, L. de Mora, and J. Walton, 2022: The Physical Climate
  at Global Warming Thresholds as Seen in the U.K. Earth System Model.
  J. Climate, 35, 29–48, https://doi.org/10.1175/JCLI-D-21-0234.1.

Example plots
-------------

.. _fig_calculate_gwl_exceedance_stats_1:
.. figure:: /recipes/figures/gwls/CMIP6_mm_mean_ssp126_1.5_tas.png
   :align:   left
   :width: 50%

   Multimodel mean of temperature under SSP1-2.6 at 1.5 degC warming.

.. _fig_calculate_gwl_exceedance_stats_2:
.. figure:: /recipes/figures/gwls/CMIP6_mm_stdev_ssp126_1.5_tas.png
   :align:   left
   :width: 50%

   Multimodel standard deviation  of temperature under SSP1-2.6 at 1.5 degC warming.
