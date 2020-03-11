.. _recipe_carvalhais2014nat:

Turnover time of carbon over land
====================================================

Overview
--------

This recipe evaluates the turnover time of carbon over land following `Carvalhais et al. (2014)`_:

* Maps of the turnover time of carbon
* Model-to-observation full factorial comparison - with global maps and density scatterplots
* Zonal distributions of turnover time
* Zonal correlations of turnover time and climate

In addition, metrics of global turnover times and correlations with observation are calculated and plotted over the figure.


.. _tau calculation:

Calculation of turnover time
----------------------------

The turnover time of carbon is defined as,
.. math::

   \tau =  \frac{cSoil + cVeg}{gpp}

where :math:`cSoil` and :math:`cVeg` are the carbon content of the soil and vegetation over the land surface, and gpp is the gross primary productivity.

This equation assumes that the total carbon content over land is the sum of the soil and vegetation components, which may be in contrast to what the original publication used. In the original publication, all the carbon pools which respired to the atmosphere were added up to calculate the total carbon stock.

Due to inherent dependence of the diagnostic on the uncertainty of the observation, the recipe is strictly dependent on the observation provided through the data portal of the Max Planck Inevertheless,
this script provides two configuration options to avoid high MVI values, but
they are not related to the original paper or any other peer-revied study and
should be used with great caution (see :ref:`user settings`).

.. _`Anav et al. (2013)`: https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-12-00417.1


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * _recipe_carvalhais2014nat.yml


Diagnostics are stored in diag_scripts/

   * land_carbon_cycle/diag_global_turnover.py
   * land_carbon_cycle/diag_zonal_turnover.py
   * land_carbon_cycle/diag_zonal_correlation.py
   * perfmetrics/main.ncl
   * perfmetrics/collect.ncl


.. _user settings:

User settings in recipe
-----------------------

#. Preprocessor

   * ``mask_fillvalues``: Mask common missing values on different datasets.
   * ``mask_landsea``: Mask land/ocean.
   * ``regrid``: Regridding.
   * ``weighting_landsea_fraction``: Land/ocean fraction weighting.

#. Script carbon_cycle/main.ncl

   * ``region``, *str*: Region to be averaged.
   * ``legend_outside``, *bool*: Plot legend in a separate file (does not
     affect errorbar plot and evolution plot)
   * ``seasonal_cycle_plot``, *bool*: Draw seasonal cycle plot.
   * ``errorbar_plot``, *bool*: Draw errorbar plot.
   * ``mean_IAV_plot``, *bool*: Draw Mean (x-axis), IAV (y-axis) plot.
   * ``evolution_plot``, *bool*: Draw time evolution of a variable comparing
     a reference dataset to multi-dataset mean; requires ref_dataset in recipe.
   * ``sort``, *bool*, optional (default: ``False``): Sort dataset in
     alphabetical order.
   * ``anav_month``, *bool*, optional (default: ``False``): Conversion of
     y-axis to PgC/month instead of /year.
   * ``evolution_plot_ref_dataset``, *str*, optional: Reference dataset for
     evolution_plot. Required when ``evolution_plot`` is ``True``.
   * ``evolution_plot_anomaly``, *str*, optional (default: ``False``): Plot
     anomalies in evolution plot.
   * ``evolution_plot_ignore``, *list*, optional: Datasets to ignore in
     evolution plot.
   * ``evolution_plot_volcanoes``, *bool*, optional (default: ``False``): Turns
     on/off lines of volcano eruptions in evolution plot.
   * ``evolution_plot_color``, *int*, optional (default: ``0``): Hue of the
     contours in the evolution plot.

#. Script carbon_cycle/mvi.ncl

   * ``region``, *str*: Region to be averaged.
   * ``reference_dataset``, *str*: Reference dataset for the MVI calculation
     specified for each variable seperately.
   * ``mean_time_range``, *list*, optional: Time period over which the mean is
     calculated (if not given, use whole time span).
   * ``trend_time_range``, *list*, optional: Time period over which the trend
     is calculated (if not given, use whole time span).
   * ``mvi_time_range``, *list*, optional: Time period over which the MVI is
     calculated (if not given, use whole time span).
   * ``stddev_threshold``, *float*, optional (default: ``1e-2``): Threshold to
     ignore low standard deviations (relative to the mean) in the MVI
     calculations. See also :ref:`mvi calculation`.
   * ``mask_below``, *float*, optional: Threshold to mask low absolute values
     (relative to the mean) in the input data (not used by default). See also
     :ref:`mvi calculation`.

#. Script carbon_cycle/two_variables.ncl

   * ``region``, *str*: Region to be averaged.

#. Script perfmetrics/main.ncl

   See :ref:`nml_perfmetrics`.

#. Script perfmetrics/collect.ncl

   See :ref:`nml_perfmetrics`.


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *gpp* (land, monthly, longitude, latitude, time)
* *cVeg* (land, monthly, longitude, latitude, time)
* *cSoil* (land, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* CRU (*tas*, *pr*)
* JMA-TRANSCOM (*nbp*, *fgco2*)
* MTE (*gpp*)
* LAI3g (*lai*)
* NDP (*cveg*)
* HWSD (*csoil*)
* HadISST (*tos*)


References
----------

* Anav, A. et al.: Evaluating the land and ocean components of the global
  carbon cycle in the CMIP5 Earth System Models, J. Climate, 26, 6901-6843,
  doi: 10.1175/JCLI-D-12-00417.1, 2013.


Example plots
-------------

.. _fig_carvalhais2014nat_1:
.. figure:: /recipes/figures/carvalhais2014nat/comparison_zonal_pearson_correlation_turnovertime_climate_Carvalhais2014.png
   :align: center
   :width: 80%

   Time series of global net biome productivity (NBP) over the period
   1901-2005. Similar to Anav et al.  (2013), Figure 5.

.. _fig_carvalhais2014nat_2:
.. figure:: /recipes/figures/carvalhais2014nat/global_comparison_matrix_models_Carvalhais2014.png
   :align: center
   :width: 80%

   Seasonal cycle plot for nothern hemisphere gross primary production (GPP)
   over the period 1986-2005. Similar to Anav et al. (2013), Figure 9.
