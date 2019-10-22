.. _recipes_anav13jclim:

Land and ocean components of the global carbon cycle
====================================================

Overview
--------

This recipe reproduces most of the figures of `Anav et al. (2013)`_:

* Timeseries plot for different regions
* Seasonal cycle plot for different regions
* Errorbar plot for different regions showing mean and standard deviation
* Scatterplot for different regions showing mean vs. interannual variability
* 3D-scatterplot for different regions showing mean vs. linear trend and the
  model variability index (MVI) as a third dimension (color coded)
* Scatterplot for different regions comparing two variable against each other
  (*cSoil* vs. *cVeg*)

In addition, performance metrics are calculated for all variables using the
performance metric diagnostics (see details in :ref:`nml_perfmetrics`).


.. _`Anav et al. (2013)`: https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-12-00417.1


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_anav13jclim.yml


Diagnostics are stored in diag_scripts/

   * carbon_cycle/main.ncl
   * carbon_cycle/mvi.ncl
   * carbon_cycle/two_variables.ncl
   * perfmetrics/main.ncl
   * perfmetrics/collect.ncl


User settings in recipe
-----------------------

#. Preprocessor

   * ``mask_landsea``: Mask land/ocean.
   * ``regrid``: Regridding.
   * ``mask_fillvalues`` Mask common missing values on different datasets.

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
* *nbp_grid* (land, monthly, longitude, latitude, time)
* *gpp_grid* (land, monthly, longitude, latitude, time)
* *lai_grid* (land, monthly, longitude, latitude, time)
* *cveg_grid* (land, monthly, longitude, latitude, time)
* *csoil_grid* (land, monthly, longitude, latitude, time)
* *tos* (ocean, monthly, longitude, latitude, time)
* *fgco2_grid* (ocean, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

* CRU (*tas*, *pr*)
* JMA-TRANSCOM (*nbp_grid*, *fgco2_grid*)
* MTE (*gpp_grid*)
* LAI3g (*lai_grid*)
* NDP (*cveg_grid*)
* HWSD (*csoil_grid*)
* HadISST (*tos*)


References
----------

* Anav, A. et al.: Evaluating the land and ocean components of the global
  carbon cycle in the CMIP5 Earth System Models, J. Climate, 26, 6901-6843,
  doi: 10.1175/JCLI-D-12-00417.1, 2013.


Example plots
-------------

.. _fig_anav13jclim_1:
.. figure:: /recipes/figures/cox18nature/temperature_anomaly_HadCRUT4.png
   :align: center
   :width: 50%

   Simulated change in global temperature from CMIP5 models (coloured lines),
   compared to the global temperature anomaly from the HadCRUT4 dataset (black
   dots). The anomalies are relative to a baseline period of 1961–1990. The model
   lines are colour-coded, with lower-sensitivity models (λ > 1
   Wm\ :sup:`-2`\ K\ :sup:`-1`\ ) shown by green lines and higher-sensitivity
   models (λ < 1 Wm\ :sup:`-2`\ K\ :sup:`-1`\ ) shown by magenta lines.

.. _fig_anav13jclim_2:
.. figure:: /recipes/figures/cox18nature/emergent_relationship_HadCRUT4.png
   :align: center
   :width: 50%

   Emergent relationship between ECS and the ψ metric. The black dot-dashed
   line shows the best-fit linear regression across the model ensemble, with
   the prediction error for the fit given by the black dashed lines. The
   vertical blue lines show the observational constraint from the HadCRUT4
   observations: the mean (dot-dashed line) and the mean plus and minus one
   standard deviation (dashed lines).

.. _fig_anav13jclim_3:
.. figure:: /recipes/figures/cox18nature/pdf_HadCRUT4.png
   :align: center
   :width: 50%

   The PDF for ECS. The orange histograms (both panels) show the prior
   distributions that arise from equal weighting of the CMIP5 models in 0.5 K
   bins.
