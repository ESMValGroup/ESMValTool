.. _recipe_surface_trace_gas:

Surface trace gases
===================

Overview
--------

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_surface_trace_gas_ch4.yml
    * recipe_surface_trace_gas_co2.yml
    * recipe_surface_trace_gas_n2o.yml

Diagnostics are stored in esmvaltool/diag_scripts/surface_trace_gases/

    * diagnostic_surface_trace_gases.py : Plots the different surface trace gas evaluations.
    * utils_surface_trace_gases.py: Utility functions used for the surface trace gas evaluation routines.

This diagnostic allows to produce 3 different types of analysis, typically
comparing station data and model data. In the case of zonal plots, the latitude
ranges considered are 60N - 90N, 30N - 60N, 30S - 30N, 90S - 30S.

    * Seasonal maps of the trace gas surface concentration overlaid with
      station data points. A complimentary scatter plot compares the surface
      concentration at the station locations with the model data extracted at
      these same locations for each season.
    * Zonal plots of surface concentration time series and seasonal cycle.
      The left column includes monthly surface concentration time series.
      The center-left column includes yearly surface concentration scatter plots.
      The center-right column includes multi-annual monthly anomalies.
      The right column includes yearly min/max months time series.
      The plot is inspired from Figure 8 of Folberth et al. (2018).
    * Zonal plots of surface concentration amplitude and growth.
      The left column includes yearly concentration amplitude time series.
      The center column includes yearly concentration growth time series.
      The right column includes yearly sensitivities of amplitude/growth.
    * Taylor diagram of trace gas surface concentration for the observations
      and the different models from the recipe.

Each analysis can be turned on or off using the argument *plots* in the recipe.
Furthermore, one can change how the station data is pre-processed before
plotting by filtering stations with a minimum amount of data for months,
seasons, and years, using the *min_mon_per_seas*, *min_seas_per_year*, and
*min_seas_per_clim* arguments in the recipe (see next section).

User settings in recipe
-----------------------

#. Preprocessor

  * ``extract_surface_from_atm`` : Extract surface from the input 3D
    atmospheric variable based on the surface pressure *ps*.

#. Script diagnostic_surface_trace_gases.py

  *Required settings for script*
  * trace_gas: The surface trace gas for which to run the evaluation.
  * min_mon_per_seas: The minimum number of months used to calculate each seasonal mean. This must be between 1 and 3. Recommended value is 3.
  * min_seas_per_year: The minimum number of seasonal means in each year. This must be between 1 and 4. Recommended value is 4.
  * min_seas_per_clim: The minimum number of seasonal means used to calculate the multiannual seasonal mean. This must be between 1 and the number of years of available NOAA GML Surface Flask data.

  *Optional settings for script*

  * plots: List of the diagnostics to plot among *seas_maps*, *timeserie_lat*,
    *sensitivity_ampl_trend* and *taylor_diag*. Default is all of them.

  *Required settings for variables*

  * None

  *Optional settings for variables*

  * None

  *Required settings for preprocessor*

  * None

  *Optional settings for preprocessor*

  * None


Variables
---------

The trace gas depends on the recipe used and can be one of the following:

* *ch4* (atmos, monthly mean, height longitude latitude time)
* *co2* (atmos, monthly mean, height longitude latitude time)
* *n2o* (atmos, monthly mean, height longitude latitude time)
* *ps* (atmos, monthly mean, longitude latitude time) as a supplementary
  variable for the `extract_surface_from_atm` preprocessor.

Observations and reformat scripts
---------------------------------

* The NOAA GML Surface Flask data is downloaded from the NOAA GML website
  using the downloaders:

  .. code-block:: yaml

        $ esmvaltool data download NOAA-GML-SURFACE-FLASK-CH4.
        $ esmvaltool data download NOAA-GML-SURFACE-FLASK-CO2.
        $ esmvaltool data download NOAA-GML-SURFACE-FLASK-N2O.

* The NOAA GML Surface Flask data is formatted using the formatters:

  .. code-block:: yaml

        $ esmvaltool data format NOAA-GML-SURFACE-FLASK-CH4.
        $ esmvaltool data format NOAA-GML-SURFACE-FLASK-CO2.
        $ esmvaltool data format NOAA-GML-SURFACE-FLASK-N2O.

References
----------
* Folberth et al.: Description and Evaluation of an Emission-Driven and Fully Coupled Methane Cycle in UKESM1, 10.1029/2021MS002982, 2018.

Example plots
-------------

.. _fig_surface_trace_gas_1:
.. figure::  /recipes/figures/surface_trace_gas/CNRM-ESM2-1_Amon_esm-hist_co2s_2000_2014_seas_map.png
   :align:   center

   Evaluation of seasonal surface concentration of CO2 from CNRM-ESM2-1 esm-hist member r1i1p1f3 against the NOAA GML climatology from ground-based observations. The multiannual seasonal mean is calculated for the model data for the period 2000-2014. The model output is overlaid with the observational climatology.

.. _fig_surface_trace_gas_2:
.. figure::  /recipes/figures/surface_trace_gas/CNRM-ESM2-1_Amon_esm-hist_co2s_2000_2014_scatter.png
   :align:   center

   Evaluation of seasonal surface concentration of CO2 from CNRM-ESM2-1 esm-hist member r1i1p1f3 against the NOAA GML climatology from ground-based observations. The multiannual seasonal mean is calculated for the model data for the period 2000-2014.

.. _fig_surface_trace_gas_3:
.. figure::  /recipes/figures/surface_trace_gas/CNRM-ESM2-1_Amon_esm-hist_co2s_2000_2014_timeseries_latitude.png
   :align:   center

   Evaluation of surface concentration time series (monthly, seasonal, annual) of CO2 from CNRM-ESM2-1 esm-hist member r1i1p1f3 against the NOAA GML climatology from ground-based observations. The multiannual seasonal mean is calculated for the model data for the period 2000-2014.

.. _fig_surface_trace_gas_4:
.. figure::  /recipes/figures/surface_trace_gas/CNRM-ESM2-1_Amon_esm-hist_co2s_2000_2014_sensitivity_ampl_growth.png
   :align:   center

   Evaluation of surface concentration time series (amplitude, growth, sensitivity) of CO2 from CNRM-ESM2-1 esm-hist member r1i1p1f3 against the NOAA GML climatology from ground-based observations. The multiannual seasonal mean is calculated for the model data for the period 2000-2014.

.. _fig_surface_trace_gas_5:
.. figure::  /recipes/figures/surface_trace_gas/trace_gas_co2_CNRM-ESM2-1_2000_2014_taylor_diag.png
   :align:   center

   Taylor diagram of surface concentration of CO2 from CNRM-ESM2-1 esm-hist member r1i1p1f3 against the NOAA GML climatology from ground-based observations. The multiannual seasonal mean is calculated for the model data for the period 2000-2014.
