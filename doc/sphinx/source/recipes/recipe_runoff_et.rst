Runoff_ET
=========

Overview
--------
A set of diagnostics to produce basin-scale comparisons of runoff (mrro), evapotranspiration
(evspsbl) and precipitation (pr) have been compiled in the runoff_et recipe. The recipe
calculates biases in the long-term annual means of the three variables for 12 large-scale
catchments areas on different continents and for different climates. For total runoff,
catchment averaged model values are compared to long-term averages of GRDC observations.
Due to the incompleteness of these station data, a year-to-year correspondence of data
cannot be achieved so only climatological data are considered, as in Hagemann et al. (2013).
Simulated precipitation is compared to catchment-averaged WATCH forcing data based on ERA-Interim
(WFDEI) data (Weedon et al., 2014) for the period 1979-2010. Here, the GPCC-corrected WFDEI
precipitation data are taken. Evapotranspiration observations are estimated using the
difference of the catchment-averaged WFDEI precipitation minus the climatological GRDC
river runoff.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_runoff_et.yml

Diagnostics are stored in diag_scripts/runoff_et/

    * catchment_analysis.py: bar and scatter plots for catchment averages of
      runoff, evapotranspiration and precipitation


User settings
-------------

#. runoff_et.yml

   *Required settings for script*

   * catchmentmask: netCDF file indicating the grid cell for a specific catchment. Modus of
     distribution not yet clearified. ESGF?

   *Optional settings for variables*

   * reference_dataset: dataset_name
     Datasets can be used as reference instead of defaults provided with the diagnostics.
     Must be identical for all variables.


Variables
---------

* evspsbl (atmos, monthly mean, time latitude longitude)
* pr      (atmos, monthly mean, time latitude longitude)
* mrro    (land,  monthly mean, time latitude longitude)


Observations and reformat scripts
---------------------------------

Default reference data based on GRDC and WFDEI are included in the diagnostic script
as catchment averages. They can be replaced with any gridded dataset by defining a
reference_dataset. A default catchment mask is (somehow) delivered together with
the diagnostic. All other datasets are remapped onto the catchment mask grid as part
of the diagnostics.


References
----------
* Duemenil Gates, L., S. Hagemann and C. Golz,
  Observed historical discharge data from major rivers for climate model validation.
  Max Planck Institute for Meteorology Report 307, Hamburg, Germany, 2000.

* Hagemann, S., A. Loew, A. Andersson,
  Combined evaluation of MPI-ESM land surface water and energy fluxes
  J. Adv. Model. Earth Syst., 5, doi:10.1029/2012MS000173, 2013.

* Weedon, G. P., G. Balsamo, N. Bellouin, S. Gomes, M. J. Best, and P. Viterbo,
  The WFDEI meteorological forcing data set: WATCH Forcing Data methodology applied
  to ERA‐Interim reanalysis data,
  Water Resour. Res., 50, 7505–7514, doi: 10.1002/2014WR015638, 2014


Example plots
-------------

.. _fig_runoff_et_1:
.. figure::  /recipes/figures/runoff_et/MPI-ESM-LR_historical_r1i1p1_bias-plot_mrro.png
   :align:   center
   :width:   14cm

   Barplot indicating the absolute and relative bias in annual runoff between MPI-ESM-LR (1970-2000)
   and long term GRDC data for specific catchments.

.. _fig_runoff_et_2:
.. figure::  /recipes/figures/runoff_et/MPI-ESM-LR_historical_r1i1p1_rocoef-vs-relprbias.png
   :align:   center
   :width:   14cm

   Biases in runoff coefficient (runoff/precipitation) and precipitation for major catchments of
   the globe. The MPI-ESM-LR historical simulation (1970-2000) is used as an example.

