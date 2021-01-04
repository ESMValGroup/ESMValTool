.. _recipes_bock20jgr:

Quantifying progress across different CMIP phases
=================================================

Overview
--------

The recipe recipe_bock20jgr.yml generates figures to quantify the progress across
different CMIP phases.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_bock20jgr.yml

Diagnostics are stored in diag_scripts/bock20jgr/

    * tsline.ncl: timeseries of global mean surface temperature anomaly
    * tsline_collect.ncl: collect different timeseries from tsline.ncl to 
      compare different models ensembels
    * model_bias.ncl: global maps of the multi-model mean and the multi-model
      mean bias


User settings in recipe
-----------------------

#. Script tsline.ncl

   *Required settings (scripts)*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings (scripts)*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined period;
     for each gird point by removing the mean for the given
     calendar month (requiring at least 50% of the data to be
     non-missing)
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_value: if true, right panel with mean values is attached
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * y_min: set min of y-axis
   * y_max: set max of y-axis
   * mean_nh_sh: if true, calculate first NH and SH mean
   * volcanoes: if true, lines of main volcanic eruptions will be added
   * header: if true, region name as header
   * write_stat: if true, write multi model statistics in nc-file

   *Required settings (variables)*

   none

   * Optional settings (variables)

   none

#. Script tsline_collect.ncl

   *Required settings (scripts)*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings (scripts)*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined period
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * region: name of domain
   * plot_units: variable unit for plotting
   * y_min: set min of y-axis
   * y_max: set max of y-axis
   * order: order in which experiments should be plotted
   * header: if true, region name as header
   * stat_shading: if true: shading of statistic range
   * ref_shading: if true: shading of reference period


   *Required settings (variables)*

   none

   * Optional settings (variables)

   none

#. Script model_bias.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * projection:    map projection, e.g., Mollweide, Mercator
   * timemean:      time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean) 

   * Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Optional settings (variables)*

   * long_name: description of variable

   *Color tables*

   * variable "tas": diag_scripts/shared/plot/rgb/ipcc-ar6_temperature_div.rgb,
   * variable "pr-mmday": diag_scripts/shared/plots/rgb/ipcc-ar6_precipitation_seq.rgb
     diag_scripts/shared/plot/rgb/ipcc-ar6_precipitation_div.rgb


Variables
---------

* tas (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* HadCRUT4 - surface temperature anomalies

  *Reformat script:* cmorizers/obs/cmorize_obs_hadcrut4.ncl

* ERA5 - reanalysis of surface temperature

  *Reformat script:* recipes/cmorizers/recipe_era5.yml

* GPCP-SG (obs4mips) - Global Precipitation Climatology Project total
  precipitation

References
----------

* Bock, L., Lauer, A., Schlund, M., Barreiro, M., Bellouin, N., Jones, C., 
  Predoi, V., Meehl, G., Roberts, M., and Eyring, V.: Quantifying progress 
  across different CMIP phases with the ESMValTool, Journal of Geophysical 
  Research: Atmospheres, 125, e2019JD032321. https://doi.org/10.1029/2019JD032321

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.


Example plots
-------------

.. _fig_1:
.. figure::  /recipes/figures/clouds/liq_h2o_path_multi.png
   :align:   center

   The 20-yr average LWP (1986-2005) from the CMIP5 historical model runs and
   the multi-model mean in comparison with the UWisc satellite climatology
   (1988-2007) based on SSM/I, TMI, and AMSR-E (O'Dell et al. 2008).

.. _fig_2:
.. figure::  /recipes/figures/clouds/liq_h2o_taylor.png
   :align:   center
   :width:   7cm

   Taylor diagram showing the 20-yr annual average performance of CMIP5 models
   for total cloud fraction as compared to MODIS satellite observations.

.. _fig_3:
.. figure::  /recipes/figures/clouds/cloud_sweffect.png
   :align:   center
   :width:   9cm

.. _fig_4:
.. figure::  /recipes/figures/clouds/cloud_lweffect.png
   :align:   center
   :width:   9cm

