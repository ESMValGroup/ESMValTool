.. _recipes_ipccwg1ar6ch3:

IPCC AR6 Chapter 3 (selected figures)
=====================================

Overview
--------

The goal of this recipe is to collect diagnostics to reproduce Chapter 3 of AR6,
so that the plots can be readily reproduced and compared to previous CMIP
versions. In this way we can next time start with what was available in the
previous round and can focus on developing more innovative methods of analysis
rather than constantly having to "re-invent the wheel".

The plots are produced collecting the diagnostics from individual recipes. The
following figures from Eyring et al. (2021) can currently be reproduced:

    * Figure 3.3 a,b,c,d: Surface Air Temperature - Model Bias

    * Figure 3.4: Anomaly Of Near-Surface Air Temperature

    * Figure 3.5: Temporal Variability Of Near-Surface Air Temperature

    * Figure 3.13: Precipitation - Model Bias

    * Figure 3.15: Precipitation Anomaly

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/ipccwg1ar6ch3/

    * recipe_ipccwg1ar6ch3_atmosphere.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    Fig. 3.3:
    * ipcc_ar5/ch12_calc_IAV_for_stippandhatch.ncl: See :ref:`here:<ch12_calc_IAV_for_stippandhatch.ncl>`.
    * ipcc_ar6/model_bias.ncl

    Fig. 3.4:
    * ipcc_ar6/tas_anom.ncl
    * ipcc_ar6/tsline_collect.ncl

    Fig. 3.5:
    * ipcc_ar6/zonal_st_dev.ncl

    Fig. 3.13:
    * ipcc_ar5/ch12_calc_IAV_for_stippandhatch.ncl: See :ref:`here:<ch12_calc_IAV_for_stippandhatch.ncl>`.
    * ipcc_ar6/model_bias.ncl

    Fig. 3.15:
    * ipcc_ar6/precip_anom.ncl

User settings in recipe
-----------------------

#. Script ipcc_ar5/ch12_calc_IAV_for_stippandhatch.ncl

   See :ref:`here<ch12_calc_IAV_for_stippandhatch.ncl>`.

#. Script ipcc_ar6/model_bias.ncl

   *Optional settings (scripts)*

   * plot_abs_diff: additionally also plot absolute differences (true, false)
   * plot_rel_diff: additionally also plot relative differences (true, false)
   * plot_rms_diff: additionally also plot root mean square differences (true, false)
   * projection: map projection, e.g., Mollweide, Mercator
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)

   *Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Color tables*

   * variable "tas" and "tos":
     diag_scripts/shared/plot/rgb/ipcc-ar6_temperature_div.rgb,
     diag_scripts/shared/plot/rgb/ipcc-ar6_temperature_10.rgb,
     diag_scripts/shared/plot/rgb/ipcc-ar6_temperature_seq.rgb
   * variable "pr": diag_scripts/shared/plots/rgb/ipcc-ar6_precipitation_seq.rgb,
     diag_scripts/shared/plot/rgb/ipcc-ar6_precipitation_10.rgb
   * variable "sos": diag_scripts/shared/plot/rgb/ipcc-ar6_misc_seq_1.rgb,
     diag_scripts/shared/plot/rgb/ipcc-ar6_misc_div.rgb

#. Script ipcc_ar6/tas_anom.ncl

   *Required settings for script*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings for script*

   * blending: if true, calculates blended surface temperature
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_value: if true, right panel with mean values is attached
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * y-min: set min of y-axis
   * y-max: set max of y-axis
   * header: if true, region name as header
   * volcanoes: if true, adds volcanoes to the plot
   * write_stat: if true, write multi model statistics in nc-file

   *Optional settings for variables*

   * reference_dataset: reference dataset; REQUIRED when calculating
     anomalies

   *Color tables*

   * e.g. diag_scripts/shared/plot/styles/cmip5.style

#. Script ipcc_ar6/tsline_collect.ncl

   *Optional settings for script*

   * blending: if true, then var="gmst" otherwise "gsat"
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * region: name of domain
   * plot_units: variable unit for plotting
   * y-min: set min of y-axis
   * y-max: set max of y-axis
   * order: order in which experiments should be plotted 
   * stat_shading: if true: shading of statistic range
   * ref_shading: if true: shading of reference period

   *Optional settings for variables*

   * reference_dataset: reference dataset; REQUIRED when calculating
     anomalies

#. Script ipcc_ar6/zonal_st_dev.ncl

   *Required settings for script*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings for script*

   * plot_legend: if true, plot legend will be plotted
   * plot_units: variable unit for plotting
   * multi_model_mean: if true, multi-model mean and uncertaintiy will be 
     plotted

   *Optional settings for variables*

   * reference_dataset: reference dataset; REQUIRED when calculating
     anomalies

--> to continue

#. Script ipcc_ar6/precip_anom.ncl

   *Required settings for script*

   * panels: list of variables plotted in each panel
   * start_year: start of time coordinate
   * end_year: end of time coordinate

   *Optional settings for script*

   * anomaly: true if anomaly should be calculated
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * header: if true, region name as header
   * stat: statistics for multi model nc-file (MinMax,5-95,10-90)
   * y_min: set min of y-axis
   * y_max: set max of y-axis

Variables
---------

* pr (atmos, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* tasa (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* BerkeleyEarth (tasa - esmvaltool/cmorizers/obs/cmorize_obs_berkeleyearth.py)
* CRU (pr - esmvaltool/cmorizers/obs/cmorize_obs_cru.py)
* ERA5 (tas - ERA5 data can be used via the native6 project)
* GHCN (pr - esmvaltool/cmorizers/obs/cmorize_obs_ghcn.ncl)
* GPCP-SG (pr - obs4MIPs)
* HadCRUT5 (tasa - esmvaltool/cmorizers/obs/cmorize_obs_hadcrut5.py)
* Kadow (tasa - esmvaltool/cmorizers/obs/cmorize_obs_kadow.ncl)
* NOAAGlobalTemp (tasa - esmvaltool/cmorizers/obs/cmorize_obs_noaaglobaltemp.ncl)


References
----------

* Eyring, V., N.P. Gillett, K.M. Achuta Rao, R. Barimalala, M. Barreiro
  Parrillo, N. Bellouin, C. Cassou, P.J. Durack, Y. Kosaka, S. McGregor, S. Min,
  O. Morgenstern, and Y. Sun, 2021: Human Influence on the Climate System. In 
  Climate Change 2021: The Physical Science Basis. Contribution of Working Group
  I to the Sixth Assessment Report of the Intergovernmental Panel on Climate 
  Change [Masson-Delmotte, V., P. Zhai, A. Pirani, S.L. Connors, C. Péan, S. 
  Berger, N. Caud, Y. Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. 
  Lonnoy, J.B.R. Matthews, T.K.  Maycock, T. Waterfield, O. Yelekçi, R. Yu, and 
  B. Zhou (eds.)]. Cambridge University Press. In Press.


Example plots
-------------

.. figure::  /recipes/figures/flato13ipcc/fig-9-2.png
   :align:   center

   Figure 9.2 a,b,c: Annual-mean surface air temperature for the period
   1980-2005. a) multi-model mean, b) bias as the difference between the
   CMIP5 multi-model mean and the climatology from ERA-Interim
   (Dee et al., 2011), c) mean absolute model error with respect to the
   climatology from ERA-Interim.

.. figure::  /recipes/figures/seaice/trend_sic_extend_Arctic_September_histogram.png
   :align:   center
   :width:   9cm

   Figure 9.24c: Sea ice extent trend distribution for the Arctic in September.
