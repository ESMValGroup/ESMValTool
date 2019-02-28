.. _nml_perfmetrics:

Performance metrics for essential climate parameters
====================================================

Overview
--------

The goal is to create a standard recipe for the calculation of performance metrics to quantify the ability of the models to reproduce the climatological mean annual cycle for selected "Essential Climate Variables" (ECVs) plus some additional corresponding diagnostics and plots to better understand and interpret the results. The recipe can be used to calculate performance metrics at different vertical levels (e.g., 5, 30, 200, 850 hPa as in Gleckler et al., 2008) and in four regions (global, tropics 20°N-20°S, northern extratropics 20°-90°N, southern extratropics 20°-90°S). As an additional reference, we consider the Righi et al. (2015) paper.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_perfmetrics_CMIP5.yml

Diagnostics are stored in diag_scripts/perfmetrics/

* main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields. The calculated fields can also be plotted as difference w.r.t. a given reference dataset. main.ncl also calculates RMSD, bias and taylor metrics. Input data have to be regridded to a common grid in the preprocessor. Each plot type is created by a separated routine, as detailed below.
* cycle.ncl: creates an annual/seasonal cycle plot.
* zonal.ncl: creates a zonal (lat-pressure) plot.
* latlon.ncl: creates a lat-lon plot.
* cycle_latlon.ncl: precalculates the metrics for a time-lat-lon field, with different options for normalization.
* collect.ncl: collects and plots the metrics previously calculated by cycle_latlon.ncl.

User settings
-------------

#. main.ncl

   *Required settings for script*

   * plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon), cycle_latlon (time, lat, lon)
   * time_avg: type of time average (opt argument of time_operations in diag_scripts/shared/statistics.ncl)
   * region: selected region (see select_region in diag_scripts/shared/latlon.ncl)
   
   *Optional settings for script*
   
   * styleset: for plot_type cycle only (as in diag_scripts/shared/plot/styles/)
   * plot_stddev: for plot_type cycle only, plots standard deviation as shading
   * legend_outside: for plot_type cycle only, plots the legend in a separate file
   * t_test: for plot_type zonal or latlon, calculates t-test in difference plots (default: False)
   * conf_level: for plot_type zonal or latlon, adds the confidence level for the t-test to the plot (default: False)
   * projection: map projection for plot_type latlon (default: CylindricalEquidistant)
   * draw_plots: draws plots (default: True)
   * plot_diff: draws difference plots (default: False)
   * calc_grading: calculates grading metrics (default: False)
   * stippling: uses stippling to mark statistically significant differences (default: False = mask out non-significant differences in gray)
   * show_global_avg: diplays the global avaerage of the input field as string at the top-right of lat-lon plots (default: False)
   * metric: chosen grading metric(s) (if calc_grading is True)
   * normalization: metric normalization (for RMSD and BIAS metrics only)
   * abs_levs: list of contour levels for absolute plot
   * diff_levs: list of contour levels for difference plot
   * zonal_cmap: for plot_type zonal only, chosen color table (default: "amwg_blueyellowred")
   * zonal_ymin: for plot_type zonal only, minimum pressure level on the y-axis (default: 5. hPa)
   * latlon_cmap: for plot_type latlon only, chosen color table (default: "amwg_blueyellowred")
   * plot_units: plotting units (if different from standard CMOR units)
   
   *Required variable_info attributes*
   
   * reference_dataset: reference dataset to compare with (usually the observations).
   
   *Optional variable_info attributes*

   * alternative_dataset: a second dataset to compare with.

These settings are passed to the other scripts by main.ncl, depending on the selected plot_type.

#. collect.ncl

   *Required settings for scripts*

   * metric: selected metric (RMSD, BIAS or taylor)
   * label_bounds: for RMSD and BIAS metrics, min and max of the labelbar
   * label_scale: for RMSD and BIAS metrics, bin width of the labelbar
   * colormap: for RMSD and BIAS metrics, color table of the labelbar
   
   *Optional settings for script*
   
   * label_lo: adds lower triange for values outside range
   * label_hi: adds upper triange for values outside range
   * cm_interval: min and max color of the color table
   * cm_reverse: reverses the color table
   * sort: sorts datasets in alphabetic order (excluding MMM)
   * title: plots title
   * scale_font: scaling factor applied to the default font size
   * disp_values: switches on/off the grading values on the plot
   * disp_rankings: switches on/off the rankings on the plot
   * rank_order: displays rankings in increasing (1) or decreasing (-1) order

Variables
---------

* clt (atmos, monthly mean, longitude latitude time)
* hus (atmos, monthly mean, longitude latitude lev time)
* od550aer, od870aer, od550abs, od550lt1aer (aero, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs, rsut, rsutcs (atmos, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude lev time)
* tas (atmos, monthly mean, longitude latitude time)
* toz (atmos, monthly mean, longitude latitude time)
* ts (atmos, monthly mean, longitude latitude time)
* ua (atmos, monthly mean, longitude latitude lev time)
* va (atmos, monthly mean, longitude latitude lev time)
* zg (atmos, monthly mean, longitude latitude lev time)

Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for non-obs4mips data for download instructions.*

* AIRS (hus – obs4mips)
* CERES-EBAF (rlut, rlutcs, rsut, rsutcs – obs4mips)
* ERA-Interim (tas, ta, ua, va, zg, hus – esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
* ESACCI-AEROSOL (od550aer, od870aer, od550abs, od550lt1aer – esmvaltool/utils/cmorizers/obs/cmorize_obs_ESACCI-AEROSOL.ncl)
* ESACCI-CLOUD (clt – esmvaltool/utils/cmorizers/obs/cmorize_obs_ESACCI-CLOUD.ncl)
* ESACCI-OZONE (toz – esmvaltool/utils/cmorizers/obs/cmorize_obs_ESACCI-OZONE.ncl)
* ESACCI-SST (ts – esmvaltool/utils/cmorizers/obs/cmorize_obs_ESACCI-SST.ncl)
* GPCP-SG (pr – obs4mips)
* HadISST (ts - esmvaltool/utils/cmorizers/obs/cmorize_obs_HadISST.ncl)
* MODIS (od550aer – obs4mips)
* NCEP (tas, ta, ua, va, zg – esmvaltool/utils/cmorizers/obs/cmorize_obs_NCEP.ncl)
* NIWA (toz – esmvaltool/utils/cmorizers/obs/cmorize_obs_NIWA.ncl)
* PATMOS (clt - esmvaltool/utils/cmorizers/obs/cmorize_obs_PATMOS.ncl)

References
----------

* Gleckler, P. J., K. E. Taylor, and C. Doutriaux, Performance metrics for climate models, J. Geophys. Res., 113, D06104, doi: 10.1029/2007JD008972 (2008).

* Righi, M., Eyring, V., Klinger, C., Frank, F., Gottschaldt, K.-D., Jöckel, P., and Cionni, I.: Quantitative evaluation of oone and selected climate parameters in a set of EMAC simulations, Geosci. Model Dev., 8, 733, doi: 10.5194/gmd-8-733-2015 (2015).

Example plots
-------------

.. centered:: |pic_permetrics1| |pic_permetrics2|

.. |pic_permetrics1| image:: /recipes/figures/perfmetrics/perfmetrics_fig_1.png
   :width: 50%

.. |pic_permetrics2| image:: /recipes/figures/perfmetrics/perfmetrics_fig_2.png
   :width: 30%

.. centered:: |pic_permetrics3| |pic_permetrics4|

.. |pic_permetrics3| image:: /recipes/figures/perfmetrics/perfmetrics_fig_3.png
   :width: 30%

.. |pic_permetrics4| image:: /recipes/figures/perfmetrics/perfmetrics_fig_4.png
   :width: 52%

.. figure:: /recipes/figures/perfmetrics/perfmetrics_fig_5.png
   :width: 75%
   :align: center

