.. _nml_perfmetrics:

Performance metrics for essential climate parameters
====================================================

Overview
--------

The goal is to create a standard recipe for the calculation of performance metrics to quantify the ability of the models to reproduce the climatological mean annual cycle for selected "Essential Climate Variables" (ECVs) plus some additional corresponding diagnostics and plots to better understand and interpret the results. The recipe can be used to calculate performance metrics at different vertical levels (e.g., 5, 30, 200, 850 hPa as in Gleckler et al., 2008) and in four regions (global, tropics 20°N-20°S, northern extratropics 20°-90°N, southern extratropics 20°-90°S). As an additional reference, we consider the Righi et al. (2015) paper.

Available recipes and diagnostics
-----------------------------------

Recipes are stored in nml/

* recipe_perfmetrics_CMIP5.xml

Diagnostics are stored in diag_scripts/

* perfmetrics_grading.ncl: calculates grades according to a given metric, with different options for normalization. It requires fields precalculated by perfmetrics_main.ncl.
* perfmetrics_grading_collect.ncl: collects results from metrics previously calculated by perfmetrics_grading.ncl and passes them to the plotting functions.
* perfmetrics_main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields from input monthly 2-d or 3-d ("T2M", "T3Ms") data. The calculated fields can be also plotted as difference w.r.t. a given reference model. They are also used as input to calculate grading metrics (see perfmetrics_grading.ncl).
* perfmetrics_taylor.ncl: calculates grades according to a given metric, with different options for normalization. It requires fields precalculated by perfmetrics_main.ncl.
* perfmetrics_taylor_collect.ncl: collects results from metrics previously calculated by perfmetrics_taylor.ncl and passes them to the plotting functions.

User settings
-------------

User setting files (cfg files) are stored in nml/cfg_perfmetrics/CMIP5/

#. perfmetrics_grading.ncl

   *diag_script_info attributes*

   * MultiModelMean: calculate multi-model mean (True, False)
   * MultiModelMedian: calculate multi-model median (True, False)
   * metric: applied metric ("RMSD" = root-mean square difference, "BIAS" = mean bias, "stddev_ratio" = ratio of standard deviations of var and ref (for Taylor diagrams only), "correlation" = pattern correlation of var and ref (for Taylor diagrams only)).
   * normalization: applied normalization ("mean" = normalization with mean, "median" = normalization with media, "centered_median" = substracting and dividing by the median, "stddev_mean" = normalization with substracting the mean and dividing by the standard deviation)

#. perfmetrics_grading_collect.ncl

   *Required diag_script_info attributes*

   * label_bounds: min and max of the labelbar
   * label_scale: bin width of the labelbar
   * disp_values: switch on/off the grading values on the plot

   *Optional diag_script_info attributes*

   * sort: sort models in alphabetic order (excluding multi-model mean)
   * title: plot title
   * scale_font: scaling factor applied to the default font size

#. perfmetrics_main.ncl

   *diag_script_info attributes*

   * plot_type: plot type ("cycle" (time), "zonal" (plev, lat), "latlon" (lat, lon), "cycle_latlon" (time, lat, lon))
   * time_avg: time averaging ("monthlyclim", "seasonalclim")
   * valid_fraction: required fraction of valid values
   * level: vertical level (hPa, "all" for no selection; set to "all" for zonal mean plots)
   * region: averaging region ("Global", "Tropics", "NH extratropics", "SH extratropics")
   * grid: regridding option ("finest", "coarsest", "ref")
   * draw_plots: create plots (True, False)
   * plot_diff: create difference plots (only for zonal and lat-lon plots) (True, False)
   * plot_stddev: plot standard deviation ("all", "none", "ref_model" or given model name)
   * legend_outside: plot legend in a separate file (only for cycle plots) (True, False)
   * styleset: plot style (only for cycle plots) ("CMIP5", "DEFAULT", "EMAC")
   * t_test: calculate t-test for difference plots (only for zonal and lat-lon plots) (True, False)
   * conf_level: confidence level for the t-test (only for zonal and lat-lon plots)

#. perfmetrics_taylor.ncl

   *Required diag_script_info attributes*

   * region: averaging region ("Global", "Tropics", "NH extratropics", "SH extratropics")
   * time_avg: time averaging ("monthlyclim", "seasonalclim")
   * metric: selected metric (required but ignored by permetrics_taylor.ncl)
   * normalization: type of metric normalization (required but ignored by permetrics_taylor.ncl)

#. perfmetrics_taylor_collect.ncl

   *diag_script_info attributes*

   * None.

Variables
---------

* hus (atmos, monthly mean, longitude latitude lev time)
* od550aer (aero, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs, rsut, rsutcs (atmos, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude lev time)
* tas (atmos, monthly mean, longitude latitude time)
* ua (atmos, monthly mean, longitude latitude lev time)
* va (atmos, monthly mean, longitude latitude lev time)
* zg (atmos, monthly mean, longitude latitude lev time)

Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.*

* AIRS L3 (hus – obs4mips)
* CERES-EBAF (rlut, rlutcs, rsut, rsutcs – obs4mips)
* ERA-Interim (tas, ta, ua, va, zg, hus – reformat_scripts/obs/reformat_obs_ERA-Interim.ncl)
* ESACCI-AEROSOL (od550aer – reformat_scripts/obs/reformat_obs_ESACCI-AEROSOL.ncl)
* GPCP-SG (pr – obs4mips)
* MODIS-L3 (od550aer – obs4mips)
* NCEP (tas, ta, ua, va, zg – reformat_scripts/obs/reformat_obs_NCEP.ncl)

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

