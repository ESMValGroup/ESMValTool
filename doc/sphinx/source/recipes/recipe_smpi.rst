Single Model Perfomance Index (SMPI)
====================================

Overview
--------

This diagnostic calculates the Single Model Performance Index (SMPI) following Reichler and Kim (2008). The SMPI (called "I\ :sup:`2`") is based on the comparison of several different climate variables (atmospheric, surface and oceanic) between climate model simulations and observations or reanalyses, and it focuses on the validation of the time-mean state of climate. For I\ :sup:`2` to be determined, the differences between the climatological mean of each model variable and observations at each of the available data grid points are calculated, and scaled to the interannual variance from the validating observations. This interannual variability is determined by performing a bootstrapping method (random selection with replacement) for the creation of a large synthetic ensemble of observational climatologies. The results are then scaled to the average error from a reference ensemble of models, and in a final step the mean over all climate variables and one model is calculated. The plot shows the I\ :sup:`2` values for each model (orange circles) and the multi-model mean (black circle), with the diameter of each circle representing the range of I\ :sup:`2` values encompassed by the 5th and 95th percentiles of the bootstrap ensemble. The I\ :sup:`2` values vary around one, with values greater than one for underperforming models, and values less than one for more accurate models. 

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipe/

* recipe_reichlerkim08bams.yml

Diagnostics are stored in diag_scripts/perfmetrics/

* main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields. The calculated fields can also be plotted as difference w.r.t. a given reference dataset. main.ncl also calculates RMSD, bias and taylor metrics. Input data have to be regridded to a common grid in the preprocessor. Each plot type is created by a separated routine, as detailed below.
* cycle_zonal.ncl: calculates single model perfomance index (Reichler and Kim, 2008). It requires fields precalculated by main.ncl.
* collect.ncl: collects the metrics previously calculated by cycle_latlon.ncl and passes them to the plotting functions.

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
   
   *Required settings for variables*
   
   * reference_dataset: reference dataset to compare with (usually the observations).
   
   *Optional settings for variables*

   * alternative_dataset: a second dataset to compare with.

These settings are passed to the other scripts by main.ncl, depending on the selected plot_type.

#. collect.ncl

   *Required settings for script*

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

* hfds (ocean, monthly mean, longitude latitude time)
* hus (atmos, monthly mean, longitude latitude lev time)
* pr (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* sic (ocean-ice, monthly mean, longitude latitude time) - not implemented yet
* ta (atmos, monthly mean, longitude latitude lev time)
* tas (atmos, monthly mean, longitude latitude time)
* tauu (atmos, monthly mean, longitude latitude time)
* tauv (atmos, monthly mean, longitude latitude time)
* tos (ocean, monthly mean, longitude latitude time)
* ua (atmos, monthly mean, longitude latitude lev time)
* va (atmos, monthly mean, longitude latitude lev time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.*

* ERA-Interim (hfds, hus, psl, ta, tas, tauu, tauv, ua, va -- obs4mips)
* HadISST (sic, tos -- reformat_scripts/obs/reformat_obs_ERA-Interim.ncl)
* GPCP monthly (pr -- reformat_scripts/obs/reformat_obs_HadISST.ncl)

References
----------

* Reichler, T. and J. Kim, How well do coupled models simulate today's climate? Bull. Amer. Meteor. Soc., 89, 303-311, doi: 10.1175/BAMS-89-3-303, 2008.

Example plots
-------------

.. figure:: /recipes/figures/smpi/reichlerkim08bams_smpi.png
   :width: 70 %
   
   Performance index I\ :sup:`2` for individual models (circles). Circle sizes indicate the length of the 95% confidence intervals. The black circle indicates the I\ :sup:`2` of the multi-model mean (similar to Reichler and Kim (2008), Figure 1).
