.. _nml_perfmetrics:

IPCC Chapter 12 selected figures
====================================================

Overview
--------

The goal is to create a standard recipe for creating selected Figures from IPCC AR5 Chapter 12 on "Long-term Climate Change: Projections, Commitments and Irreversibility". These include maps showing the change in a variable between historical and future periods, zonal means for 3D variables, timeseries showing the change in certain variables from historical to future periods for multiple scenarios, and maps visualizing change in variables normalized by global mean temperature change (pattern scaling) as in Collins et al., 2013.


Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_collins13ipcc.xml

Diagnostics are stored in diag_scripts/ipcc_ar5/

* ch12_map_diff_each_model_fig12-9.ncl: 
* ch12_ts_line_mean_spread.ncl: 
* ch12_plot_ts_line_mean_spread.ncl: 
* ch12_calc_IAV_for_stippandhatch.ncl: 
* ch12_calc_map_diff_mmm_stippandhatch.ncl: 
* ch12_plot_map_diff_mmm_stipp.ncl:
* ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl: 
* ch12_plot_zonal_diff_mmm_stipp.ncl:
* ch12_calc_map_diff_scaleT_mmm_stipp.ncl:
* ch12_snw_area_change_fig12-32.ncl:

User settings
-------------

#. ch12_map_diff_each_model_fig12-9.ncl

   *Required attributes*

   * time_avg: time averaging ("annualclim", "seasonalclim")
   * experiment: IPCC Scenario, used to pair historical and rcp runs from same model

   *Optional attributes*

   * projection: map projection, any valid ncl projection, default = Robinson
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * title: plot title
   * colormap: alternative colormap, path to rgb file or ncl name
   * diff_levs: list with contour levels for plots
   * span: span whole colormap? (True, False, default = False)

#. ch12_ts_line_mean_spread.ncl

   *Required attributes*

   * scenarios: list with scenarios included in figure
   * syears: list with start years in time periods (e.g. start of historical period and rcps)
   * eyears: list with end years in time periods (end year of historical runs and rcps)
   * begin_ref_year: start year of reference period (e.g. 1986)
   * end_ref_year: end year of reference period (e.g 2005)
   * label: list with labels to use in legend depending on scenarios

   *Optional attributes*

   * spread: how many standard deviations to calculate the spread with? default is 1., ipcc tas used 1.64
   * model_nr: save number of model runs per period and scenario in netcdf to print in plot? (True, False, default = False)
   * ts_minlat: minimum latitude if not global
   * ts_maxlat: maximum latitude if not global
   * ts_minlon: minimum longitude if not global
   * ts_maxlon: maximum longitude if not global

#. ch12_plot_ts_line_mean_spread.ncl: 

   *Required attributes*

   * ancestors: variable and diagnostics that calculated data to be plotted

   *Optional attributes*

   * title: specify plot title
   * yaxis: specify y-axis title
   * ymin: minimim value on y-axis, default calculated from data
   * ymax: maximum value on y-axis
   * colormap: alternative colormap, path to rgb file or ncl name

#. ch12_calc_IAV_for_stippandhatch.ncl:

   *Required attributes*

   * time_avg: time averaging ("annualclim", "seasonalclim"), needs to be consistent with calculation in ch12_calc_map_diff_mmm_stippandhatch.ncl

   *Optional attributes*

   * periodlength: length of period in years to calculate variability over, default is total time period
   * iavmode: calculate IAV from multi-model mean or save individual models ("each": save individual models, "mmm": multi-model mean, default), needs to be consistent with ch12_calc_map_diff_mmm_stippandhatch.ncl

#. ch12_calc_map_diff_mmm_stippandhatch.ncl:

   *Required attributes*

   * ancestors: variable and diagnostics that calculated interannual variability for stippling and hatching
   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional attributes*

   * seasons: list with seasons index if time_avg "seasonalclim" (then required),  DJF:0, MAM:1, JJA:2, SON:3
   * iavmode: calculate IAV from multi-model mean or save individual models ("each": save individual models, "mmm": multi-model mean, default), needs to be consistent with ch12_calc_IAV_for_stippandhatch.ncl
   * percent: determines if difference expressed in percent (0, 1, default = 0)

#. ch12_plot_map_diff_mmm_stipp.ncl:

   *Required attributes*

   * ancestors: variable and diagnostics that calculated field to be plotted

   *Optional attributes*

   * projection: map projection, any valid ncl projection, default = Robinson
   * diff_levs: list with explicit levels for all contour plots
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * model_nr: save number of model runs per period and scenario in netcdf to print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * span: span whole colormap? (True, False, default = True)
   * sig: plot stippling for significance? (True, False)
   * not_sig: plot hatching for uncertainty? (True, False)
   * pltname: alternative name for output plot, default is diagnostic + varname + time_avg
   * plotmask: apply a mask when plotting? ("None" (default), "ocean", "land")
   * units: units written next to colorbar, e.g (~F35~J~F~C)

#. ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl:

   *Required attributes*

   * ancestors: variable and diagnostics that calculated interannual variability for stippling and hatching
   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional attributes*

   * base_cn: if want contours of base period as contour lines, need to save base period field (True, False)
   * seasons: list with seasons index if time_avg "seasonalclim" (then required),  DJF:0, MAM:1, JJA:2, SON:3
   * iavmode: calculate IAV from multi-model mean or save individual models ("each": save individual models, "mmm": multi-model mean, default), needs to be consistent with ch12_calc_IAV_for_stippandhatch.ncl
   * percent: determines if difference expressed in percent (0, 1, default = 0)

#. ch12_plot_zonal_diff_mmm_stipp.ncl:

   *Required attributes*

   * ancestors: variable and diagnostics that calculated field to be plotted

   *Optional attributes*

   * diff_levs: list with explicit levels for all contour plots
   * max_vert: maximum number of plots in vertical
   * max_hori: maximum number of plots in horizontal
   * model_nr: save number of model runs per period and scenario in netcdf to print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * span: span whole colormap? (True, False, default = True)
   * sig: plot stippling for significance? (True, False)
   * not_sig: plot hatching for uncertainty? (True, False)
   * pltname: alternative name for output plot, default is diagnostic + varname + time_avg
   * units: units written next to colorbar in ncl strings, e.g (m s~S~-1~N~)
   * if base_cn: True in ch12_calc_zonal_cont_diff_mmm_stippandhatch.ncl further settings to control contour lines
        base_cnLevelSpacing: spacing between contour levels
        base_cnMinLevel: minimum contour line
        base_cnMaxLevel: maximum contour line

#. ch12_calc_map_diff_scaleT_mmm_stipp.ncl:

   *Required attributes*

   * time_avg: time averaging ("annualclim", "seasonalclim")
   * scenarios: list with scenarios to be included
   * periods: list with start years of periods to be included
   * label: list with labels to use in legend depending on scenarios

   *Optional attributes*

   * seasons: list with seasons index if time_avg "seasonalclim" (then required),  DJF:0, MAM:1, JJA:2, SON:3
   * percent: determines if difference expressed in percent (0, 1, default = 0)

#. ch12_snw_area_change_fig12-32.ncl:

   *Required attributes*

   * scenarios: list with scenarios included in figure
   * syears: list with start years in time periods (e.g. start of historical period and rcps)
   * eyears: list with end years in time periods (end year of historical runs and rcps)
   * begin_ref_year: start year of reference period (e.g. 1986)
   * end_ref_year: end year of reference period (e.g 2005)
   * months: first letters of  months included in analysis? e.g. for MA (March + April) for Northern Hemisphere
   * label: list with labels to use in legend depending on scenarios

   *Optional attributes*

   * spread: how many standard deviations to calculate the spread with? default is 1., ipcc tas used 1.64
   * model_nr: save number of model runs per period and scenario in netcdf to print in plot? (True, False, default = False)
   * colormap: alternative colormap, path to rgb file or ncl name
   * ts_minlat: minimum latitude if not global
   * ts_maxlat: maximum latitude if not global
   * ts_minlon: minimum longitude if not global
   * ts_maxlon: maximum longitude if not global

Variables
---------

*Note: These are the variables tested and used in IPCC AR5. However, the code is flexible and in theory other variables of the same kind can be used.*

* tas (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rsut, rtmt (atmos, monthly mean, longitude latitude time)
* hurs (atmos, monthly mean, longitude latitude time)
* clt (atmos, monthly mean, longitude latitude time)
* psl (atmos, monthly mean, longitude latitude time)
* evspsbl (atmos, monthly mean, longitude latitude time)
* mrsos (land, monthly mean, longitude latitude time)
* mrro (land, monthly mean, longitude latitude time)
* sos (ocean, monthly mean, longitude latitude time)
* ta (atmos, monthly mean, longitude latitude lev time)
* ua (atmos, monthly mean, longitude latitude lev time)
* thetao (ocean, monthly mean, longitude latitude lev time)
* snw (land, monthly mean, longitude latitude time)

Observations and reformat scripts
---------------------------------

*Note: No observations are used since the comparison is between historical and scenario runs.*

References
----------

* Collins, M., R. Knutti, J. Arblaster, J.-L. Dufresne, T. Fichefet, P. Friedlingstein, X. Gao, W.J. Gutowski, T. Johns, G. Krinner, M. Shongwe, C. Tebaldi, A.J. Weaver and M. Wehner, 2013: Long-term Climate Change: Projections, Commitments and Irreversibility. In: Climate Change 2013: The Physical Science Basis. Contribution of Working Group I to the Fifth Assessment Report of the Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K. Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom and New York, NY, USA.


Example plots
-------------

.. centered:: |pic_collins1| |pic_collins2|

.. |pic_collins1| image:: /recipes/figures/collins13ipcc/collins13ipcc_fig_1.png
   :width: 50%

.. |pic_collins2| image:: /recipes/figures/collins13ipcc/collins13ipcc_fig_2.png
   :width: 30%

.. centered:: |pic_collins3| |pic_collins4|

.. |pic_collins3| image:: /recipes/figures/collins13ipcc/collins13ipcc_fig_3.png
   :width: 30%

.. |pic_collins4| image:: /recipes/figures/collins13ipcc/collins13ipcc_fig_4.png
   :width: 52%

.. figure:: /recipes/figures/collins13ipcc/collins13ipcc_fig_5.png
   :width: 75%
   :align: center

