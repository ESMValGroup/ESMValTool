.. _recipes_flato13ipcc:

IPCC AR5 Chapter 9 (selected figures)
=====================================

Overview
--------

The goal of this recipe is to collect diagnostics to reproduce Chapter 9 of AR5,
so that the plots can be readily reproduced and compared to previous CMIP
versions. In this way we can next time start with what was available in the
previous round and can focus on developing more innovative methods of analysis
rather than constantly having to "re-invent the wheel".

The plots are produced collecting the diagnostics from individual recipes. The
following figures from Flato et al. (2013) can currently be reproduced:

    * Figure 9.2 a,b,c: Annual-mean surface air temperature for the period
      1980-2005. a) multi-model mean, b) bias as the difference between the
      CMIP5 multi-model mean and the climatology from ERA-Interim
      (Dee et al., 2011), c) mean absolute model error with respect to the
      climatology from ERA-Interim.

    * Figure 9.4: Annual-mean precipitation rate (mm day-1) for the period
      1980-2005. a) multi-model mean, b) bias as the difference between the
      CMIP5 multi-model mean and the climatology from the Global Precipitation
      Climatology Project (Adler et al., 2003), c) difference between the
      multi-model mean and the ECMWF reanalysis of the seasonality, and d)
      difference between the multi-model mean and the ERA-Interim absolute
      seasonality.

    * Figure 9.5: Climatological (1985-2005) annual-mean cloud radiative
      effects in Wm-2 for the CMIP5 models against CERES EBAF (2001-2011) in
      Wm-2. Top row shows the shortwave effect; middle row the longwave effect,
      and bottom row the net effect. Multi-model-mean biases against CERES
      EBAF 2.6 are shown on the left, whereas the right panels show zonal
      averages from CERES EBAF 2.6 (black), the individual CMIP5 models (thin
      gray lines), and the multi-model mean (thick red line).

    * Figure 9.8: Observed and simulated time series of the anomalies in annual
      and global mean surface temperature. All anomalies are differences from
      the 1961-1990 time-mean of each individual time series. The reference
      period 1961-1990 is indicated by yellow shading; vertical dashed grey
      lines represent times of major volcanic eruptions. Single simulations
      for CMIP5 models (thin lines); multi-model mean (thick red line);
      different observations (thick black lines). Dataset pre-processing like
      described in Jones et al., 2013.

    * Figure 9.14: Sea surface temperature plots for zonal mean error, equatorial
      (5 deg north to 5 deg south) mean error, and multi model mean for zonal error
      and equatorial mean.

    * Figure 9.24: Time series of (a) Arctic and (b) Antarctic sea ice extent;
      trend distributions of (c) September Arctic and (d) February Antarctic
      sea ice extent.

    * Figure 9.42a: Equilibrium climate sensitivity (ECS) against the global
      mean surface air temperature, both for the period 1961-1990 and for the
      pre-industrial control runs.

    * Figure 9.42b: Transient climate response (TCR) against equilibrium climate
      sensitivity (ECS).

    * Figure 9.45a: Scatterplot of springtime snow-albedo effect values in climate
      change vs. springtime d(alpha\ :sub:`s`\)/d(T\ :sub:`s`\) values in the seasonal
      cycle in transient climate change experiments (Hall and Qu, 2006).

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_flato13ipcc.yml

Diagnostics are stored in esmvaltool/diag_scripts/

    * climate_metrics/ecs.py: See :ref:`recipes_ecs`.
    * clouds/clouds_bias.ncl: global maps of the multi-model mean and the multi-model
      mean bias (Fig. 9.2, 9.4)
    * clouds/clouds_isccp: global maps of multi-model mean minus observations + zonal
      averages of individual models, multi-model mean and observations (Fig. 9.5)
    * ipcc_ar5/tsline.ncl: time series of the global mean (anomaly) (Fig. 9.8)
    * ipcc_ar5/ch09_fig09_14.py: Zonally averaged and equatorial SST (Fig. 9.14)
    * seaice/seaice_tsline.ncl: Time series of sea ice extent (Fig. 9.24a/b)
    * seaice/seaice_trends.ncl: Trend distributions of sea ice extent (Fig 9.24c/d)
    * ipcc_ar5/ch09_fig09_42a.py: ECS vs. surface air temperature (Fig. 9.42a)
    * ipcc_ar5/ch09_fig09_42b.py: TCR vs. ECS (Fig. 9.42b)
    * emergent_constraints/snowalbedo.ncl: snow-albedo effect (Fig. 9.45a)

User settings in recipe
-----------------------

#. Script climate_metrics/ecs.py

   See :ref:`recipes_ecs`.

#. Script clouds/clouds_bias.ncl

#. Script clouds_bias.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * plot_abs_diff: additionally also plot absolute differences (true, false)
   * plot_rel_diff: additionally also plot relative differences (true, false)
   * projection: map projection, e.g., Mollweide, Mercator
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)

   * Required settings (variables)*

   * reference_dataset: name of reference datatset

   *Optional settings (variables)*

   * long_name: description of variable

   *Color tables*

   * variable "tas": diag_scripts/shared/plot/rgb/ipcc-tas.rgb,
     diag_scripts/shared/plot/rgb/ipcc-tas-delta.rgb
   * variable "pr-mmday": diag_scripts/shared/plots/rgb/ipcc-precip.rgb,
     diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script clouds/clouds_ipcc.ncl

   *Required settings (scripts)*

   none

   *Optional settings (scripts)*

   * explicit_cn_levels: contour levels
   * mask_ts_sea_ice: true = mask T < 272 K as sea ice (only for variable "ts");
     false = no additional grid cells masked for variable "ts"
   * projection: map projection, e.g., Mollweide, Mercator
   * styleset: style set for zonal mean plot ("CMIP5", "DEFAULT")
   * timemean: time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
     "annualclim" (annual mean)
   * valid_fraction: used for creating sea ice mask (mask_ts_sea_ice = true):
     fraction of valid time steps required to mask grid cell as valid data

   *Required settings (variables)*

   * reference_dataset:  name of reference data set

   *Optional settings (variables)*

   * long_name: description of variable
   * units: variable units

   *Color tables*

   * variables "pr", "pr-mmday": diag_scripts/shared/plot/rgb/ipcc-precip-delta.rgb

#. Script ipcc_ar5/tsline.ncl

   *Required settings for script*

   * styleset: as in diag_scripts/shared/plot/style.ncl functions

   *Optional settings for script*

   * time_avg: type of time average (currently only "yearly" and "monthly" are
     available).
   * ts_anomaly: calculates anomalies with respect to the defined period; for
     each gird point by removing the mean for the given calendar month
     (requiring at least 50% of the data to be non-missing)
   * ref_start: start year of reference period for anomalies
   * ref_end: end year of reference period for anomalies
   * ref_value: if true, right panel with mean values is attached
   * ref_mask: if true, model fields will be masked by reference fields
   * region: name of domain
   * plot_units: variable unit for plotting
   * y-min: set min of y-axis
   * y-max: set max of y-axis
   * mean_nh_sh: if true, calculate first NH and SH mean
   * volcanoes: if true, lines of main volcanic eruptions will be added
   * run_ave: if not equal 0 than calculate running mean over this number of
     years
   * header: if true, region name as header

   *Required settings for variables*

   none

   *Optional settings for variables*

   * reference_dataset: reference dataset; REQUIRED when calculating
     anomalies

   *Color tables*

   * e.g. diag_scripts/shared/plot/styles/cmip5.style

#. Script seaice/seaice_trends.ncl

   *Required settings (scripts)*

   * month: selected month (1, 2, ..., 12) or annual mean ("A")
   * region: region to be analyzed ( "Arctic" or "Antarctic")

   *Optional settings (scripts)*

   * fill_pole_hole: fill observational hole at North pole, Default: False

   *Optional settings (variables)*

   * ref_model: array of references plotted as vertical lines

#. Script seaice/seaice_tsline.ncl

   *Required settings (scripts)*

   * region: Arctic, Antarctic
   * month: annual mean (A), or month number (3 = March, for Antarctic; 9 = September for Arctic)

   *Optional settings (scripts)*

   * styleset: for plot_type cycle only (cmip5, cmip6, default)
   * multi_model_mean: plot multi-model mean and standard deviation (default: False)
   * EMs_in_lg: create a legend label for individual ensemble members (default: False)
   * fill_pole_hole: fill polar hole (typically in satellite data) with sic = 1 (default: False)

#. Script ipcc_ar5/ch09_fig09_42a.py

   *Required settings for script*

   none

   *Optional settings for script*

   * axes_functions: :obj:`dict` containing methods executed for the plot's
     :class:`matplotlib.axes.Axes` object.
   * dataset_style: name of the style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
   * matplotlib_style: name of the matplotlib style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python.matplotlib`).
   * save: :obj:`dict` containing keyword arguments for the function
     :func:`matplotlib.pyplot.savefig`.
   * seaborn_settings: Options for seaborn's ``set()`` method (affects all
     plots), see https://seaborn.pydata.org/generated/seaborn.set.html.

#. Script ipcc_ar5/ch09_fig09_42b.py

   *Required settings for script*

   none

   *Optional settings for script*

   * dataset_style: name of the style file (located in
     :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
   * log_x: Apply logarithm to X axis (ECS).
   * log_y: Apply logarithm to Y axis (TCR).
   * seaborn_settings: Options for seaborn's ``set()`` method (affects all
     plots), see https://seaborn.pydata.org/generated/seaborn.set.html.

#. Script emergent_constraints/snowalbedo.ncl

   *Required settings for script*

   * exp_presentday: name of present-day experiment (e.g. "historical")
   * exp_future: name of climate change experiment (e.g. "rcp45")

   *Optional settings for script*

   * diagminmax: observational uncertainty (min and max)
   * legend_outside: create extra file with legend (true, false)
   * styleset: e.g. "CMIP5" (if not set, this diagnostic will create its own
     color table and symbols for plotting)
   * suffix: string to be added to output filenames
   * xmax: upper limit of x-axis (default = automatic)
   * xmin: lower limit of x-axis (default = automatic)
   * ymax: upper limit of y-axis (default = automatic)
   * ymin: lower limit of y-axis (default = automatic)

   *Required settings for variables*

   * ref_model: name of reference data set

   *Optional settings for variables*

   none

Variables
---------

* areacello (fx, longitude latitude)
* pr (atmos, monthly mean, longitude latitude time)
* rlut, rlutcs (atmos, monthly mean, longitude latitude time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsuscs, rsdscs (atmos, monthly mean, longitude latitude time)
* rsut, rsutcs (atmos, monthly mean, longitude latitude time)
* sic (ocean-ice, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time)
* tos (ocean, monthly mean, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4mips data for download
instructions.*

* CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4mips)
* ERA-Interim (tas, ta, ua, va, zg, hus - esmvaltool/utils/cmorizers/obs/cmorize_obs_ERA-Interim.ncl)
* GPCP-SG (pr - obs4mips)
* HadCRUT4 (tas - esmvaltool/utils/cmorizers/obs/cmorize_obs_hadcrut4.ncl)
* HadISST (sic, tos - esmvaltool/utils/cmorizers/obs/cmorize_obs_hadisst.ncl)
* ISCCP-FH (rsuscs, rsdscs, rsdt - esmvaltool/utils/cmorizers/obs/cmorize_obs_isccp_fh.ncl)


References
----------

* Flato, G., J. Marotzke, B. Abiodun, P. Braconnot, S.C. Chou, W. Collins, P.
  Cox, F. Driouech, S. Emori, V. Eyring, C. Forest, P. Gleckler, E. Guilyardi,
  C. Jakob, V. Kattsov, C. Reason and M. Rummukainen, 2013: Evaluation of
  Climate Models. In: Climate Change 2013: The Physical Science Basis.
  Contribution of Working Group I to the Fifth Assessment Report of the
  Intergovernmental Panel on Climate Change [Stocker, T.F., D. Qin, G.-K.
  Plattner, M. Tignor, S.K. Allen, J. Boschung, A. Nauels, Y. Xia, V. Bex and
  P.M. Midgley (eds.)]. Cambridge University Press, Cambridge, United Kingdom
  and New York, NY, USA.

* Hall, A., and X. Qu, 2006: Using the current seasonal cycle to constrain
  snow albedo feedback in future climate change, Geophys. Res. Lett., 33,
  L03502, doi:10.1029/2005GL025127.

* Jones et al., 2013: Attribution of observed historical near-surface temperature
  variations to anthropogenic and natural causes using CMIP5 simulations. Journal
  of Geophysical Research: Atmosphere, 118, 4001-4024, doi:10.1002/jgrd.50239.


Example plots
-------------

.. _fig_flato13ipcc_1:
.. figure::  /recipes/figures/flato13ipcc/fig-9-2.png
   :align:   center

   Figure 9.2 a,b,c: Annual-mean surface air temperature for the period
   1980-2005. a) multi-model mean, b) bias as the difference between the
   CMIP5 multi-model mean and the climatology from ERA-Interim
   (Dee et al., 2011), c) mean absolute model error with respect to the
   climatology from ERA-Interim.

.. _fig_flato13ipcc_2:
.. figure::  /recipes/figures/flato13ipcc/fig-9-4.png
   :align:   center

   Figure 9.4: Annual-mean precipitation rate (mm day-1) for the period
   1980-2005. a) multi-model mean, b) bias as the difference between the
   CMIP5 multi-model mean and the climatology from the Global Precipitation
   Climatology Project (Adler et al., 2003), c) difference between the
   multi-model mean and the ECMWF reanalysis of the seasonality, and d)
   difference between the multi-model mean and the ERA-Interim absolute
   seasonality.

.. _fig_flato13ipcc_3:
.. figure::  /recipes/figures/flato13ipcc/fig-9-5.png
   :align:   center

   Figure 9.5: Climatological (1985-2005) annual-mean cloud radiative
   effects in Wm-2 for the CMIP5 models against CERES EBAF (2001-2011) in
   Wm-2. Top row shows the shortwave effect; middle row the longwave effect,
   and bottom row the net effect. Multi-model-mean biases against CERES
   EBAF 2.6 are shown on the left, whereas the right panels show zonal
   averages from CERES EBAF 2.6 (black), the individual CMIP5 models (thin
   gray lines), and the multi-model mean (thick red line).

.. _fig_flato13ipcc_4:
.. figure::  /recipes/figures/flato13ipcc/fig-9-8.png
   :align:   center

   Figure 9.8: Observed and simulated time series of the anomalies in annual
   and global mean surface temperature. All anomalies are differences from
   the 1961-1990 time-mean of each individual time series. The reference
   period 1961-1990 is indicated by yellow shading; vertical dashed grey
   lines represent times of major volcanic eruptions. Single simulations
   for CMIP5 models (thin lines); multi-model mean (thick red line);
   different observations (thick black lines). Dataset pre-processing like
   described in Jones et al., 2013.

.. _fig_flato13ipcc_5:
.. figure:: /recipes/figures/flato13ipcc/fig-9-14.png
   :align: center

   Figure 9.14: (a) Zonally averaged sea surface temperature (SST) error
   in CMIP5 models. (b) Equatorial SST error in CMIP5 models. (c) Zonally
   averaged multi-model mean SST error for CMIP5 together with
   inter-model standard deviation (shading). (d) Equatorial multi-model
   mean SST in CMIP5 together with inter-model standard deviation
   (shading) and observations (black).  Model climatologies are derived
   from the 1979-1999 mean of the historical simulations. The Hadley
   Centre Sea Ice and Sea Surface Temperature (HadISST) (Rayner et
   al., 2003) observational climatology for 1979-1999 is used as a
   reference for the error calculation (a), (b), and (c); and for
   observations in (d).

.. figure::  /recipes/figures/seaice/trend_sic_extend_Arctic_September_histogram.png
   :align:   center
   :width:   9cm

   Figure 9.24c: Sea ice extent trend distribution for the Arctic in September.

.. figure::  /recipes/figures/seaice/extent_sic_Arctic_September_1960-2005.png
   :align:   center
   :width:   12cm

   Figure 9.24a: Time series of total sea ice area and extent (accumulated) for the Arctic
   in September including multi-model mean and standard deviation.

.. _fig_flato13ipcc_6:
.. figure:: /recipes/figures/flato13ipcc/figure-9-42a.png
   :align: center

   Figure 9.42a: Equilibrium climate sensitivity (ECS) against the global mean
   surface air temperature of CMIP5 models, both for the period 1961-1990
   (larger symbols) and for the pre-industrial control runs (smaller symbols).

.. _fig_flato13ipcc_7:
.. figure:: /recipes/figures/flato13ipcc/fig-9-42b.png
   :align: center

   Figure 9.42b: Transient climate response (TCR) against equilibrium climate
   sensitivity (ECS) for CMIP5 models.

.. figure:: /recipes/figures/flato13ipcc/fig-9-45a.png
   :align: center

   Figure 9.45a: Scatterplot of springtime snow-albedo effect values in climate
   change vs. springtime :math:`\Delta \alpha_s`/:math:`\Delta T_s` values in
   the seasonal cycle in transient climate change experiments (CMIP5 historical
   experiments: 1901-2000, RCP4.5 experiments: 2101-2200).

