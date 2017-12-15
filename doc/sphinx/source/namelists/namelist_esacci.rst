.. _nml_esacci:


ESA CCI
=======

Overview
--------

Namelist for creating the 20 figures from Lauer et al. (2017) using the European Space Agency’s Climate Change Initiative (ESA CCI) data sets for sea surface temperature, sea ice, cloud, soil moisture, land cover, aerosol, ozone, and greenhouse gases (CO\ :sub:`2`\).
This namelist demonstrates the value of the ESA CCI data for model evaluation providing an overview on the possible applications of the new data to evaluating CMIP models.


Available namelists and diagnostics
-----------------------------------

Namelists are stored in nml/

   * namelist_lauer17rse.xml

Diagnostics are stored in diag_scripts/

   * aerosol_stations.ncl: comparison of ESA CCI aerosol with AeroNet and MODIS
   * clouds.ncl: global maps of (multi-year) annual means including multi-model mean
   * clouds_interannual.ncl: global maps of interannual variability of cloud properties
   * clouds_ipcc.ncl: maps of multi-model mean bias and zonal averages
   * clouds_taylor.ncl: taylor diagrams
   * eyring13jgr_fig01.ncl: calculates seasonal cycles of zonally averaged total ozone columns.
   * eyring13jgr_fig02.ncl: time series of area-weighted total ozone from 1960-2005 for the annual mean averaged over the global domain (90°S-90°N), Tropics (25°S-25°N), northern mid-latitudes (35°N-60°N), southern mid-latitudes (35°S-60°S), and the March and October mean averaged over the Arctic (60°N-90°N) and the Antarctic (60°S-90°S).
   * eyring13jgr_fig04.ncl: climatological annual mean tropospheric ozone columns (geographical distribution).
   * lc_ESACCI.py: ESA CCI land cover diagnostics including global maps of grass and cropland cover and of forest and shrub cover
   * perfmetrics_grading.ncl: calculates grades according to a given metric with different options for normalization. It requires fields precalculated by perfmetrics_main.ncl (see :ref:`nml_perfmetrics`).
   * perfmetrics_grading_collect.ncl: collects results from metrics previously calculated by perfmetrics_grading.ncl and passes them to the plotting functions (see :ref:`nml_perfmetrics`).
   * perfmetrics_main.ncl: calculates and (optionally) plots annual/seasonal cycles, zonal means, lat-lon fields and time-lat-lon fields from input monthly 2-d or 3-d ("T2M", "T3Ms") data. The calculated fields can be also plotted as difference w.r.t. a given reference model. They are also used as input to calculate grading metrics (see perfmetrics_grading.ncl) (see :ref:`nml_perfmetrics`).
   * SeaIce_polcon.ncl: polar stereographic plots of sea ice concentration (= sea ice area fraction) and extent (grid cells with a sea ice concentration of at least 15%) for individual models or observational data sets, for Arctic and Antarctic regions with flexible paneling of the individual plots. The edges of sea ice extent can be highlighted via an optional red line.
   * SeaIce_polcon_diff.ncl: polar stereographic plots of sea ice area concentration difference between individual models and reference data (e.g., an observational data set) for both Arctic and Antarctic with flexible paneling of the individual plots. All data are regridded to a common grid (1°x1°) before comparison.
   * SeaIce_tsline.ncl: time series line plots of total sea ice area and extent (accumulated) for northern and southern hemispheres with optional multi-model mean and standard deviation. One value is used per model per year, either annual mean or the mean value of a selected month.
   * sm_ESACCI.py: ESA CCI soil moisture diagnostics including global maps of temporal trend in soil moisture and percentile maps for soil moisture
   * sst_ESACCI.py: ESA CCI SST diagnostics including global maps of absolute and relative SST bias and time series of mean SST for different ocean basins
   * tsline.ncl: time line plots of annual means for spatial averages
   * vpline.ncl: produces vertical profiles according to Eyring et al. (2006) Figure 5 upper panels following eyring06jgr_fig05.ncl.


User settings
-------------

User setting files (cfg files) are stored in nml/cfg_lauer17rse/

   * cfg_aerosol_stations_AERONET.ncl
   * cfg_carbon_line_3030.ncl
   * cfg_carbon_line_3060.ncl
   * cfg_carbon_line_6030.ncl
   * cfg_carbon_line_h.ncl
   * cfg_clouds_err.ncl
   * cfg_clouds_interannual_esa.ncl
   * cfg_clouds_ipcc.ncl
   * cfg_clouds.ncl
   * cfg_clouds_taylor_esa.ncl
   * cfg_clouds_taylor_esa-sic.ncl
   * cfg_dummy.conf
   * cfg_esacci_vpline.ncl
   * cfg_eyring13jgr_fig01.ncl
   * cfg_eyring13jgr_fig01_NIWA.ncl
   * cfg_eyring13jgr_fig02.ncl
   * cfg_eyring13jgr_fig04.ncl
   * cfg_lc_ESACCI.py
   * cfg_perfmetrics_grading_collect.ncl
   * cfg_perfmetrics_grading_RMSD_200_glob.ncl
   * cfg_perfmetrics_grading_RMSD_400_glob.ncl
   * cfg_perfmetrics_grading_RMSD_500_glob.ncl
   * cfg_perfmetrics_grading_RMSD_850_glob.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob_aero.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob_sm.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob_toz.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob_ts.ncl
   * cfg_perfmetrics_grading_RMSD_all_glob_xco2.ncl
   * cfg_perfmetrics_grading_RMSD_all_NHpolar_sic.ncl
   * cfg_perfmetrics_grading_RMSD_all_SHpolar_sic.ncl
   * cfg_perfmetrics_grading_RMSD_all_SHpolar_toz.ncl
   * cfg_perfmetrics_latlon_annualclim_all_glob_aerosol.ncl
   * cfg_perfmetrics_latlon_annualclim_all_glob.ncl
   * cfg_SeaIce_NH.ncl
   * cfg_SeaIce_SH.ncl
   * cfg_sm_ESACCI.py
   * cfg_sst_ESACCI_fig3.py
   * cfg_sst_ESACCI_fig4.py


Variables
---------

* abs550aer
* clt, cltStderr
* grassNcropFrac
* hus
* LW_CRE
* od550aer, od550aerStderr
* od550lt1aer
* od870aer, od870aerStderr
* pr
* rlut, rsut
* shrubNtreeFrac
* sic
* sm, smStderr
* SW_CRE
* ta
* tas
* tos
* toz, tozStderr
* tro3prof
* ts, tsStderr
* ua, va
* xco2, xco2Stderr
* zg


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.*

   * AIRS (hus): obs4mips
   * BDBP (tro3prof): reformat_scripts/obs/reformat_obs_BDBP.ncl
   * CERES-EBAF (LW_CRE, rlut, rsut, SW_CRE): obs4mips
   * CLARA-A2 (clt): *contact ESMValtool development team*
   * ERA-Interim (hus, ta, tas, ua, va, zg): reformat_scripts/obs/reformat_obs_ERA-Interim.ncl, reformat_scripts/obs/reformat_obs_ERA-Interim_surffluxes.ncl
   * ESACCI-AEROSOL (abs550aer, od550aer, od550aerStderr, od550lt1aer, od870aer, od870aerStder): reformat_scripts/obs/reformat_obs_ESACCI-AEROSOL.ncl
   * ESACCI-CLOUD (clt, cltStderr): reformat_scripts/obs/reformat_obs_ESACCI-CLOUD.ncl
   * ESACCI-GHG (xco2, xco2Stderr): reformat_scripts/obs/reformat_obs_ESACCI-GHG.csh
   * ESACCI-LANDCOVER (grassNcropFrac, shrubNtreeFrac): reformat_scripts/obs/reformat_obs_ESACCI-LANDCOVER.py
   * ESACCI-OZONE (toz, tozStderr, tro3prof): reformat_scripts/obs/reformat_obs_ESACCI-OZONE.ncl, reformat_scripts/obs/reformat_obs_ESACCI-OZONE_LP.ncl
   * ESACCI-SIC (sic): reformat_scripts/obs/reformat_obs_ESACCI-sic.ncl
   * ESACCI-SOILMOISTURE (sm, smStderr): reformat_scripts/obs/reformat_obs_ESACCI-SOILMOISTURE.ncl
   * ESACCI-SST (tos, ts, tsStderr): reformat_scripts/obs/reformat_obs_ESACCI-SST.ncl
   * GPCP-SG (pr): obs4mips
   * HadISST (ts): reformat_scripts/obs/reformat_obs_HadISST.ncl
   * MODIS-L3-C6 (clt, od550aer): reformat_scripts/obs/reformat_obs_MODIS-L3-C6.ncl
   * NCEP (ta, tas, ua, va, zg): reformat_scripts/obs/reformat_obs_NCEP.ncl
   * NIWA (toz): reformat_scripts/obs/reformat_obs_NIWA.ncl
   * NSIDC-NT (sic): reformat_scripts/obs/reformat_obs_NSIDC.ncl
   * PATMOS (clt): *contact ESMValtool development team*


References
----------

* Lauer, A., V. Eyring, M. Righi, M. Buchwitz, P. Defourny, M. Evaldsson, P. Friedlingstein, R. de Jeuf, G. de Leeuw, A. Loew, C. J. Merchant, B. Müller, T. Popp, M. Reuter, S. Sandven, D. Senftleben, M. Stengel, M. Van Roozendael, S. Wenzel, and U. Willén: Benchmarking CMIP5 models with a subset of ESA CCI Phase 2 data using the ESMValTool, Remote Sensing of Environment, http://dx.doi.org/10.1016/j.rse.2017.01.007, 2017.


Example plots
-------------

.. _fig_esacci_1:
.. figure::  /namelists/figures/esacci/Lauer17_fig01.png
   :align:   center
   :width:   14cm

   Relative space-time root-mean-square deviation (RMSD) calculated from the climatological seasonal cycle of the CMIP5 simulations (Lauer et al. 2017, Fig. 1).

.. _fig_esacci_2:
.. figure::  /namelists/figures/esacci/Lauer17_fig02.png
   :align:   center

   Extended Taylor diagrams showing the multi-year annual average performance of CMIP5 models in comparison with ESA CCI data (Lauer et al. 2017, Fig. 2).

.. _fig_esacci_3:
.. figure::  /namelists/figures/esacci/Lauer17_fig03.png
   :align:   center
   :width:   14cm

   Temporal means of SST in K for the ESA CCI data set (top right) and the CMIP5 model MPI-ESM (top left) as well as absolute (bottom left) and relative differences (bottom right) (Lauer et al. 2017, Fig. 3).

.. _fig_esacci_4:
.. figure::  /namelists/figures/esacci/Lauer17_fig04.png
   :align:   center
   :width:   11cm

   Time series of SST for different ocean basins from 7 CMIP5 models compared with the ESA CCI SST data (Lauer et al. 2017, Fig. 4).

.. _fig_esacci_5:
.. figure::  /namelists/figures/esacci/Lauer17_fig05.png
   :align:   center
   :width:   10cm

   Evolution (1960-2020) of September Arctic sea ice extent in million km\ :sup:`2`\  from the CMIP5 models (colored lines) and from observations (thick black lines) (Lauer et al. 2017, Fig. 5).

.. _fig_esacci_6:
.. figure::  /namelists/figures/esacci/Lauer17_fig06.png
   :align:   center
   :width:   12cm

   Polar-stereographic map of Arctic September (upper row) and Antarctic March (lower row) sea ice concentration from ESA CCI SI SSM/I (left column) and NSIDC-NT (middle column) observations averaged over the years 1992-2008. The right column depicts the differences between the CMIP5 multi-model mean and the ESA CCI SI SSM/I observations averaged over the years 1992-2005 (Lauer et al. 2017, Fig. 6).

.. _fig_esacci_7:
.. figure::  /namelists/figures/esacci/Lauer17_fig07.png
   :align:   center

   Maps of the multi-years seasonal mean of total cloud cover, 1-sigma uncertainty from ESA CCI cloud, the differences between the ESA CCI data and the CMIP5 multi-model mean, and zonal means (Lauer et al. 2017, Fig. 7).

.. _fig_esacci_8:
.. figure::  /namelists/figures/esacci/Lauer17_fig08.png
   :align:   center

   Interannual variability in total cloud cover estimate from relative temporal standard deviation of the deseasonalized monthly means time series from 1982 to 2014 (Lauer et al. 2017, Fig. 8).

.. _fig_esacci_9:
.. figure::  /namelists/figures/esacci/Lauer17_fig09.png
   :align:   center

   Temporal mean fields of volumetric soil moisture from the CNRM-CM5 model (top left), the ESA CCI soil moisture data set (top right) as well as their absolute (bottom left) and relative differenecs (bottom right) (Lauer et al. 2017, Fig. 9).

.. _fig_esacci_10:
.. figure::  /namelists/figures/esacci/Lauer17_fig10.png
   :align:   center

   Temporal trend in soil moisture over the period 1988-2008 as derived from the CNRM-CM5 model (left) and the ESA CCI soil moisture data sets (right) (Lauer et al. 2017, Fig. 10).

.. _fig_esacci_11:
.. figure::  /namelists/figures/esacci/Lauer17_fig11.png
   :align:   center
   :width:   11cm

   Percentile maps for ESA CCI soil moisture (left column) and soil moisture from CNRM-CM5 (right column) (Lauer et al. 2017, Fig. 11).

.. _fig_esacci_12:
.. figure::  /namelists/figures/esacci/Lauer17_fig12.png
   :align:   center
   :width:   11cm

   Area fraction (%) of forest and shrub cover in the MPI-ESM-MR model (top left) and the ESA CCI land cover data set (top right) and absolute (bottom left) and relative differences (bottom right) (Lauer et al. 2017, Fig. 12).

.. _fig_esacci_13:
.. figure::  /namelists/figures/esacci/Lauer17_fig14.png
   :align:   center
   :width:   14cm

   Climatological mean AOD (left column), fine mode optical depth (middle) and absorption optical depth (right column) at 550 nm averaged over the period 1997-2011. The first row shows the the observations (ESA CCI ATSR SU v4.21), the other rows the differences between selected CMIP5 models with interactive aerosols and the ESA CCI data (Lauer et al. 2017, Fig. 14).

.. _fig_esacci_14:
.. figure::  /namelists/figures/esacci/Lauer17_fig15.png
   :align:   center
   :width:   12cm

   Comparison of AOD at 550 nm from the ESA CCI ATSR SU v4.21 and the MODIS Terra C6 satellite products against the AERONET ground-based measurements for the period 2003-2011. The top row shows the AERONET values as open circles plotted on top of the satellite data averaged over the same time period. The bottom row shows scatter plots of spatially and temporally collocated measurements on a monthly-mean basis (Lauer et al. 2017, Fig. 15).

.. _fig_esacci_15:
.. figure::  /namelists/figures/esacci/Lauer17_fig16.png
   :align:   center
   :width:   12cm

   Time series of area-weighted total column ozone from 1960 to 2010 for a) global annual mean (90°S-90°N) and b) Antarctic October mean (60°S-90°S). The figure shows the multi-model mean (black line) and standard deviation (gray shading) as well as individual CMIP5 models with interactive chemistry (colored lines) compared with ESA CCI (filled circles) and NIWA (open triangles) data (Lauer et al. 2017, Fig. 16).

.. _fig_esacci_16:
.. figure::  /namelists/figures/esacci/Lauer17_fig17.png
   :align:   center
   :width:   14cm

   Vertical ozone profile climatologies (2007-2008) at a) 80°N in March, b) the equator in March, and c) at 80°S in October from individual CMIP5 models with interactive chemistry (colored lines) and the ESA CCI ozone data set (solid black line). The multi-model mean (MMM) is shown as a red solid line with one standard deviation of the inter-model spread shown as the light-blue shaded area (Lauer et al. 2017, Fig. 17).

.. _fig_esacci_17:
.. figure::  /namelists/figures/esacci/Lauer17_fig18.png
   :align:   center
   :width:   12cm

   Total column ozone climatologies (1997-2010) for (upper row, from left to right) the multi-model mean of CMIP5 models with interactive chemistry, the ESA CCI ozone data set, and the differences between the CMIP5 multi-model mean and the ESA CCI ozone data. The lower row shows the same plots but for the NIWA combined total column ozone data (Lauer et al. 2017, Fig. 18).

.. raw:: latex

    \clearpage 

.. _fig_esacci_18:
.. figure::  /namelists/figures/esacci/Lauer17_fig19.png
   :align:   center
   :width:   13cm

   Time series of column averaged carbon dioxide (XCO\ :sub:`2`\) from 2003 to 2014 from the CMIP5 emission driven simulations for the historical period (2003 to 2005) extended with RCP8.5 simulations (from 2006 to 2014) in comparison with the ESA CCI GHG XCO\ :sub:`2` data (Lauer et al. 2017, Fig. 19).

.. _fig_esacci_19:
.. figure::  /namelists/figures/esacci/Lauer17_fig20.png
   :align:   center
   :width:   14cm

   Annual mean XCO\ :sub:`2` climatologies averaged over the years 2003-2008 (top row) and over the years 2009-2014 (bottom row). Shown are deviations from the global annual mean (printed in the right above each panel) for (left) the CMIP5 multi-model mean and (middle) ESA CCI XCO\ :sub:`2`\. The right panels show the absolute differences between the CMIP5 multi-model mean and ESA CCI XCO\ :sub:`2` data (Lauer et al. 2017, Fig. 20).

