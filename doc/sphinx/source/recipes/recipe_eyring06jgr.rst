.. _recipes_<mynewrecipe>:

Diagnostics of stratospheric dynamics and chemistry reproducing selected figures from Eyring et al. JGR (2006).

=====

Overview
--------

This recipe reproduces the figures of `Eyring et al. (2006)`_
The following plots are reproduced:
 * Vertical profile climatological mean bias of climatological mean for selected seasons and latitudinal region. 
 * Transition to easterlies at 60S, climatological mean for selected seasons and latitudinal region.
 * Scatter plot heat fluxes vs. temperature at selected latitude and month.
 * Time series of monthly anomalies respect to a reference period at selected level, trend climatological mean for selected seasons and latitudinal region.
 * Vertical and latitudinal profile of climatological mean for selected seasons this figure and setting is valid for figure 5 (CH4) figure 6 (H2O) figure 11(HCL) figure 13 (tro3).
 * Climatological mean at selected level (100hPa)  and latitudinal region (Equator) a) temperature and b)water vapor.
 * Time-height sections of water vapor mixing ratio shown as the deviation (in parts per million by volume) from the time mean profile, averaged between 10°S and 10°.
 * Tape recorder phase and amplitude climatological mean for selected latitudinal region at all levels.
 * Latitudinal profiles of mean age of air climatological mean for selected seasons and levels.
 * Climatological mean vertical profiles (altitude Km)  at selected latitude and month for chemical tracers and time series at selected level and month for chemical tracers.
 * Annual cycle for total ozone (month-lat).
 * Total ozone anomalies at different latitudinal band and seasons.

.. _`Eyring et al. (2006)`: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006JD007327

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_eyring06jgr.yml

Diagnostics are stored in esmvaltool/diag_scripts/eyring06jgr/

    * eyring06jgr_fig01.ncl  
    * eyring06jgr_fig02.ncl
    * eyring06jgr_fig03.ncl
    * eyring06jgr_fig04.ncl
    * eyring06jgr_fig05a.ncl
    * eyring06jgr_fig05b.ncl
    * eyring06jgr_fig07.ncl
    * eyring06jgr_fig08.ncl
    * eyring06jgr_fig09.ncl
    * eyring06jgr_fig12a.ncl
    * eyring06jgr_fig12b.ncl
    * eyring06jgr_fig14.ncl
    * eyring06jgr_fig15.ncl
 


User settings in recipe
-----------------------
#. Preprocessor

   * ``regrid_interp_lev``: Regridding and interpolation reference_dataset levels used by eyring06jgr_fig01, eyring06jgr_fig03 eyring06jgr_fig05, eyring06jgr_fig10, eyring06jgr_fig12a
   * ``regrid_interp_lev_mean_lat``: Regridding and interpolation reference_dataset levels selection of lat/lon area and average  used by eyring06jgr_fig02, eyring06jgr_fig08 and by eyring06jgr_fig09
   * ``regrid_extract_lev``: Regridding on  reference_dataset  grid and extraction of a Pa level used by used by eyring06jgr_fig04, eyring06jgr_fig12b
   * ``regrid_extract_lev_mean_lat`` : Regridding on  reference_dataset  grid and extraction of a Pa level and area selection and mean used by eyring06jgr_fig07
   * ``interp_lev_mean_lat`` : Level interpolation region extraction and mean.
   * ``zonal mean`` : Regridding and zonal mean used by eyring06jgr_fig14 and eyring06jgr_fig15


#. Script <eyring06jgr_fig01.ncl>

   *Required settings for script*

   * ``e06fig01_latmin``: array of float, min lat where variable is averaged, i.e. [ 60., 60., -90., -90. ]
   * ``e06fig01_latmax``: array of float,and max lat where variable is averaged, i.e. [ 90., 90., -60., -60. ]
   * ``e06fig01_season``: array of string., season when variable is averaged, i.e. ["DJF", "MAM","JJA","SON"]
   * ``e06fig01_XMin``: array of float, min limit X axis [-30., -30.,-30.,-30.]
   * ``e06fig01_XMax``: array of float, max limit X axis [20.,20.,20.,20.]
   * ``e06fig01_levmin``: array of float, min limit Y axis [1., 1.,1.,1.]
   * ``e06fig01_levmax``: array of float, max limit Y axis [350., 350.,350.,350.]


   *Optional settings for script*
   * ``e06fig01_start_year``: int,  year when start the climatology calculation [1980] (default max among the models start year).
   * ``e06fig01_end_year``:int, year when end  the climatology calculation [1999] (default min among the models end year).
   * ``e06fig01_multimean``: bool, calculate multi-model mean, (i.e. False/True) (default False).

   *Required settings for variables*
   * ``preprocessor``: regrid_interp_lev.
   * ``reference_dataset``: name of the reference model or observation for regridding and bias calculation (e.g. ERA-Interim").
   *  ``mip``:  Amon.



Variables
---------

*  ta (atmos, monthly mean, longitude latitude level time)



Example plots
-------------

.. _fig_mynewdiag_1:
.. figure::  /recipes/figures/<mynewdiagnostic>/awesome1.png
   :align:   center

   Add figure caption here.
