.. _nml_seaice:

Sea Ice
====================================================

Overview
--------
The sea ice diagnostics plot time series of Arctic and Antarctic sea ice area and extent (calculated as the total area (km2) of grid cells with sea ice concentrations (sic) of at least 15%).

Available recipes and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_seaice.yml

Diagnostics are stored in diag_scripts/seaice/

* seaice_tsline.ncl: creates line plot for sea ice area and extent versus time.
* seaice_aux.ncl: contains a function for calculating sea ice area or extent from sea ice concentration.

User settings in recipe
-----------------------

#. Script seaice_tsline.ncl

   *Required settings (scripts)*

   * region: Arctic, Antarctic
   * month: annual (A), or month number (typically, 3 for March, 9 for September)
    
   *Optional settings (scripts)*
   
   * styleset: for plot_type cycle only (cmip5, cmip6, default)
   * multi_model_mean: plot multi-model mean and standard deviation (default: False)
   * EMs_in_lg: create a legend label for individual ensemble members (default: False)
   * fill_pole_hole: fill polar hole (typically in satellite data) with sic = 1 (default: False)
  

Variables
---------

* sic (ocean-ice, monthly mean, longitude latitude time)
* areacello (fx, longitude latitude)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of cmorization scripts (in esmvaltool/utils/cmorizers/obs) for non-obs4mips data for download instructions.*

* HadISST (sic - esmvaltool/utils/cmorizers/obs/cmorize_obs_HadISST.ncl)


References
----------

* Stroeve, J. et al., Geophys. Res. Lett., 34, L09501, doi:10.1029/2007GL029703, 2007.


Example plots
-------------

