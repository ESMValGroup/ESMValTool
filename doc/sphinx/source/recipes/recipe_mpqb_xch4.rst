.. _recipe_mpqb_xch4:

Diagnostics of integrated atmospheric methane (XCH4)
====================================================

Overview
--------

This recipe ``recipe_mpqb_xch4.yml`` allows the comparison of integrated atmospheric methane
between CMIP6 model simulations and observations, and produces lineplots of monthly mean
methane values, annual cycles and annual growth rates:

* Monthly mean time series of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)
* Annual cycles of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)
* Annual growth rates of XCH4 for pre-defined regions (global, Northern Hemisphere, Southern Hemisphere)

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/mpqb/

* recipe_mpqb_xch4.yml

Diagnostics are stored in esmvaltool/diag_scripts/mpqb/

* mpqb_lineplot.py
* mpqb_lineplot_anncyc.py
* mpqb_lineplot_growthrate.py

Observations and reformat scripts
---------------------------------

Observations used in this recipe are:

* CDS-XCH4 (ESA CCI dataset served on the Copernicus Climate data store)

A cmorizing script for this dataset is available (``cmorize_obs_cds_xch4.ncl``).

XCH4 is a derived variable that needs to be calculated from four different variables (ch4, hus, zg, ps).
A derivation script is included in the ESMValCore.


User settings in recipe
-----------------------
#. Preprocessor

   * ``pp_lineplots_xx_mon``: Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv] to obtain one time series of monthly mean values for the selected region (global, southern hemisphere, northern hemisphere)
   * ``pp_lineplots_xx_ann``: Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv] to obtain one time series of annual mean values for the selected region (global, southern hemisphere, northern hemisphere)   
   * ``pp_lineplots_anncyc_xx:`` : Regridding, masking all missing values from all used datasets, area-mean ('xx' can ge replaced by 'gl'=global, 'sh'=southern hemisphere, 'nh'=northern hemisphere), units converted to [ppbv], monthly climate statistics applied to one annual cycle for the whole chosen time period and for the selected region (global, southern hemisphere, northern hemisphere)
   * ``xch4_def_xx``: defining the time period over which the analysis should be calculated; options are "cmip6" which overlapping period of the observations and the CMIP6 historical simulations, and "future" which covers the time period of CMIP6 scenarios

#. Additional needed files
   
   * ``mpqb_cfg_xch4.yml``: In this file additional information for the used datasets are defined and stored, e.g. alias of the dataset name and the color that is used to display the dataset in the figures
   * ``mpqb_utils.yml``: In this file the preparations for the dataset displays are made.

#. Script <mpqb_lineplot.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required

#. Script <mpqb_lineplot_anncyc.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required

#. Script <mpqb_lineplot_growthrate.py>

   *Required settings for script*

   * no additional settings required

   *Optional settings for script*
   
   * no optional settings available

   *Required settings for variables*
   
   * no settings for the variable required


Variables
---------

* ch4 (atmos, monthly mean, longitude latitude level time)
* hus (atmos, monthly mean, longitude latitude level time)
* zg (atmos, monthly mean, longitude latitude level time)
* ps (atmos, monthly mean, longitude latitude time)

All variables are necessary to calculate the derived variable xch4.


Example plots
-------------

.. _lineplot_xch4_2003-2014_monmean: 
.. figure::  /recipes/figures/mpqb/lineplot_xch4_2003-2014_monmean.png
   :align:   center

   Monthly mean time series of XCH4, calculated over the whole globe, for individual CMIP6 model simulations.
   
   
   
