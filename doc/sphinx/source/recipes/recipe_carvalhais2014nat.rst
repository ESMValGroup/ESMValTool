.. _recipe_carvalhais2014nat:

Turnover time of carbon over land
====================================================

Overview
--------

This recipe evaluates the following aspects of turnover time of carbon over
land (tau_ctotal) as in `Carvalhais et al. (2014)` :

* Global distributions of tau_ctotal from all models against observations and
  other models
* Variation of tau_ctotal along the latitude (zonal distributions)
* Variation of association of tau_ctotal and climate (zonal correlations)

Additionally, metrics of global turnover times and correlations with
observation are plotted over the maps/global distributions.

.. _tau calculation:

Calculation of turnover time
----------------------------

First, the total carbon content over the land surface is calculated as,
.. math::

{ctotal} =  {cSoil + cVeg}

where :math:`cSoil` and :math:`cVeg` are the carbon content of the soil and
vegetation over the land surface. **Note that this may be inconsistent with the
original publication, as all the carbon pools which respired to the atmosphere
were added up to calculate the total carbon content over land surface**.


The turnover time of carbon is then calculated as,
.. math::

\tau_{ctotal} =  \frac{ctotal}{gpp}


where `gpp` is the gross primary productivity. **Note that the above equation
is valid for steady state, and is only applicable when both the ctotal and gpp
are long-term averages.** That is why the recipe includes the preprocessor for
calculating the long term  averages from the monthly time series.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * _recipe_carvalhais2014nat.yml


Diagnostics are stored in diag_scripts/

   * land_carbon_cycle/diag_global_turnover.py
   * land_carbon_cycle/diag_zonal_turnover.py
   * land_carbon_cycle/diag_zonal_correlation.py

.. _user settings:

User settings in recipe
-----------------------

#. Observation related settings

    * obs_data_subdir: 'data_ESMValTool_Carvalhais2014'
    * source_label: 'Carvalhais2014'
    * variant_label: 'BE'
    * grid_label: 'gn'
       *  'gr1': 2.812x2.813, bcc-csm1-1
       *  'gr2': 1.25x0.937, CCSM4
       *  'gr3': 2.812x2.813, CanESM2
       *  'gr4': 2.5x2.0, GFDL-ESM2G
       *  'gr5': 1.875x1.241, HadGEM2-ES
       *  'gr6': 2.0x1.5, inmcm4
       *  'gr7': 2.5x1.259, IPSL-CM5A-MR
       *  'gr8': 2.812x2.813, MIROC-ESM
       *  'gr9': 1.875x1.875, MPI-ESM-LR
       *  'gr': 2.5x1.875, NorESM1-M
    * frequency: 'fx'

#. Preprocessor

   * ``climate_statistics``: mean over full time period
   * ``regrid``: nearest neighbor regridding to the resolution of the
     observation data
   * ``mask_landsea``: mask out all the data points from sea
   * ``multi_model_statistics``: calculate the multimodel median


#. Script land_carbon_cycle/diag_global_turnover.py

  * Required settings:

    * ``obs_variable``: list of the variables to be read from the observation
      files

  * Optional settings:

    * ``ax_fs``: fontsize for the figure (default: 7.1)
    * ``fill_value``: fill value to be used in analysis (default: NaN)
    * ``x0``: X - coordinate of the lower left corner of the matrix of maps
      (between 0 and 1)
    * ``y0``: Y - coordinate of the lower left corner of the matrix of maps
      (between 0 and 1)
    * ``wp``: width of each map
    * ``hp``: height of each map
    * ``xsp``: spacing betweeen maps in X - direction
    * ``ysp``: spacing between maps in Y -direction
    * ``aspect_map``: aspect of the map
    * ``xsp_sca``: spacing between the scatter plots in X - direction
    * ``ysp_sca``: spacing between the scatter plots in X - direction
    * ``hcolo``: height of the colorbar (thickness for horizontal orientation)
    * ``wcolo``: width of the colorbar (length)
    * ``cb_off_y``: distance of colorbar from top the maps
    * ``x_colo_d``: x-coordinate of the colorbar for maps along the diagonal
    * ``x_colo_r``: x-coordinate of the colorbar for ratio maps above the
      diagonal
    * ``y_colo_single``: y-coordinate of the colorbar in the maps per model
    * ``correlation_method``: correlation method to be used while calculating
      the correlation displayed in the scatter plots ({spearman} | pearson)
    * ``tx_y_corr``: relative y-coordinate of the text of correlation
    * ``valrange_sc``: The range of values for x- and y- axes limits in the
      scatterplot
    * ``obs_global``: global turnover time, provided as additional info for the
      map of the observation. {23}. For models, they are calculated by the
      diagnostic.
    * ``gpp_threshold``: The threshold of gpp in gC m-2 yr-1 for consideration
      in the calculation. A default value of 10, as in the observation, is set.

#. Script land_carbon_cycle/diag_zonal_turnover.py

  * Required settings:

    * ``obs_variable``: list of the variables to be read from the observation
      files

  * Optional settings:

    * ``ax_fs``: fontsize for the figure (default 7.1)
    * ``multimodel``: False
    * ``fill_value``: fill value to be used in analysis (default``: NaN)
    * ``min_points_frac``:
    * ``valrange_x``: (2, 256)
    * ``valrange_y``: (2, 256)
    * ``bandsize``: 1.0
    * ``gpp_threshold``: 10  in gC m-2 yr -1


#. Script land_carbon_cycle/diag_zonal_correlation.py

  * Required settings:

    * ``obs_variable``: list of the variables to be read from the observation
      files

  * Optional settings:

    * ``ax_fs``: fontsize for the figure (default 7.1)
    * ``multimodel``: False
    * ``fill_value``: fill value to be used in analysis (default: NaN)
    * ``correlation_method``: 'spearman'
    * ``min_points_frac``:
    * ``valrange_x``: (2, 256)
    * ``valrange_y``: (2, 256)
    * ``bandsize``: 1.0
    * ``gpp_threshold``: 10  #gC m-2 yr -1


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)
* *pr* (atmos, monthly, longitude, latitude, time)
* *gpp* (land, monthly, longitude, latitude, time)
* *cVeg* (land, monthly, longitude, latitude, time)
* *cSoil* (land, monthly, longitude, latitude, time)


Observations
------------

The observations needed in the diagnostics are publicly available for download
from the .. _Data Portal of the Max Planck Institute for Biogeochemistry:
https://www.bgc-jena.mpg.de/geodb/projects/Home.php.

Due to inherent dependence of the diagnostic on the uncertainty estimates,
the data needed for each diagnostic script are provided as preprocessed netCDF
files.



References
----------

* Carvalhais, N., et al. (2014), Global covariation of carbon turnover times
  with climate in terrestrial ecosystems, Nature, 514(7521), 213-217,
  doi: 10.1038/nature13731.


Example plots
-------------

.. _fig_carvalhais2014nat_1:
.. figure:: /recipes/figures/carvalhais2014nat/comparison_zonal_pearson_correlation_turnovertime_climate_Carvalhais2014.png
   :align: center
   :width: 80%

   Time series of global net biome productivity (NBP) over the period
   1901-2005. Similar to Anav et al.  (2013), Figure 5.

.. _fig_carvalhais2014nat_2:
.. figure:: /recipes/figures/carvalhais2014nat/global_comparison_matrix_models_Carvalhais2014.png
   :align: center
   :width: 80%

   Seasonal cycle plot for nothern hemisphere gross primary production (GPP)
   over the period 1986-2005. Similar to Anav et al. (2013), Figure 9.
