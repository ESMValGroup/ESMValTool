Emergent constraint on carbon cycle concentration feedback
===========================================================

Overview
--------

Figures from Wenzel et al. (2016) are reproduced with recipe_wenzel16nature.yml. Gross primary productivity (gpp) and atmospheric CO2 (co2) at the surface are analyzed for the carbon cycle - concentration feedback in the historical (esmHistorical) and uncoupled (esmFixCLim1, here the carbon cycle is uncoupled to the climate response) simulations. The standard namelist includes a set of routines to diagnose the long-term carbon cycle - climate feedback parameter (beta) from an ensemble of CMIP5 models. Also included in the recipe is a comparison of the co2 seasonal cycle amplitude for historical simulations used to diagnose the observable sensitivity of the seasonal cycle amplitude to rising [CO2] levels. As a key figure of this recipe, the diagnosed values from the models beta vs. sensitivity of co2 amplitude are compared in a scatter plot constituting an emergent constraint.


Available recipe and diagnostics
-----------------------------------

Recipes are stored in recipes/

* recipe_wenzel16nature.yml

Diagnostics are stored in diag_scripts/

* carbon_beta: scatter plot of ; diagnosing beta for later use in EC plot
* carbon_co2_cycle.ncl: scatter plot of time line plots of annual means for spatial averages


User settings
-------------

User setting files (cfg files) are stored in nml/cfg_carbon/

#. carbon_beta 

   *Required Settings (scripts)*

   * styleset: CMIP5

   *Optional Settings (scripts)*

   * cl_mean: if true calculates mean of beta
   * bc_xmax_year: end year to calculate beta, else end_year is used
   * bc_xmin_year: start year to calculate beta, else start_year is used

   *Required settings (variables)*

   * reference_dataset: name of reference data set

#. carbon_co2_cycle.ncl 

   *Required Settings (scripts)*

   * bc_xmax_year: same end year as in carbon_beta
   * bc_xmin_year: same start year as in carbon_beta
   * styleset: CMIP5
   * nc_infile: path were to find file with beta values


Variables
---------

* co2 (atmos, monthly mean, plev longitude latitude time)
* gpp (land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* ESRL: Earth System Research Laboratory, ground mased co2 measuremends


References
----------

* Wenzel, S., Cox, P., Eyring, V. et al., 2016, Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO2. Nature 538, 499501, doi: doi.org/10.1038/nature19772


Example plots
-------------

.. figure:: /recipes/figures/wenzel14jgr/tas_Global_CMIP5_1pctCO2_anom__1-1999.png
   :width: 10 cm 
   :align: center
   
   XXXX Time series of tropical (30S to 30N) mean near surface temperature (tas) change between year 30 and year 110 for the CMIP5 models simulated with prescribed CO2 (1%/yr CO2 increase) coupled simulation (1pctCO2).
   
   
.. figure:: /recipes/figures/wenzel14jgr/corr_tas-nbp_anom_1960-2005.png
   :width: 10 cm 
   :align: center
   
   XXXX Correlations between the interannual variability of global co2flux (nbp+fgco2) and tropical temperature for the individual CMIP5 models using esmHistorical simulations, and for observations.


.. figure:: /recipes/figures/wenzel14jgr/constr_tas-nbp_30-1960.000001.png
   :scale: 50 %
   :align: center

   XXXX Carbon cycle-climate feedback of tropical land carbon vs. the sensitivity of co2flux to interannual temperature variability in the tropics (30S to 30N). The red line shows the linear best fit of the regression together with the prediction error (orange shading) and the gray shading shows the observed range.
      
