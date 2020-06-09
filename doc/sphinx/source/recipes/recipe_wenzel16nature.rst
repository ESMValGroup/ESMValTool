.. _recipes_wenzel16nature:

Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO2
=============================================================================================

Overview
--------

Figures from `Wenzel et al. (2016)`_ are reproduced with recipe_wenzel16nature.yml. Gross primary productivity (gpp) and atmospheric CO2 concentrations at the surface  (co2s) are analyzed for the carbon cycle - concentration feedback in the historical (esmHistorical) and uncoupled (esmFixCLim1, here the carbon cycle is uncoupled to the climate response) simulations. The standard namelist includes a set of routines to diagnose the long-term carbon cycle - concentration feedback parameter (beta) from an ensemble of CMIP5 models and the observable change in the [CO2] seasonal cycle amplitude due to rising atmospheric CO2 levels. As a key figure of this recipe, the diagnosed values from the models beta vs. the change in CO2 amplitude are compared in a scatter plot constituting an emergent constraint.

.. _`Wenzel et al. (2016)`: https://www.nature.com/articles/nature19772

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

* co2s (atmos, monthly mean, plev longitude latitude time)
* gpp (land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* ESRL: Earth System Research Laboratory, ground mased co2 measuremends


References
----------

* Wenzel, S., Cox, P., Eyring, V. et al., 2016, Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO2. Nature 538, 499501, doi: doi.org/10.1038/nature19772


Example plots
-------------

.. figure:: /recipes/figures/wenzel16nature/
   :width: 10 cm 
   :align: center
   
   XXXX Comparison of CO2 seasonal amplitudes for CMIP5 historical simulations and observations showing Annual mean atmospheric CO2 versus the amplitudes of the CO2 seasonal cycle at Pt. Barrow, Alaska 
      
.. figure:: /recipes/figures/wenzel16nature/
   :width: 10 cm 
   :align: center
   
   XXXX Histogram showing the gradient of the linear correlations for the comparison of CO2 seasonal amplitudes for CMIP5 historical for at Pt. Barrow, Alaska 

.. figure:: /recipes/figures/wenzel16nature/
   :scale: 50 %
   :align: center

   XXXX Emergent constraints on the relative increase of large-scale GPP for a doubling of CO2, showing the correlations between the sensitivity of the CO2 amplitude to annual mean CO2 increases at Pt. Barrow (x axis) and the high-latitude (60 N - 90 N) CO2 fertilization on GPP at 2 x CO2. The red line shows the linear best fit of the regression together with the prediction error (orange shading) and the gray shading shows the observed range.
      
