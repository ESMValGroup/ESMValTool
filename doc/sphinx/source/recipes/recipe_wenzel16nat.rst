.. _recipes_wenzel16nat:

Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO\ :sub:`2`
======================================================================================================

Overview
--------

Selected figures from `Wenzel et al. (2016)`_ are reproduced with recipe_wenzel16nat.yml. Gross primary productivity (gpp) and atmospheric CO\ :sub:`2` concentrations at the surface  (co2s) are analyzed for the carbon cycle - concentration feedback in the historical (esmHistorical) and uncoupled (esmFixCLim1, here the carbon cycle is uncoupled to the climate response) simulations. The recipe includes a set of routines to diagnose the long-term carbon cycle - concentration feedback parameter (beta) from an ensemble of CMIP5 models and the observable change in the CO\ :sub:`2` seasonal cycle amplitude due to rising atmospheric CO\ :sub:`2` levels. As a key figure of this recipe, the diagnosed values from the models beta vs. the change in CO\ :sub:`2` amplitude are compared in a scatter plot constituting an emergent constraint.

.. _`Wenzel et al. (2016)`: https://www.nature.com/articles/nature19772

Available recipe and diagnostics
-----------------------------------

Recipes are stored in recipes/

    * recipe_wenzel16nat.yml

Diagnostics are stored in diag_scripts/carbon_ec/

    * carbon_beta: (1) scatter plot of annual gpp vs. annual CO\ :sub:`2` and
      (2) barchart of gpp(2xCO\ :sub:`2`)/gpp(1xCO\ :sub:`2`); calculates beta
      for emergent constraint (carbon_co2_cycle.ncl)
    * carbon_co2_cycle.ncl: (1) scatter plot of CO\ :sub:`2` amplitude vs.
      annual CO\ :sub:`2`, (2) barchart of sensitivity of CO\ :sub:`2` amplitude
      to CO\ :sub:`2`, (3) emergent constraint:
      gpp(2xCO\ :sub:`2`)/gpp(1xCO\ :sub:`2`) vs. sensitivity of CO\ :sub:`2`
      amplitude to CO\ :sub:`2`, (4) probability density function of constrained
      and unconstrained sensitivity of CO\ :sub:`2` amplitude to CO\ :sub:`2`


User settings
-------------

#. Script carbon_beta.ncl

   *Required Settings (scripts)*

   * styleset: project style for lines, colors and symbols

   *Optional Settings (scripts)*

   * bc_xmax_year: end year to calculate beta (default: use last available year of all models)
   * bc_xmin_year: start year to calculate beta (default: use first available year of all models)

   *Required settings (variables)*

   none

   *Optional settings (variables)*

   none

#. Script carbon_co2_cycle.ncl 

   *Required Settings (scripts)*

   * nc_infile: path of netCDF file containing beta (output from carbon_beta.ncl)
   * styleset: project style for lines, colors and symbols

   *Optional Settings (scripts)*

   * bc_xmax_year: end year (default = last year of all model datasets available)
   * bc_xmin_year: start year (default = first year of all model datasets available)

   *Required settings (variables)*

   * reference_dataset: name of reference datatset (observations)

   *Optional settings (variables)*

   none


Variables
---------

* co2s (atmos, monthly mean, plev longitude latitude time)
* gpp (land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

* ESRL: Earth System Research Laboratory, ground-based CO\ :sub:`2` measurements


References
----------

* Wenzel, S., Cox, P., Eyring, V. et al., 2016, Projected land photosynthesis constrained by changes in the seasonal cycle of atmospheric CO\ :sub:`2`. Nature 538, 499501, doi: doi.org/10.1038/nature19772


Example plots
-------------

.. figure:: /recipes/figures/wenzel16nat/fig_1.png
   :width: 12 cm 
   :align: center
   
   Comparison of CO\ :sub:`2` seasonal amplitudes for CMIP5 historical simulations and observations showing annual mean atmospheric CO\ :sub:`2` versus the amplitudes of the CO\ :sub:`2` seasonal cycle at Pt. Barrow, Alaska (produced with carbon_co2_cycle.ncl, similar to Fig. 1a from Wenzel et al. (2016)).
      
.. figure:: /recipes/figures/wenzel16nat/fig_2.png
   :width: 12 cm 
   :align: center
   
   Barchart showing the gradient of the linear correlations for the comparison of CO\ :sub:`2` seasonal amplitudes for CMIP5 historical for at Pt. Barrow, Alaska (produced with carbon_co2_cycle.ncl, similar to Fig. 1b from Wenzel et al. (2016)).

.. figure:: /recipes/figures/wenzel16nat/fig_3.png
   :width: 12 cm
   :align: center

   Emergent constraint on the relative increase of large-scale GPP for a doubling of CO\ :sub:`2`, showing the correlations between the sensitivity of the CO\ :sub:`2` amplitude to annual mean CO\ :sub:`2` increases at Pt. Barrow (x-axis) and the high-latitude (60N - 90N) CO\ :sub:`2` fertilization on GPP at 2xCO\ :sub:`2`. The red line shows the linear best fit of the regression together with the prediction error (orange shading), the gray shading shows the observed range (produced with carbon_co2_cycle.ncl, similar to Fig. 3a from Wenzel et al. (2016)).
