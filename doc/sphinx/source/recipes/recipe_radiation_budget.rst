.. _recipes_radiation_budget:

Radiation Budget
================

Overview
--------

The aim of monitoring the energy budget is to understand the (im)balance
of energy flux between the atmosphere and the surface of a model due to its
link with the hydrological cycle and climate change.

This diagnostic analyses the radiation budget by separating top-of-atmosphere
fluxes into clear-sky and cloud forcing components, and surface fluxes into
downwelling and upwelling components. Model predictions are compared against
three observational estimates, one of which (Stephens et al. 2012) includes
uncertainty estimates. When the black error bars overlap the zero line, the
model is consistent with observations according to Stephens et al. (2012).

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_radiation_budget.yml

Diagnostics are stored in esmvaltool/diag_scripts/radiation_budget/

    * radiation_budget.py: Plot the global radiation budget.
    * seasonal_radiation_budget.py: Write the global climatological seasonal radiation budget to a text file.



User settings in recipe
-----------------------

None


Variables
---------

* rss (atmos, monthly mean, longitude latitude time)
* rsdt (atmos, monthly mean, longitude latitude time)
* rsut (atmos, monthly mean, longitude latitude time)
* rsutcs (atmos, monthly mean, longitude latitude time)
* rsds (atmos, monthly mean, longitude latitude time)
* rls (atmos, monthly mean, longitude latitude time)
* rlut (atmos, monthly mean, longitude latitude time)
* rlutcs (atmos, monthly mean, longitude latitude time)
* rlds (atmos, monthly mean, longitude latitude time)
* hfss (atmos, monthly mean, longitude latitude time)
* hfls (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* CERES-EBAF (rlut, rlutcs, rsut, rsutcs - obs4MIPs)
* Demory observations can be found in esmvaltool/diag_scripts/radiation_budget/Demory_et_al_2014_obs_Energy_Budget.yml and are from Figure 2 in Demory et al. (2014).
* Stephens observations can be found in esmvaltool/diag_scripts/radiation_budget/Stephens_et_al_2012_obs_Energy_Budget.yml from figure 1b in Stephens et al. (2012).


References
----------

* Demory, ME., Vidale, P.L., Roberts, M.J. et al. The role of horizontal resolution in simulating drivers of the global hydrological cycle. Clim Dyn 42, 2201–2225 (2014). https://doi.org/10.1007/s00382-013-1924-4
* Stephens, G., Li, J., Wild, M. et al. An update on Earth's energy balance in light of the latest global observations. Nature Geosci 5, 691–696 (2012). https://doi.org/10.1038/ngeo1580


Example plots
-------------

.. _fig_radiation_budget_1:
.. figure::  /recipes/figures/radiation_budget/UKESM1-0-LL.png
   :align:   center

   Radiation budget for UKESM1-0-LL
