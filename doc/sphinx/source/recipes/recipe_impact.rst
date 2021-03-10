.. _recipes_impact:

Quick insights for climate impact researchers
=============================================

Overview
--------

Many impact researchers do not have the time and finances to use a large
ensemble of climate model runs for their impact analysis. To get an idea of the
range of impacts of climate change it also suffices to use a small number of
climate model runs. In case a system is only sensitive to annual temperature,
one can select a run with a high change and one with a low change of annual
temperature, preferably both with a low bias.

This recipe calculates the bias with respect to observations, and the change
with respect to a reference period, for a wide range of (CMIP) models. These
metrics are tabulated and also visualized in a diagram.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_impact.yml

Diagnostics are stored in esmvaltool/diag_scripts/<mynewdiag>/

    * impact/quick_insights.py: tabulate and visualize bias and change.


User settings in recipe
-----------------------

#. Script ``impact.py``

   *Required settings for variables*

   * tag: 'model' or 'observation', so the diagnostic script knows which datasets to use for the bias calculation. This must be specified for each dataset.

   *Optional settings for preprocessor*

   * Region and time settings (both for the future and reference period) can be changed at will.


Variables
---------

* tas (atmos, mon, longitude latitude time)
* pr (atmos, mon, longitude latitude time)
* any other variables of interest


Observations and reformat scripts
---------------------------------

* ERA5 data can be used via the native6 project.

References
----------

* None

Example plots
-------------

.. _fig_impact_1:
.. figure::  /recipes/figures/impact/bias_vs_change.png
   :align:   center

   "Bias and change for each variable"
