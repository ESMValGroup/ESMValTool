.. _recipe_climate_change_hotspot.rst:

Climate Change Hotspot
======================

Overview
--------



Plots:



Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * _recipe_climate_change_hotspot.yml

Diagnostics are stored in esmvaltool/diag_scripts/cos22esd/

    * .py:
    * .py: descripció

User settings in recipe
-----------------------

#. Script .py

   *Required settings for script*

   script: cos22esd/climate_change_hotspot.py
   baseline_period: &baseline [1986, 2005]
   future_periods: &future ["2041-2060", "2081-2100"]
   region: &region [-10, 40, 30, 45]
   region_name:



#. Script .py

   *Required settings for script*

   * area: must equal land_surface_permafrost to select this diagnostic
   * control_model: name of model to be used as control in metrics plot
   * exp_model: name of model to be used as experiment in metrics plot
   * title: string to use as plot title



Variables
---------

* tas (atmos, monthly mean, longitude latitude time)
* pr (atmos, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

None


References
----------

* cos2022


Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_north_america_ACCESS-CM2.png
   :scale: 50 %
   :alt: pf_extent_north_america_ACCESS-CM2.png

   Permafrost extent and zero degC isotherm, showing North America

.. figure:: /recipes/figures/autoassess_landsurface/pf_extent_asia_ACCESS-CM2.png
   :scale: 50 %
   :alt: pf_extent_asia_ACCESS-CM2.png

   Permafrost extent and zero degC isotherm, showing Asia and Europe

.. figure:: /recipes/figures/autoassess_landsurface/Permafrost_Metrics.png
   :scale: 50 %
   :alt: Permafrost_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation


Additional notes on usage
-------------------------
region
