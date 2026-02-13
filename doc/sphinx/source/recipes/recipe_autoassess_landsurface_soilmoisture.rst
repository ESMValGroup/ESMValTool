.. _recipe_autoassess_landsurface_soilmoisture.rst:

Land-surface Soil Moisture - Autoassess diagnostics
===================================================

Overview
--------

Soil moisture is a critical component of the land system, controling surface
energy fluxes in many areas of the world. This recipe provides metrics that
evaluate the skill of models' spatial and seasonal distribution of soil
moisture against the ESA CCI soil moisture ECV.

Performance metrics:

* median absolute error (model minus observations)

Metrics are calculated using model and observation multi-year climatologies (seasonal means)
for meteorological seasons:

* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)

Plots:

* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_soilmoisture.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * land_surface_soilmoisture/soilmoisture.py: script to calculate soil moisture
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script soilmoisture.py

   *Required settings for script*

   * area: must equal land_surface_soilmoisture to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment

   *Optional settings for script*

   none

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


#. Script plot_autoassess_metrics.py

   *Required settings for script*

   * area: must equal land_surface_soilmoisture to select this diagnostic
   * control_model: name of model to be used as control in metrics plot
   * exp_model: name of model to be used as experiment in metrics plot
   * title: string to use as plot title

   *Optional settings for script*

   none

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


Variables
---------

* mrsos (from models: land, monthly mean, longitude latitude time)
* sm (from observations: land, monthly mean, longitude latitude time)


Observations and reformat scripts
---------------------------------

1999-2008 climatologies (seasonal means) from ESA ECV Soil Moisture Dataset v1.
Produced by the ESA CCI soil moisture project: https://www.esa-soilmoisture-cci.org/node/93


References
----------
* Dorigo, W.A., Wagner, W., Albergel, C., Albrecht, F.,  Balsamo, G., Brocca, L., Chung, D., Ertl, M., Forkel, M., Gruber, A., Haas, E., Hamer, D. P. Hirschi, M., Ikonen, J., De Jeu, R. Kidd, R.  Lahoz, W., Liu, Y.Y., Miralles, D., Lecomte, P. (2017).  ESA CCI Soil Moisture for improved Earth system understanding: State-of-the art and future directions. In Remote Sensing of Environment, 2017,  ISSN 0034-4257, https://doi.org/10.1016/j.rse.2017.07.001.

* Gruber, A., Scanlon, T., van der Schalie, R., Wagner, W., Dorigo, W. (2019). Evolution of the ESA CCI Soil Moisture Climate Data Records and their underlying merging methodology. Earth System Science Data 11, 717-739, https://doi.org/10.5194/essd-11-717-2019


Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/Soilmoisture_Metrics.png
   :scale: 50 %
   :alt: Soilmoisture_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation
