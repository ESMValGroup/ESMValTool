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

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_soilmoisture/soilmoisture.py: script to calculate soil moisture
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_soilmoisture to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment
   * start: date (YYYY/MM/DD) at which period begins (see note on time gating)
   * end: date (YYYY/MM/DD) at which period ends (see note on time gating)
   * climfiles_root: path to observation climatologies

   *Optional settings for script*

   * title: arbitrary string with name of diagnostic
   * obs_models: unused for this recipe

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

* mrsos (land, monthly mean, longitude latitude time)


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


Additional notes on usage
-------------------------
The ``landsurface_soilmoisture`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
and, as any other ``autoassess`` metric, it uses the ``autoassess_area_base.py`` as general purpose
wrapper. This wrapper accepts a number of input arguments that are read through from the recipe.

This recipe is part of the larger group of Autoassess metrics ported to ESMValTool
from the native Autoassess package from the UK's Met Office. The ``diagnostics`` settings
are almost the same as for the other Autoassess metrics.

.. note::

   **Time gating for autoassess metrics.**

   To preserve the native Autoassess functionalities,
   data loading and selection on time is done somewhat
   differently for ESMValTool's autoassess metrics: the
   time selection is done in the preprocessor as per usual but
   a further time selection is performed as part of the diagnostic.
   For this purpose the user will specify a ``start:`` and ``end:``
   pair of arguments of ``scripts: autoassess_script`` (see below
   for example). These are formatted as ``YYYY/MM/DD``; this is
   necessary since the Autoassess metrics are computed from 1-Dec
   through 1-Dec rather than 1-Jan through 1-Jan. This is a temporary
   implementation to fully replicate the native Autoassess functionality
   and a minor user inconvenience since they need to set an extra set of
   ``start`` and ``end`` arguments in the diagnostic; this will be phased
   when all the native Autoassess metrics have been ported to ESMValTool
   review has completed.


An example of standard inputs as read by ``autoassess_area_base.py`` and passed
over to the diagnostic/metric is listed below.


.. code-block:: yaml

    scripts:
      autoassess_landsurf_soilmoisture: &autoassess_landsurf_soilmoisture_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Soilmoisture Diagnostic"
        area: land_surface_soilmoisture
        control_model: IPSL-CM5A-LR
        exp_model: inmcm4
        obs_models: []
        start: 1997/12/01
        end: 2002/12/01
        climfiles_root: '/gws/nopw/j04/esmeval/autoassess_specific_files/files'  # on JASMIN
