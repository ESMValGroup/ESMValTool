.. _recipe_autoassess_landsurface_surfrad.rst:

Land-surface Surface Radiation - Autoassess diagnostics
=======================================================

Overview
--------

The simulation of surface radiation is central to all aspects of model
performance, and can often reveal compensating errors which are hidden within
top of atmosphere fluxes. This recipe provides metrics that evaluate the skill
of models' spatial and seasonal distribution of surface shortwave and longwave
radiation against the CERES EBAF satellite dataset.

Performance metrics:

* median absolute error (model minus observations) net surface shortwave (SW) radiation
* median absolute error (model minus observations) net surface longwave (LW) radiation

Metrics are calculated using model and observation multi-year climatologies (seasonal means)
for meteorological seasons:
* December-January-February (djf)
* March-April-May (mam)
* June-July-August (jja)
* September-October-November (son)
* Annual mean (ann)


Plots:

* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_surfrad.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_surfrad/surfrad.py: script to calculate surface radiation
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_surfrad to select this diagnostic
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

   * area: must equal land_surface_surfrad to select this diagnostic
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

* rsns (atmos, monthly mean, longitude latitude time)
* rlns (atmos, monthly mean, longitude latitude time)
* sftlf (mask, fixed, longitude latitude)


Observations and reformat scripts
---------------------------------

2001-2012 climatologies (seasonal means) from CERES-EBAF Ed2.7.


References
----------
* Loeb, N. G., D. R. Doelling, H. Wang, W. Su, C. Nguyen, J. G. Corbett, L. Liang, C. Mitrescu, F. G. Rose, and S. Kato, 2018: Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) Top-of-Atmosphere (TOA) Edition-4.0 Data Product. J. Climate, 31, 895-918, doi: 10.1175/JCLI-D-17-0208.1.

* Kato, S., F. G. Rose, D. A. Rutan, T. E. Thorsen, N. G. Loeb, D. R. Doelling, X. Huang, W. L. Smith, W. Su, and S.-H. Ham, 2018: Surface irradiances of Edition 4.0 Clouds and the Earth's Radiant Energy System (CERES) Energy Balanced and Filled (EBAF) data product, J. Climate, 31, 4501-4527, doi: 10.1175/JCLI-D-17-0523.1



Example plots
-------------

.. figure:: /recipes/figures/autoassess_landsurface/Surfrad_Metrics.png
   :scale: 50 %
   :alt: Surfrad_Metrics.png

   Normalised metrics plot comparing a control and experiment simulation



Inputs and usage
----------------
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
      autoassess_landsurf_surfrad: &autoassess_landsurf_surfrad_settings
        script: autoassess/autoassess_area_base.py
        title: "Autoassess Land-Surface Diagnostic Surfrad Metric"
        area: land_surface_surfrad
        control_model: UKESM1-0-LL
        exp_model: UKESM1-0-LL
        obs_models: [CERES-EBAF]
        obs_type: obs4MIPs
        start: 1997/12/01
        end: 2002/12/01
