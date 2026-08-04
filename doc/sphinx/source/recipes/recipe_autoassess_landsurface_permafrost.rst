.. _recipe_autoassess_landsurface_permafrost.rst:

Land-surface Permafrost - Autoassess diagnostics
================================================

Overview
--------

Permafrost thaw is an important impact of climate change, and is the source of
a potentially strong Earth system feedback through the release of soil carbon
into the atmosphere. This recipe provides metrics that evaluate the
climatological performance of models in simulating soil temperatures that
control permafrost. Performance metrics (with observation-based estimates in brackets):

* permafrost area (17.46 million square km)
* fractional area of permafrost northwards of zero degree isotherm (0.47)
* soil temperature at 1m minus soil temperature at surface (-0.53 degrees C)
* soil temperature at surface minus air temperature (6.15 degrees C)
* annual amplitude at 1m / annual amplitude at the surface (0.40 unitless)
* annual amplitude at the surface / annual air temperature (0.57 unitless)


Plots:

* Maps of permafrost extent and zero degC isotherm
* Normalised assessment metrics plot comparing control and experiment

The recipe takes as input a control model and experimental model, comparisons being made
with these two models.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_autoassess_landsurface_permafrost.yml

Diagnostics are stored in esmvaltool/diag_scripts/autoassess/

    * autoassess_area_base.py: wrapper for autoassess scripts
    * land_surface_permafrost/permafrost.py: script to calculate permafrost
      metrics
    * plot_autoassess_metrics.py: plot normalised assessment metrics


User settings in recipe
-----------------------

#. Script autoassess_area_base.py

   *Required settings for script*

   * area: must equal land_surface_permafrost to select this diagnostic
   * control_model: name of model to be used as control
   * exp_model: name of model to be used as experiment
   * start: date (YYYY/MM/DD) at which period begins (see note on time gating)
   * end: date (YYYY/MM/DD) at which period ends (see note on time gating)

   *Optional settings for script*

   * title: arbitrary string with name of diagnostic
   * obs_models: unused for this recipe

   *Required settings for variables*

   none

   *Optional settings for variables*

   none


#. Script plot_autoassess_metrics.py

   *Required settings for script*

   * area: must equal land_surface_permafrost to select this diagnostic
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

* tas (atmos, monthly mean, longitude latitude time)
* tsl (land, monthly mean, longitude latitude time)
* mrsos (land, monthly mean, longitude latitude time)
* sftlf (mask, fixed, longitude latitude)


Observations and reformat scripts
---------------------------------

None


References
----------

* Observed permafrost extent is from https://nsidc.org/data/ggd318/versions/2: Brown, J.,
  O. Ferrians, J. A. Heginbottom, and E. Melnikov. 2002. Circum-Arctic Map of
  Permafrost and Ground-Ice Conditions, Version 2. Boulder, Colorado USA. NSIDC:
  National Snow and Ice Data Center.  When calculating the global area of
  permafrost the grid cells are weighted by the proportion of permafrost within
  them.

* Annual mean air temperature is from: Legates, D. R., and C. J. Willmott, 1990:
  Mean seasonal and spatial variability in global surface air temperature. Theor.
  Appl. Climatol., 41, 11-21.  The annual mean is calculated from the seasonal
  mean data available at the Met Office.

* The soil temperature metrics are calcuated following: Charles D. Koven, William
  J. Riley, and Alex Stern, 2013: Analysis of Permafrost Thermal Dynamics and
  Response to Climate Change in the CMIP5 Earth System Models. J. Climate, 26.
  (Table 3) http://dx.doi.org/10.1175/JCLI-D-12-00228.1 The
  locations used for Table 3 were extracted from the model and the modelled
  metrics calculated.


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
The ``landsurface_permafrost`` area metric is part of the ``esmvaltool/diag_scripts/autoassess`` diagnostics,
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
      plot_landsurf_permafrost: &plot_landsurf_permafrost_settings
        <<: *autoassess_landsurf_permafrost_settings
        control_model: MPI-ESM-LR
        exp_model: MPI-ESM-MR
        script: autoassess/plot_autoassess_metrics.py
        ancestors: ['*/autoassess_landsurf_permafrost']
        title: "Plot Land-Surface Permafrost Metrics"
        plot_name: "Permafrost_Metrics"
        diag_tag: aa_landsurf_permafrost
        diag_name: autoassess_landsurf_permafrost
