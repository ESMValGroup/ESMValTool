.. _broken-recipe-list:

Broken recipe list
==================

This table gives an overview of the recipes that are known to have issues.
The table is always valid for the latest stable release of ESMValTool.
More details can be found in the :ref:`broken recipe policy
<broken-recipe-policy>`.

.. list-table:: Broken recipes
   :widths: 25 25 25 25
   :header-rows: 1

   * - Broken recipe
     - Affected diagnostics
     - Problem
     - GitHub issue
   * - `recipe_check_obs.yml`
     - `ERA5_native6`
     - Derivation of custom variables `rlus` and `rsus`
     - `#1388 <https://github.com/ESMValGroup/ESMValCore/issues/1388>`_
   * - :ref:`recipe_julia.yml <recipe_examples>`
     - `example`
     - fill values are not interpreted, resulting in an unusable plot
     - `#2595 <https://github.com/ESMValGroup/ESMValTool/issues/2595>`_
   * - :ref:`recipe_seaice_drift.yml <recipes_seaice_drift>`
     - `sea_ice_drift_SCICEX`
     - ``shapely 2`` issue
     - `#3243 <https://github.com/ESMValGroup/ESMValTool/issues/3243>`_
   * - :ref:`recipe_pysplot.yml <recipes_psyplot_diag>`
     - `plot_map`
     - ``shapely 2`` issue
     - `#3483 <https://github.com/ESMValGroup/ESMValTool/issues/3483>`_
