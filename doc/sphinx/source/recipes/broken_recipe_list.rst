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
   * - :ref:`recipe_carvalhais14nat.yml <recipe_carvalhais14nat>`
     - `global_turnover_time`
     - ``cartopy`` issue
     - `#3281 <https://github.com/ESMValGroup/ESMValTool/issues/3281>`_
   * - `recipe_check_obs.yml`
     - `ERA5_native6`
     - Derivation of custom variables `rlus` and `rsus`
     - `#1388 <https://github.com/ESMValGroup/ESMValCore/issues/1388>`_
   * - :ref:`recipe_seaice_drift.yml <recipes_seaice_drift>`
     - `sea_ice_drift_SCICEX`
     - ``shapely`` issue
     - `#3243 <https://github.com/ESMValGroup/ESMValTool/issues/3243>`_
