.. _broken-recipe-list:

Broken recipe list
==================

This table gives an overview of the recipes that are known to have issues; the table is always updated for the
latest stable release of ESMValTool.
More details can be found in the :ref:`broken recipe policy
<broken-recipe-policy>`.

.. list-table:: Broken recipes
   :widths: 25 25 25 25
   :header-rows: 1

   * - Broken recipe
     - Affected diagnostics
     - Problem
     - GitHub issue 
   * - :ref:`recipe_autoassess_landsurface_soilmoisture.yml <recipe_autoassess_landsurface_soilmoisture.rst>`
     - All
     - Dependency on some external climatology files
     - `#2309 <https://github.com/ESMValGroup/ESMValTool/issues/2309>`_
   * - `recipe_check_obs.yml`
     - `ERA5_native6`
     - Derivation of custom variables `rlus` and `rsus`
     - `#1388 <https://github.com/ESMValGroup/ESMValCore/issues/1388>`_
