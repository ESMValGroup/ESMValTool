.. _broken-recipe-list:

Broken recipe list
==================

This table gives an overview of the recipes that are known to have issues.
The table is always valid for the latest stable release of ESMValTool.
More details can be found in the :ref:`broken recipe policy
<broken-recipe-policy>`.

.. list-table:: Broken recipes
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Broken recipe
     - Affected diagnostics
     - Broken since release
     - Problem
     - GitHub issue
   * - :ref:`recipe_climwip_brunner2019_med.yml`
     - All (preprocessor issue)
     - v2.11.0
     - Failed to run preprocessor function 'fix_metadata' on the data: Unable to convert units
     - `#3694 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_collins13ipcc.yml`
     - All (preprocessor issue)
     - v2.10.0
     - Failed to run preprocessor function 'save' on the data: HDF error
     - `#3702 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_easy_ipcc.yml`
     - All
     - v2.11.0
     - Failed to download data
     - `#3703 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_ipccwg1ar6ch3_atmosphere.yml`
     - All (preprocessor issue)
     - v2.10.0
     - Failed to run preprocessor function 'save' on the data: HDF error
     - `#3702 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_ocean_amoc.yml`
     - `diag_timeseries_amoc`, `diag_transects`
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_preprocessor_derive_test.yml`
     - All (preprocessor issue)
     - v2.11.0
     - Failed to run preprocessor function 'save' on the data: HDF error
     - `#3702 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_tebaldi21esd.yml`
     - All (preprocessor issue)
     - v2.10.0
     - Failed to run preprocessor function 'save' on the data: HDF error
     - `#3702 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_russell18jgr.yml`
     - `Figure_4`
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_wenzel14jgr.yml`
     - `diag_tsline_Fig2d`
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
