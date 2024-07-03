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
   * - :ref:`recipe_julia.yml <recipe_examples>`
     - `example`
     - v2.5.0
     - Fill values are not interpreted, resulting in an unusable plot
     - `#2595 <https://github.com/ESMValGroup/ESMValTool/issues/2595>`_
   * - :ref:`recipe_climwip_brunner2019_med.yml <recipe_climwip>`
     - All (preprocessor issue)
     - v2.11.0
     - Failed to run preprocessor function ``fix_metadata`` on the data: Unable to convert units
     - `#3694 <https://github.com/ESMValGroup/ESMValTool/issues/3694>`_
   * - :ref:`recipe_ocean_amoc.yml <XML_oceans>`
     - ``diag_timeseries_amoc``, ``diag_transects``
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3693>`_
   * - :ref:`recipe_preprocessor_derive_test.yml <recipe_examples>`
     - ``cmip6/toz``
     - v2.11.0
     - Failed to run preprocessor function ``derive`` on the data: Unable to convert units
     - `#3709 <https://github.com/ESMValGroup/ESMValTool/issues/3709>`_
   * - :ref:`recipe_russell18jgr.yml <nml_oceanmetrics>`
     - ``Figure_4``
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3693>`_
   * - :ref:`recipe_wenzel14jgr.yml <recipe_wenzel14jgr>`
     - ``diag_tsline_Fig2d``
     - v2.11.0
     - CESM1 CMIP5 Omon data no longer available
     - `#3693 <https://github.com/ESMValGroup/ESMValTool/issues/3693>`_
