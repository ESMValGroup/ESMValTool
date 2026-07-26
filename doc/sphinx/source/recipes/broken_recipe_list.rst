:html_theme.sidebar_secondary.remove:

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
   * - ``recipe_kcs.yml``
     - ``kcs/local_resampling.py``
     - v2.14.0
     - Diagnostic error (related to datetimes)
     - :issue:`4353`
   * - ``recipe_weathertyping_CMIP6.yml``
     -
     - v2.14.0
     - Missing data (Observations, daily)
     - :issue:`4533`
   * - ``recipe_aod_aeronet_assess.yml``
     -
     - v2.15.0
     - Missing data (version for AERONET)
     - :issue:`4541`
   * - ``recipe_miles_regimes.yml``, ``recipe_miles_eof.yml``, ``recipe_miles_block.yml``
     - ``miles/miles_regimes.R``, ``miles/miles_eof.R``, ``miles/miles_block.R``
     - v2.15.0
     - Diagnostic error (hangs)
     - :issue:`4542`
   * - ``recipe_flato13ipcc_figures926_927.yml``
     - ``carbon_cycle/main.ncl``
     - v2.15.0
     - Diagnostic error (``fig09-27bottom`` hangs)
     - :issue:`4543`
   * - ``recipe_ipccwg1ar6ch3_atmosphere.yml``
     - ``ipcc_ar6/precip_anom.ncl``
     - v2.15.0
     - Diagnostic error (with MulitModelMean)
     - :issue:`4544`
   * - ``recipe_extreme_events.yml``
     - ``extreme_events/ extreme_events.R``
     - v2.15.0
     - Diagnostic error (open filename)
     - :issue:`4546`
   * - ``recipe_climate_patterns.yml``
     -
     - v2.15.0
     - Data/Fix error
     - :issue:`4558`
