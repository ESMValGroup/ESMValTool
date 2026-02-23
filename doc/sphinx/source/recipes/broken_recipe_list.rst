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
   * - ``recipe_check_obs.yml``
     - All
     - v2.12.0
     - Various missing datasets
     - :issue:`3939`
   * - ``recipe_kcs.yml``
     - ``kcs/local_resampling.py``
     - v2.14.0
     - Diagnostic error (related to datetimes)
     - :issue:`4353`
